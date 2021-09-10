RNA-seq analysis
================
April
9/9/2021

## 1. retrieve samples & download files (on command line)

``` r
#download reads for this sample
prefetch SRR5293974
fastq-dump SRR5293974.sra

#download genome 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz
gunzip Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz

#download annotation file
wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz
gunzip Homo_sapiens.GRCh38.101.gtf.gz

#download GO gene list
#go to https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/7.2/c5.go.bp.v7.2.symbols.gmt to download, will save file called 'c5.go.bp.v7.2.symbols.gmt'

#create list that contains ensembl ID and gene symbol using gtf file
zcat Homo_sapiens.GRCh38.101.gtf | awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[5]"\t"$7}' | sed 's/gene_id "//' | sed 's/gene_id "//' | sed 's/gene_biotype "//'| sed 's/gene_name "//' | sed 's/gene_biotype "//' | sed 's/"//g' | sed 's/ //g' | sed '1igene_id\tGeneSymbol\tChromosome\tClass\tStrand' > Homo_sapiens.GRCh38.101_gene_annotation_table.txt
```

## 2. QC/trim reads (on command line)

``` r
trim_galore --fastqc --gzip --cores 7 --output_dir . SRR5293974.fastq
```

## 3. alignment (on command line)

``` r
#create index for hisat
hisat2-build Homo_sapiens.GRCh38.dna_sm.toplevel.fa hisat_idx
#align samples
hisat2 -x hisat_idx --rna-strandness RF -U SRR5293974_trimmed.fq.gz -S SRR5293974.sam
samtools view -bS SRR5293974.sam > SRR5293974.bam
samtools sort SRR5293974.bam -o SRR5293974.sorted.bam
samtools index SRR5293974.sorted.bam
```

# 4. quantification (on command line)

``` r
htseq-count -f bam --stranded=reverse SRR5293974.sorted.bam Homo_sapiens.GRCh38.101.gtf > SRR5293974_hisat_counts.txt
```

# 5 DE analysis (in R)

``` r
#load packages
library("DESeq2") #used for differential expression analysis
library(plyr) #used for data manipulation
library(dplyr) #used for data manipulation
library(fgsea) #to do geneset enrichment analysis 
library(tidyverse) #used for data manipulation
library(ggplot2) #used for plotting
```

``` r
sample_info<-read.delim('sample_info.txt',stringsAsFactors=FALSE,sep = ' ')
rownames(sample_info)<-sample_info$sample
sample_info$group<-factor(sample_info$group)
```

``` r
#read in samples that were aligned with hisat2 and quantified with htseq
files <- list.files('.', full.names=FALSE)
files<-files[grep(x=files,pattern = '.*_hisat_counts.txt')]
files<-grep(paste(sample_info$sample,collapse="|"), 
                        files, value=TRUE)

hisat<-read.delim(files[1],stringsAsFactors = FALSE,header = FALSE,
                  col.names = c('gene_id',gsub(x=files[1],pattern='_hisat_counts.txt',replacement = '')))

for(file in files[-1])
{
  temp<-read.delim(file,stringsAsFactors = FALSE,header = FALSE,
             col.names = c('gene_id',gsub(x=file,pattern='_hisat_counts.txt',replacement = '')))
  hisat<-join(hisat,temp)
}

rownames(hisat)<-hisat$gene_id
hisat$gene_id<-NULL
hisat<-hisat[-((nrow(hisat)-4):nrow(hisat)),] #remove last 4 lines b/c they are sample wide stats not genes
#filter out lowly expressed genes (keep only genes that have more than 5 reads in atleast 2 samples)
trans2keep<-which(rowSums(hisat > 5) >= 2)
hisat<-hisat[trans2keep,]
```

``` r
hisat_counts<-as.matrix(hisat)
all(rownames(sample_info) == colnames(hisat_counts))
dds <- DESeqDataSetFromMatrix(countData = hisat_counts,
                              colData = sample_info,
                              design = ~ group)
dds <- DESeq(dds)
```

``` r
resC_CDS <- results(dds,
                     contrast=c("group","Si_RNF219","SCR"))
#number of DE genes at pvalue adjusted <0.05 and -0.5 < logfold change > 0.5
nrow(filter(as.data.frame(resC_CDS), padj<0.05 & abs(log2FoldChange)>0.5))
```

    ## [1] 1222

``` r
gene_list<-read.delim('Homo_sapiens.GRCh38.101_gene_annotation_table.txt',stringsAsFactors = FALSE)
resC_CDS<-as.data.frame(resC_CDS)
resC_CDS$ensembl_id<-rownames(resC_CDS)
resC_CDS <- inner_join(resC_CDS, gene_list, by=c("ensembl_id"="gene_id"))
head(resC_CDS)
```

    ##     baseMean log2FoldChange     lfcSE       stat      pvalue       padj
    ## 1  135.58133    -1.19000725 0.4341063 -2.7412809 0.006120016 0.05627599
    ## 2  351.66045     0.20976040 0.3454823  0.6071524 0.543749774 0.80837323
    ## 3   53.28986    -0.10635681 0.5339038 -0.1992059 0.842101648 0.94840620
    ## 4  153.74288    -0.39509211 0.4456623 -0.8865280 0.375333095 0.69053832
    ## 5 2227.75110    -0.41851605 0.2967653 -1.4102596 0.158463042 0.45695961
    ## 6  136.28823     0.08361291 0.3975082  0.2103426 0.833400293 0.94446250
    ##        ensembl_id GeneSymbol            Chromosome          Class Strand
    ## 1 ENSG00000000003     TSPAN6 X:100627108-100639991 protein_coding      -
    ## 2 ENSG00000000419       DPM1  20:50934867-50958555 protein_coding      -
    ## 3 ENSG00000000457      SCYL3 1:169849631-169894267 protein_coding      -
    ## 4 ENSG00000000460   C1orf112 1:169662007-169854080 protein_coding      +
    ## 5 ENSG00000001036      FUCA2 6:143494812-143511720 protein_coding      -
    ## 6 ENSG00000001084       GCLC   6:53497341-53616970 protein_coding      -

``` r
resC_CDS_DE<-mutate(resC_CDS,
       DE= case_when(
         padj  <0.05 & log2FoldChange >0.5 ~ "up",
         padj  <0.05 & log2FoldChange < (-0.5) ~ "down",
         padj  >=0.05 | is.na(padj) ~ "not DE")
       )

ggplot(resC_CDS_DE) +
        geom_point(aes(x = log2FoldChange, y = -log10(padj), color = DE)) +
        xlab("log2 fold change") + 
        ylab("-log10 adjusted p-value") +
        theme_minimal()+
        theme(legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25)))  +
  scale_color_manual(values = c("#E31A1C", "black","#1F78B4"))
```

![](RNAseq_analysis_code_9_9_21_github_files/figure-markdown_github/unnamed-chunk-11-1.png)

# 5 GSEA analysis (in R)

``` r
#read in gene list
GO.BP.pathways<-gmtPathways("c5.go.bp.v7.2.symbols.gmt")

ranks_C_CDS <- resC_CDS %>% 
  dplyr::select(GeneSymbol, log2FoldChange) %>% 
  na.omit() %>% 
  distinct()%>%
  dplyr::group_by(GeneSymbol) %>% 
  dplyr::summarize(stat=mean(log2FoldChange))

ranks_C_CDS <- ranks_C_CDS %>% arrange(stat)
ranks_C_CDS <- deframe(ranks_C_CDS)

fgsea_C_CDS <- fgseaMultilevel(GO.BP.pathways,ranks_C_CDS, minSize = 10)
fgsea_C_CDS_sig<- filter(fgsea_C_CDS, padj<0.05)
fgsea_C_CDS_sig<-fgsea_C_CDS_sig[order(fgsea_C_CDS_sig$NES,-fgsea_C_CDS_sig$padj),]
```

``` r
#clean up pathway names for plot
fgsea_C_CDS_sig$pathway<-gsub(fgsea_C_CDS_sig$pathway, pattern='GO_', replacement = '')
fgsea_C_CDS_sig$pathway<-gsub(fgsea_C_CDS_sig$pathway, pattern='_', replacement = ' ')
fgsea_C_CDS_sig$pathway<-str_to_title(fgsea_C_CDS_sig$pathway)

#GSEA plot
ggplot(rbind(head(fgsea_C_CDS_sig,10),tail(fgsea_C_CDS_sig,10)), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES>0)) +
  coord_flip() +
  labs(x="GO BP pathways", y="Normalized Enrichment Score") + 
  theme_minimal()+
  scale_fill_manual(values=c("#E31A1C","#1F78B4"))+
  theme(legend.position = "none")
```

![](RNAseq_analysis_code_9_9_21_github_files/figure-markdown_github/unnamed-chunk-13-1.png)
