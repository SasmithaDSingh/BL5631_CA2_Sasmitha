#Install the required packages
BiocManager::install('IRanges')
BiocManager::install('GenomicFeatures')
BiocManager::install('GRanges')

#Loading packages into R
library(GenomicRanges)
library(IRanges)
library(rtracklayer)
library(AnnotationHub)
library(oligo)
library(ggplot2)

#[method1] Importing the human transcripts for GRCh38 
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
head(seqlevels(txdb),25)
transcripts <- transcripts(txdb)
transcripts

#[method2]Importing the human transcripts for GRCh38 
ahub <- AnnotationHub()
query(ahub, c("GRanges","Homo sapiens","GRCh38"))
GRCh38 <- ahub[['AH105317']]
seqlevels(GRCh38)
txdb_GRCh38 <- makeTxDbFromGRanges(GRCh38)
transcripts_GRCh38 <- transcripts(txdb_GRCh38, use.names=TRUE)
transcripts_GRCh38

#Creating a Granges object for the promoters of all protein-coding transcripts, defined for this problem set as 1500 bp upstream of the TSS and 500 bp downstream of the TSS.
#To get protein coding transcripts the coding sequences are got using the cds command
Protein_coding <- cds(txdb)
Protein_coding
Promoters <- promoters(Protein_coding, upstream = 1500, downstream = 500)
Promoters <- GRanges(Promoters)

#Creating a Granges object of all CpG islands for the human genome (had to import chain for conversion from hg19 to hg38)
ahub <- AnnotationHub()
query(ahub, c("GRanges", "Homo sapiens", "CpG Island"))
CpG_GR<- ahub[['AH5086']]
download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", "hg19ToHg38.over.chain.gz")
system("gzip -d hg19ToHg38.over.chain.gz")
Chain <- import.chain("hg19ToHg38.over.chain")
CpG_38<-liftOver(CpG_GR,Chain)
CpGislands_38 <- unlist(CpG_38)
genome(CpGislands_38) = "hg38"
CpGislands_38

#Calculating the fraction of CpG island annotations that overlap a promoter
Overlaps <- findOverlaps(CpGislands_38, Promoters)
Overlaps
Fraction <- length(unique(queryHits(Overlaps)))/length(CpGislands_38)
Fraction
Percentage_Overlap <- Fraction*100
Percentage_Overlap

##Plotting the length distributions for CpG islands that do and do not overlap a promoter
Overlapping <- subsetByOverlaps(CpGislands_38, Promoters)
Overlapping
R_OL <- as.data.frame(ranges(unique(Overlapping)))
#Length distribution density curve for overlapping
ggplot(R_OL, aes(x=log10(width))) + geom_density() + labs(title="Density curve for CpG islands overlapping promoters")
#Histogram overlaid with kernel density curve for overlapping
ggplot(R_OL, aes(x = log10(width)))+
  geom_histogram(aes(y = ..density..), colour = 1, fill = "white") + 
  geom_density(colour=22) + labs(title="Length/density distribution for CpG islands overlapping promoters")

Non_overlapping <- subsetByOverlaps(CpGislands_38, Promoters, invert = TRUE)
Non_overlapping
R_NOL <- as.data.frame(ranges(Non_overlapping))
#Length distribution density curve for non-overlapping
ggplot(R_NOL, aes(x=log10(width))) + geom_density() + labs(title="Density curve for CpG islands not overlapping promoters")
#Histogram overlaid with kernel density curve for non-overlapping
ggplot(R_NOL, aes(x = log10(width)))+
  geom_histogram(aes(y = ..density..), colour = 1, fill = "white") + 
  geom_density(colour=22) + labs(title="Length/density distribution for CpG islands not overlapping promoters")

#(optional) Plotting overlapping vs non-overlapping ranges and vice versa
plotRanges <- function(x, xlim=x, main=deparse(substitute(x)),
                       col="black", sep=0.5, ...) {
  height <- 1
  if (is(xlim, "IntegerRanges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col, ...)
  title(main)
  axis(1)}
plotRanges(ranges(unique(Overlapping)),ranges(Non_overlapping))
plotRanges(ranges(Non_overlapping), ranges(unique(Overlapping)))
