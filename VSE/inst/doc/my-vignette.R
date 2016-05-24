## ----warning=FALSE, message=FALSE, tidy=TRUE-----------------------------
library("VSE")
bca.ld <- loadLd(file.path(system.file("extdata", 
                                       "ld_BCa_raggr.csv", 
                                       package="VSE")), 
                 type="raggr")
#bca.ld is a GRanges object
bca.ld

## ----warning=FALSE, message=FALSE, tidy=TRUE, fig.height=5, fig.cap="LD Sizes of BCa AVS"----
bca.avs <- makeAVS(bca.ld)
#Check the size of each LD block
avs.size <- avsSize(bca.avs)
head(avs.size)
library(ggplot2)
ggplot(avs.size, aes(x=tagID, y=Size)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_classic() + 
  theme(axis.text.y = element_text(size=5))

## ----warning=FALSE, message=FALSE, tidy=TRUE-----------------------------
bca.avs

## ----warning=FALSE, message=FALSE, tidy=TRUE-----------------------------
# Load MRVS
load(file.path(system.file("extdata", 
                           "bca.mrvs.200.Rda", 
                           package="VSE")))
class(bca.mrvs.200)

## ----warning=FALSE, message=FALSE, tidy=TRUE-----------------------------
#Downloading sample regions
sampleSheet_path <- loadSampleRegions()
#Loading sample sheet
samples <- read.csv(sampleSheet_path, header=TRUE)
samples

## ----warning=FALSE, message=FALSE, tidy=TRUE, fig.height=5, fig.width=5----
bca.intersect <- intersectMatrix(bca.avs, 
                                 regions=samples, 
                                 col=c("white","grey10"), 
                                 scale="none", 
                                 margins=c(5,5), 
                                 cexRow = 1, 
                                 cexCol = 0.5, 
                                 Rowv=NA, 
                                 Colv=NA)
bca.intersect

## ----warning=FALSE, message=FALSE, tidy=TRUE-----------------------------
bca.vse <- variantSetEnrichment(bca.avs, bca.mrvs.200, samples) 

## ----warning=FALSE, message=FALSE, tidy=TRUE, fig.height=7, fig.width=7, fig.cap="QQ plot of null distribution across all genomic features"----
par.original <- par(no.readonly=TRUE);
par(mfrow=c(ceiling(length(samples$Peaks)/3), 3), 
    mai=c(1,1,0.5,0.1))
VSEqq(bca.vse)
par(par.original)

## ----warning=FALSE, message=FALSE, tidy=TRUE-----------------------------
bca.vse.res <- VSESummary(bca.vse)
bca.vse.res

## ----warning=FALSE, message=FALSE, tidy=TRUE, fig.cap="VSE analysis of BCa AVS across 6 genomic features. The points in red denotes significantly enriched region (Bonferroni adjusted P-value < 0.01)"----
VSEplot(bca.vse, 
        las=2, 
        pch=20, 
        cex=1, 
        cex.main=0.6, 
        padj = 0.05, 
        main="BCa AVS in MCF7 genomic features")

