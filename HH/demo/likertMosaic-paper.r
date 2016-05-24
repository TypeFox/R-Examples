## This file requires the HH package, version 2.3-32 or newer.
## Many of the examples will require stretching of the display window.
## Recommended width by height values are included in the comments.

## If you want to create pdf files, then manually set
##     PrintPDF <- TRUE
## before running this script.
## If the variable PrintPDF is not already set, it will be set to FALSE
if (class(try(PrintPDF, silent=TRUE))=="try-error")
  PrintPDF <- FALSE

library(HH)
library(vcd)

data(ProfChal)

PCQ <- as.character(ProfChal$Question)
PCQ[6] <- "Other (including retired, students,\nnot employed, etc.)" ## insert line break
ProfChal$Question <- factor(PCQ, levels=PCQ)
attributes(ProfChal)$names.dimnames <- c("Characteristic","ResponseLevel")

## Figure 1

if (PrintPDF) {
pdf("ProfChalMosaicPct.pdf", width=7, height=7)
}

ProfChalMosaicPctPlot <- ## 7in x 7in  ## ProfChalMosaicPct.pdf
likertMosaic(Question ~ . | Subtable, ProfChal,
       as.percent=TRUE,
       main="Is your job professionally challenging?",
       positive.order=TRUE
       )
ProfChalMosaicPctPlot

if (PrintPDF) {
dev.off()
}


## Table 1
EmpRows <- ProfChal$Subtable == "Employment sector"
ProfChal2 <- ProfChal[EmpRows, 1:5]
rownames(ProfChal2) <- substr(ProfChal[EmpRows, "Question"], 1, 5)
ProfChal2
as.likert(ProfChal2, padding=TRUE, reverse.left=FALSE)


## Table 2
rownames(ProfChal2) <- ProfChal$Question[EmpRows]
ProfChal2[1:5]


## Figure 2
if (PrintPDF) {
pdf("PC2CMosaic.pdf", width=9, height=3.5)
}
PC2CMosaic <-
likertMosaic(Question ~ . , data=ProfChal[EmpRows,],  ## 9in x 3in  ## PC2CMosaic.pdf
       main="Is your job professionally challenging?")
PC2CMosaic

if (PrintPDF) {
dev.off()
}



## Figure 3
if (PrintPDF) {
pdf("PC2Cpct.pdf", width=9, height=3.5)
}
PC2Cpct <-
likertMosaic(Question ~ . , data=ProfChal[EmpRows,],  ## 9in x 3in  ## PC2Cpct.pdf
       as.percent=TRUE,
       main="Is your job professionally challenging?")
PC2Cpct

if (PrintPDF) {
dev.off()
}


## Figure 4
if (PrintPDF) {
pdf("PC2Cpctpo.pdf", width=9, height=3)
}
PC2Cpctpo <-
likertMosaic(Question ~ . , data=ProfChal[EmpRows,],  ## 9in x 3in  ## PC2Cpctpo.pdf
       as.percent=TRUE,
       main="Is your job professionally challenging?",
       positive.order=TRUE)
PC2Cpctpo

if (PrintPDF) {
dev.off()
}




## Figure 5
## ProfChal6_131.pdf
## 9in x 9in
## ProfChal6_132.pdf
## 9in x 9in
if (PrintPDF) {
pdf("ProfChal6_AA.pdf", width=9, height=2.5)
}
likertMosaic(Question ~ . , ProfChal[1,],
       main=levels(ProfChal$Subtable)[1],
       as.percent=TRUE,
       positive.order=TRUE)
if (PrintPDF) {
dev.off()
}
if (PrintPDF) {
pdf("ProfChal6_BB.pdf", width=9, height=5)
}
likertMosaic(Question ~ . , ProfChal[2:6,],
       main=levels(ProfChal$Subtable)[2],
       as.percent=TRUE,
       positive.order=TRUE)
if (PrintPDF) {
dev.off()
}
if (PrintPDF) {
pdf("ProfChal6_CC.pdf", width=9, height=5)
}
likertMosaic(Question ~ . , ProfChal[7:10,],
       main=levels(ProfChal$Subtable)[3],
       as.percent=TRUE,
       positive.order=TRUE)
if (PrintPDF) {
dev.off()
}
if (PrintPDF) {
pdf("ProfChal6_DD.pdf", width=9, height=3)
}
likertMosaic(Question ~ . , ProfChal[11:12,],
       main=levels(ProfChal$Subtable)[4],
       as.percent=TRUE,
       positive.order=TRUE)
if (PrintPDF) {
dev.off()
}
if (PrintPDF) {
pdf("ProfChal6_EE.pdf", width=9, height=3)
}
likertMosaic(Question ~ . , ProfChal[13:14,],
       main=levels(ProfChal$Subtable)[5],
       as.percent=TRUE,
       positive.order=TRUE)
if (PrintPDF) {
dev.off()
}
if (PrintPDF) {
pdf("ProfChal6_FF.pdf", width=9, height=3)
}
likertMosaic(Question ~ . , ProfChal[15:16,],
       main=levels(ProfChal$Subtable)[6],
       as.percent=TRUE,
       positive.order=TRUE)
if (PrintPDF) {
dev.off()
}

## Figure 6
## Construct the two-column plot by printing the two one-column plots
##    ProfChalMosaicPct.pdf above
##    ProfChalMosaicCountPlot here
## to pdf and then place them adjacent to each other in the document.
if (PrintPDF) {
pdf("ProfChalMosaicCount.pdf", width=7, height=7)
}

ProfChalC <- ProfChal[match(rownames(ProfChalMosaicPctPlot), levels(ProfChal$Question)), ]
ProfChalC$Question <- factor(ProfChalC$Question, levels=rownames(ProfChalMosaicPctPlot))

ProfChalMosaicCountPlot <- ## 7in x 7in  ## ProfChalMosaicCount.pdf
likertMosaic(Question ~ . | Subtable, ProfChalC,
       main="Is your job professionally challenging?",
       )
ProfChalMosaicCountPlot

if (PrintPDF) {
dev.off()
}


## Figure 7
## http://www.morst.govt.nz/Documents/publications/researchreports/
## Staying-in-Science-summary.pdf

## New Zealand students still taking science in Year 13.
## "Students' interest in science, and views about secondary science teaching"
## "Students' feelings about making tertiary study decisions"
data(NZScienceTeaching)

if (PrintPDF) {
pdf("NZscienceteaching.pdf", width=9, height=5)
}
NZscienceteaching.structable <-
  likertMosaic(Question ~ . | Subtable, NZScienceTeaching, ## 10in x 5in  ## NZscienceteaching.pdf
               main="New Zealand Students Still Taking Science in Year 13",
               margins=c(3,2,4,24)
               )
NZscienceteaching.structable

if (PrintPDF) {
dev.off()
}



## Figure 8
## 9in x 7.5in   ## SFF8121.pdf
if (PrintPDF) {
pdf("SFF8121likert.pdf", width=9, height=10)
}
data(SFF8121)
SFF8121.likert <-
likertMosaic(aperm(SFF8121, c(3,1,2)),
             main="Student Feedback Forms, Spring 2010")
SFF8121.likert
if (PrintPDF) {
dev.off()
}
## > names(dimnames(aperm(SFF8121, c(3,1,2))))
## [1] "SDCLU"         "Question"      "ResponseLevel"


## Figure 9
## PopUSA1939-1979.pdf
## this needs a wx4 window
data(USAge.table) ## from package:latticeExtra

if (PrintPDF) {
pdf("PL5.pdf", width=7, height=4.5)
}
PL5 <-
  likertMosaic(aperm(USAge.table[75:1, 2:1, seq(40,80,10)], c(3,1,2))/1000000,
               main="Population of United States (ages 0-74)",
               sub="Look for the Baby Boom",
               margins=c(3,2,4,4),
               legend.y=.15
               )
PL5
if (PrintPDF) {
dev.off()
}



if (PrintPDF) {
pdf("PL.pdf", width=5.8, height=4.5)
}
PL <- likertMosaic( ~ . , data.frame(USAge.table[75:1, 2:1, "1979"]/1000000),
                   main="Population of United States 1979 (ages 0-74)",
                   margins=c(3,2,4,4),
                   legend.y=.15
                )
PL

if (PrintPDF) {
dev.off()
}


## Figure 10
## popUSA1979.pdf
## not currently available in Mosaic setting


## Figure 11
## top is historical
## bottom
## ProfitDividend.pdf
## 4.5in x 2.75in
data(ProfDiv)

if (PrintPDF) {
pdf("ProfDiv.pdf", width=8.8, height=4.5)
}
ProfDivPlot <-
  likertMosaic( ~ ., data.frame(ProfDiv),
               main=paste("\n","Profit-and-Dividend Status of 348 Corporations",
                 "in the United States\nfor the period from 1929 to 1935."),
               sub="Dun's Review, April 1938",
               margins=c(3,2,2,4),
               legend.y=.10
      )
ProfDivPlot
if (PrintPDF) {
dev.off()
}


## Figure 12
## PoorChildren
## Colors taken from NY Times figure using Snagit editor
PCWPpalette <- c("#A8C2C0", "#C9E9E6", "#D0E39A", "#A3B37B")
data(PoorChildren)

if (PrintPDF) {
pdf("PoorChildrenPlot.pdf", width=7.5, height=5)
}
PoorChildrenPlot <-
  likertMosaic( ~ . , PoorChildren,
               col=PCWPpalette, as.percent=TRUE,
               main="Poor Children, Working Parents",
               margins=c(1,2,3,6)
       )
PoorChildrenPlot
if (PrintPDF) {
dev.off()
}

if (PrintPDF) {
pdf("PoorChildrenPlotVariableWidth.pdf", width=7.5, height=5)
}
PoorChildrenPlotVariableWidth <-
  likertMosaic( ~ . , PoorChildren,
               col=PCWPpalette, as.percent=TRUE,
               variable.width=TRUE,
               main="Poor Children, Working Parents",
               margins=c(1,2,3,6)
               )
PoorChildrenPlotVariableWidth
if (PrintPDF) {
dev.off()
}



## Figure 13
tmp4 <- colSums(PoorChildren[,c(2,1,3,4)])
ByWP <- rbind(NWP=tmp4, '1+WP'=tmp4)
ByWP[cbind(c(2,2,1,1), 1:4)] <- 0
ByWP <- data.frame(ByWP,
                   Number=rowSums(ByWP),
                   NumberLabel=format(rowSums(ByWP), big.mark=","),
                   WP=factor(rownames(ByWP), levels=rownames(ByWP)),
                   check.names=FALSE)
ByWP
## note the sequencing: Mod, Ext, Mod, Ext
## this, together with Ref=0 in the likertMosaic() call,
## gets entire graph on the positive side in the intended order.

if (PrintPDF) {
  pdf("PClikNWC.pdf", width=8, height=3.5)
}
PClikNWC <-
  likertMosaic(NumberLabel ~ . , data=ByWP[,-5], Ref=0,
               col=PCWPpalette[c(2,1,3,4)], as.percent="noRightAxis",
               main="Poor Children, Working Parents",
               margins=c(1,2,3,6),
               x.legend=list(text=list(names(ByWP[c(4,3,1,2)])), columns=4,
                 space="bottom", size=2, cex=.8, between=1,
                 rect=list(col=PCWPpalette[c(4,3,2,1)], border="white")))
PClikNWC
if (PrintPDF) {
  dev.off()
}


## Figure 13 with width control by controlling panel width
if (PrintPDF) {
pdf("PClikWC.pdf", width=8, height=3.5)
}
PClikWC <-
  likertMosaic(NumberLabel ~ . | WP, data=ByWP[,-5], Ref=0,
               col=PCWPpalette[c(2,1,3,4)], as.percent="noRightAxis",
               variable.width=TRUE,
               main="Poor Children, Working Parents",
               margins=c(1,2,3,6),
               x.legend=list(text=list(names(ByWP[c(4,3,1,2)])), columns=4,
                 space="bottom", size=2, cex=.8, between=1,
                 rect=list(col=PCWPpalette[c(4,3,2,1)], border="white")))
PClikWC
if (PrintPDF) {
dev.off()
}




## Figure 14
## AudiencePercent.pdf
## 7in x 4in
data(AudiencePercent)

if (PrintPDF) {
pdf("AudiencePercentPlot.pdf", width=8, height=3)
}
AudiencePercentPlot <-
  likertMosaic( ~ . , data.frame(AudiencePercent),
               positive.order=TRUE,
               sub=paste("Percentage of audience younger than 35",
                 "(left of zero) and older than 35 (right of zero)"),
               main="Brand A has the most even distribution of ages",
               col=likertColor(nc=5, colorFunction="sequential_hcl")[2:5],
               margins=c(1,2,3,2),
               legend.y=.20)
AudiencePercentPlot
if (PrintPDF) {
dev.off()
}

## Grouped Bar Chart  ## we are not recommending this plot for this data
## Figure 15
## We do not construct it using Mosaic


## Heat Map   ## we are not recommending this plot for this data
## Figure 16 ## 9in x 5in
## PC2CpctpoHM.pdf
if (PrintPDF) {
pdf("PC2CpctpoHM.pdf", width=9, height=3.5)
}
PC2CpctpoHM <-
likertMosaic(Question ~ . , data=ProfChal[EmpRows,], ## 8in x 3in  ## PC2CpctpoHM.pdf
       as.percent=TRUE,
       main="Is your job professionally challenging?",
       ReferenceZero=0,
       col=likertColor(5),  ## col= needed because ReferenceZero=0 sets all colors to positive
       box.ratio=20)
PC2CpctpoHM
if (PrintPDF) {
dev.off()
}


## Figure 17 in the paper is
## Figure 1 from this file.


## Figure 18
## screen shot from Tableau worksheet JSS_figures6-1.twbx
## TableauWorksheet.pdf

## Figure 19
## pdf file drawn by Tableau
## tableau_redone2.pdf
