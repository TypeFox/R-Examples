## This file requires the HH package, version 2.3-30 or newer.
## Many of the examples will require stretching of the display window.
## Recommended width by height values are included in the comments.

## If you want to create pdf files, then manually set
##     PrintPDF <- TRUE
## before running this script.
## If the variable PrintPDF is not already set, it will be set to FALSE
if (class(try(PrintPDF, silent=TRUE))=="try-error")
  PrintPDF <- FALSE

library(HH)

data(ProfChal)

PCQ <- as.character(ProfChal$Question)
PCQ[6] <- "Other (including retired, students,\nnot employed, etc.)" ## insert line break
ProfChal$Question <- factor(PCQ, levels=PCQ)
attributes(ProfChal)$names.dimnames <- c("Characteristic","ResponseLevel")

## Figure 1
ProfChalPctPlot <- ## 8in x 9in  ## ProfChal.pdf
likert(Question ~ . | Subtable, ProfChal,
       as.percent=TRUE,    ## implies display Row Count Totals
       ylab=NULL,
       main=list("Is your job professionally challenging?", x=unit(.65, "npc")),
       strip.left=strip.custom(bg="gray85"),
       strip=FALSE,
       par.strip.text=list(cex=.6, lines=5),
       positive.order=TRUE,
       layout=c(1,6),
       scales=list(y=list(relation="free"))  ## implies resizePanels
       )
ProfChalPctPlot

if (PrintPDF) {
pdf("ProfChal.pdf", width=8, height=9)
print(ProfChalPctPlot)
dev.off()
}


## Table 1
EmpRows <- ProfChal$Subtable == "Employment sector"
ProfChal2 <- ProfChal[EmpRows, 1:5]
rownames(ProfChal2) <- substr(ProfChal[EmpRows, "Question"], 1, 5)
ProfChal2
as.likert(ProfChal2)


## Table 2
rownames(ProfChal2) <- ProfChal$Question[EmpRows]
ProfChal2[1:5]


## Figure 2
PC2C <-
likert(Question ~ . , data=ProfChal[EmpRows,],  ## 9in x 3in  ## PC2C.pdf
       ylab=NULL,
       main="Is your job professionally challenging?")
PC2C

if (PrintPDF) {
pdf("PC2C.pdf", width=9, height=3)
print(PC2C)
dev.off()
}



## Figure 3
PC2Cpct <-
likert(Question ~ . , data=ProfChal[EmpRows,],  ## 9in x 3in  ## PC2Cpct.pdf
       as.percent=TRUE,
       ylab=NULL,
       main="Is your job professionally challenging?")
PC2Cpct

if (PrintPDF) {
pdf("PC2Cpct.pdf", width=9, height=3)
print(PC2Cpct)
dev.off()
}


## Figure 4
PC2Cpctpo <-
likert(Question ~ . , data=ProfChal[EmpRows,],  ## 9in x 3in  ## PC2Cpctpo.pdf
       as.percent=TRUE,
       ylab=NULL,
       main="Is your job professionally challenging?",
       positive.order=TRUE)
PC2Cpctpo

if (PrintPDF) {
pdf("PC2Cpctpo.pdf", width=9, height=3)
print(PC2Cpctpo)
dev.off()
}




## Figure 5
## ProfChal6_131.pdf
## 9in x 9in
## ProfChal6_132.pdf
## 9in x 9in
data(ProfChal)
##
AA <- likert(Question ~ . , ProfChal[1,],
             ylab=NULL,
             main=levels(ProfChal$Subtable)[1],
             as.percent=TRUE, box.width=unit(.4,"cm"),
             positive.order=TRUE)
BB <- likert(Question ~ . , ProfChal[2:6,],
             ylab=NULL,
             main=levels(ProfChal$Subtable)[2],
             as.percent=TRUE, box.width=unit(.4,"cm"),
             positive.order=TRUE)
CC <- likert(Question ~ . , ProfChal[7:10,],
             ylab=NULL,
             main=levels(ProfChal$Subtable)[3],
             as.percent=TRUE, box.width=unit(.4,"cm"),
             positive.order=TRUE)
DD <- likert(Question ~ . , ProfChal[11:12,],
             ylab=NULL,
             main=levels(ProfChal$Subtable)[4],
             as.percent=TRUE, box.width=unit(.4,"cm"),
             positive.order=TRUE)
EE <- likert(Question ~ . , ProfChal[13:14,],
             ylab=NULL,
             main=levels(ProfChal$Subtable)[5],
             as.percent=TRUE, box.width=unit(.4,"cm"),
             positive.order=TRUE)
FF <- likert(Question ~ . , ProfChal[15:16,],
             ylab=NULL,
             main=levels(ProfChal$Subtable)[6],
             as.percent=TRUE, box.width=unit(.4,"cm"),
             positive.order=TRUE)
##
print(AA, more=TRUE,  split=c(1,1,1,3))
print(BB, more=TRUE,  split=c(1,2,1,3))
print(CC, more=FALSE,  split=c(1,3,1,3))
## ProfChal6_131.pdf
print(DD, more=TRUE,  split=c(1,1,1,3))
print(EE, more=TRUE,  split=c(1,2,1,3))
print(FF, more=FALSE, split=c(1,3,1,3))
## ProfChal6_132.pdf

if (PrintPDF) {
pdf("ProfChal6_131.pdf", width=9, height=9)
print(AA, more=TRUE,  split=c(1,1,1,3))
print(BB, more=TRUE,  split=c(1,2,1,3))
print(CC, more=FALSE,  split=c(1,3,1,3))
dev.off()
pdf("ProfChal6_132.pdf", width=9, height=9)
print(DD, more=TRUE,  split=c(1,1,1,3))
print(EE, more=TRUE,  split=c(1,2,1,3))
print(FF, more=FALSE, split=c(1,3,1,3))
dev.off()
}

## Figure 1, reprise
ProfChalReprise <- resizePanels(c(AA,BB,CC,DD,EE,FF, layout=c(1,6), x.same=TRUE), h=c(1,5,4,2,2,2))
ProfChalReprise$condlevels <- list(ListNames=levels(ProfChal$Subtable))
ProfChalReprise <- update(ProfChalReprise, main=ProfChalPctPlot$main, strip.left=strip.custom(bg="gray97"),
                          par.strip.text=list(cex=.6, lines=5))
ProfChalReprise$y.limits <- lapply(list(AA,BB,CC,DD,EE,FF), `[[`, "y.limits")
ProfChalReprise

if (PrintPDF) {
pdf("ProfChalReprise.pdf", width=9, height=7)
print(ProfChalReprise)
dev.off()
}

## Figure 6
## ProfChalPctPlot ## from above

ProfChalCountPlot <- ## 8in x 9in  ## ProfChal.pdf
  likert(Question ~ . | Subtable, ProfChal,
         ylab=NULL,
         rightAxis=TRUE,
         main=list("Is your job professionally challenging?", x=unit(.65, "npc")),
         strip.left=strip.custom(bg="gray85"),
         strip=FALSE,
         par.strip.text=list(cex=.6, lines=5),
         positive.order=TRUE,
         layout=c(1,6),
         scales=list(y=list(relation="free"))  ## implies resizePanels
         )
ProfChalCountPlot
##
ProfChalPercentCount <- ## 11in x 9in  ## twocolumn.pdf
  as.TwoTrellisColumns5(update(ProfChalPctPlot,
                               rightAxis=FALSE,
                               sub="This plot needs an 11in x 9in plotting surface"),
                        update(ProfChalCountPlot, sub=" "),
                        pw=c(.28, .39, .01, .22, .10))
ProfChalPercentCount

if (PrintPDF) {
pdf("twocolumn.pdf", width=11, height=9)
print(ProfChalPercentCount)
dev.off()
}


## Figure 7
## http://www.morst.govt.nz/Documents/publications/researchreports/
## Staying-in-Science-summary.pdf

## New Zealand students still taking science in Year 13.
## "Students' interest in science, and views about secondary science teaching"
## "Students' feelings about making tertiary study decisions"
data(NZScienceTeaching)

NZscienceteaching <-
likert(Question ~ . | Subtable, NZScienceTeaching, ## 10in x 5in  ## NZscienceteaching.pdf
       main="New Zealand Students Still Taking Science in Year 13",
       layout=c(1,2),
       xlab="Percent",
       ylab=NULL,
       scales=list(y=list(relation="free")),
       strip=FALSE,
       strip.left=strip.custom(bg="gray97"),
       par.strip.text=list(cex=1.2, lines=1.2)
       )
NZscienceteaching

if (PrintPDF) {
pdf("NZscienceteaching.pdf", width=10, height=5)
print(NZscienceteaching)
dev.off()
}



## Figure 8
## 9in x 7.5in   ## SFF8121.pdf
data(SFF8121)
SFF8121.df <- as.likertDataFrame(SFF8121)
## as.MatrixList(SFF8121)
SFF8121.likert <-
likert(Question ~ . | Subtable, SFF8121.df,
       layout=c(2,1),
       xlab="Percent",
       ylab="",
       main="Student Feedback Forms, Spring 2010",
       strip=strip.custom(bg="gray97"),
       par.strip.text=list(cex=.6))
SFF8121.likert
## warnings -- wrong bg in pdf in Mac 15.2.2

if (PrintPDF) {
pdf("SFF8121.pdf", width=9, height=7.5)
print(SFF8121.likert)
dev.off()
}


## Figure 9
## PopUSA1939-1979.pdf
## this needs a 17x8 window
data(USAge.table) ## from package:latticeExtra
tmp.df <- as.likertDataFrame(USAge.table[1:75, 2:1, seq(40,80,10)]/1000000)
names(tmp.df)[3] <- "Age"
names(tmp.df)[4] <- "Year"
PL5 <-
likert(Age ~ . | Year, tmp.df,
       main="Population of United States (ages 0-74)",
       xlab="Count in Millions",
       ylab="Age",
       sub="Look for the Baby Boom",
       reverse=FALSE,
       scales=list(
         y=list(
           limits=c(0,77),
           at=seq(1,76,5),
           labels=seq(0,75,5),
           alternating=3,
           tck=1),
         x=list(alternating=FALSE, at=-2:2)),
       auto.key=list(title=NULL),
       strip=strip.custom(bg="gray97"),
       par.strip.text=list(cex=1),
       layout=c(5,1), between=list(x=.5),
       yscale.components=yscale.components.default)
PL5

if (PrintPDF) {
pdf("PL5.pdf", width=17, height=8)
print(PL5)
dev.off()
}


## 7x8 window
PL <- likert( ~ . , data.frame(USAge.table[1:75, 2:1, "1979"]/1000000),
             main="Population of United States 1979 (ages 0-74)",
             xlab="Count in Millions",
             ylab="Age",
             reverse=FALSE,
             scales=list(
               y=list(
                 limits=c(0,77),
                 at=seq(1,76,5),
                 labels=seq(0,75,5),
                 tck=.5))
                )
PL

if (PrintPDF) {
pdf("PL.pdf", width=7, height=8)
print(PL)
dev.off()
}


## Figure 10
## popUSA1979.pdf
as.pyramidLikert(PL)

if (PrintPDF) {
pdf("PLpyr.pdf", width=7, height=8)
print(as.pyramidLikert(PL))
dev.off()
}


## Figure 11
## top is historical
## bottom
## ProfitDividend.pdf
## 4.5in x 2.75in
data(ProfDiv)
ProfDivPlot <-
likert( ~ ., data.frame(ProfDiv), horizontal=FALSE, positive.order=FALSE,
       auto.key=list(size=4, padding.text=1.5, title=names(dimnames(ProfDiv))[2], cex.title=.9),
       ylab="Per Cent",
       xlab=names(dimnames(ProfDiv))[1],
       main=paste("Profit-and-Dividend Status of 348 Corporations",
         "in the United States\nfor the period from 1929 to 1935."),
       sub="Dun's Review, April 1938"
      )
ProfDivPlot

if (PrintPDF) {
pdf("ProfDiv.pdf", width=7, height=4)
print(ProfDivPlot)
dev.off()
}


## Figure 12
## PoorChildren
## Colors taken from NY Times figure using Snagit editor
PCWPpalette <- c("#A8C2C0", "#C9E9E6", "#D0E39A", "#A3B37B")
data(PoorChildren)
rowCounts <- rowSums(PoorChildren[,1:4])

PoorChildrenPlot <-
likert( ~ . | PercentPoorInArea, PoorChildren,
       col=PCWPpalette, as.percent=TRUE,
       ylab="Percent of poor households in area",
       xlab="Percent of Children",
       xlab.top=c("No Working Parents", "One or more Working Parents"),
       ylab.right="Row Count Totals",
       main="Poor Children, Working Parents",
       strip=FALSE,
       rightAxisLabels=format(rowCounts, big.mark=","),
       par.settings=list(
         axis.line=list(col="transparent"),
         layout.widths=list(ylab.right=1, axis.key.padding=10)
         ),
       box.ratio=500,
       layout=c(1,6),
       between=list(y=0),
       scales=list(y=list(relation="free", tck=c(0,3))),
       h.resizePanels=rowCounts
       )
PoorChildrenPlot


if (PrintPDF) {
pdf("PoorChildrenPlot.pdf", width=8, height=7)
print(PoorChildrenPlot)
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
## this, together with Ref=0 in the likert() call,
## gets entire graph on the positive side in the intended order.

## Figure 13, no box.width control because
##    box.width=rep(rowSums(ByWP[,1:4]), each=2)/8000000
## is applied by barchart to each item in the stacked bar, not to the entire bar.
##
PClikNWC <-
  likert(NumberLabel ~ ., data=ByWP[,-5], horizontal=FALSE, Ref=0,
         col=PCWPpalette[c(2,1,3,4)], as.percent="noRightAxis",
         xlab="Number of Children",
         ylab="Percent of Children, All Areas",
         xlab.top=c("No Working Parents", "One or More Working Parents"),
         ylab.right=list(c("Moderately\nPoor", "Extremely\nPoor"), rot=0),
         main="Poor Children, Working Parents",
         auto.key=list(
           space="bottom", columns=4, between=1,
           text=names(ByWP[c(4,3,1,2)]),
           rect=list(col=PCWPpalette[c(4,3,2,1)], border="white")))
PClikNWC

if (PrintPDF) {
pdf("PClikNWC.pdf", width=8, height=7)
print(PClikNWC)
dev.off()
}


## Figure 13 with width control by controlling panel width
PClikWC <-
  likert(NumberLabel ~ . | WP, data=ByWP[-5], horizontal=FALSE, Ref=0,
         col=PCWPpalette[c(2,1,3,4)], as.percent="noRightAxis",
         xlab="Number of Children",
         ylab="Percent of Children, All Areas",
         xlab.top=c("No Working Parents", "One or More Working Parents"),
         ylab.right=list(c("Moderately\nPoor", "Extremely\nPoor"), rot=0),
         main="Poor Children, Working Parents",
         auto.key=list(
           space="bottom", columns=4, between=1,
           text=names(ByWP[c(4,3,1,2)]),
           rect=list(col=PCWPpalette[c(4,3,2,1)], border="white")),
         ## arguments for multi-panel width control
         box.ratio=1000,
         strip=FALSE,
         layout=c(2,1),
         scales=list(x=list(relation="free")),
         par.settings=list(axis.line=list(col="transparent")),
         w.resizePanels=ByWP$Number ## this line makes the panel widths
         ## hence the effective box.width,
         ## proportional to the Number of persons
         )
PClikWC

if (PrintPDF) {
pdf("PClikWC.pdf", width=8, height=7)
print(PClikWC)
dev.off()
}



## Figure 14
## AudiencePercent.pdf
## 7in x 4in
data(AudiencePercent)
AudiencePercentPlot <-
likert( ~ . , data.frame(AudiencePercent, check.names=FALSE),
       positive.order=TRUE,
       auto.key=list(between=1, between.columns=2),
       xlab=paste("Percentage of audience younger than 35",
         "(left of zero) and older than 35 (right of zero)"),
       ylab="Brands",
       main="Brand A has the most even distribution of ages",
       scales=list(x=list(at=seq(-90,60,10),
                     labels=as.vector(rbind("", abs(seq(-80,60,20)))))),
       col=likertColor(nc=5, colorFunction="sequential_hcl")[2:5])
AudiencePercentPlot$x.scales$at <- seq(-90,60,10)   ## I don't know why these lines are needed
AudiencePercentPlot$x.scales$labels <- as.vector(rbind("", abs(seq(-80,60,20))))
AudiencePercentPlot

if (PrintPDF) {
pdf("AudiencePercentPlot.pdf", width=7, height=3.7)
print(AudiencePercentPlot)
dev.off()
}


## Grouped Bar Chart  ## we are not recommending this plot for this data
## Figure 15
tmp <- ProfChal2
dimnames(tmp)[[1]][5] <- "Other (including retired, students,\nnot employed, etc.)"
tmp <- data.frame(stack(tmp),
                  Employment=factor(row.names(tmp), levels=row.names(tmp)))
names(tmp)[1:2] <- c("Count", "Agreement")
tmp$Agreement <- factor(tmp$Agreement, unique(tmp$Agreement))
tmp$Percent <- unlist(ProfChal2 / rowSums(ProfChal2))

## Figure 15a
## PctAEv.pdf ## 16in x 6in
PctAEv <-
barchart(Percent ~ Agreement | Employment, data=tmp,
         col=likertColor(5),
         horizontal=FALSE, layout=c(5,1), between=list(x=1), box.ratio=10,
         origin=0,
         par.strip.text=list(cex=.7, lines=2.5),  ## cex too large on quartz
         scales=list(x=list(rot=90)),
         main="Percent ~ Agreement | Employment")
PctAEv

if (PrintPDF) {
pdf("PctAEv.pdf", width=11.5, height=4)
print(PctAEv)
dev.off()
}


## Figure 15b
## PctAEh.pdf ## 16in x 4in
PctAEh <-
barchart(Agreement ~ Percent | Employment, data=tmp,
         col=likertColor(5),
         horizontal=TRUE, layout=c(5,1), between=list(x=1), box.ratio=10,
         origin=0,
         par.strip.text=list(cex=.7, lines=2.5),  ## cex too large on quartz
         scales=list(x=list(alternating=1)),
         main="Agreement ~ Percent | Employment")
PctAEh

if (PrintPDF) {
pdf("PctAEh.pdf", width=12, height=3.5)
print(PctAEh)
dev.off()
}



## Heat Map   ## we are not recommending this plot for this data
## Figure 16 ## 9in x 5in
## PC2CpctpoHM.pdf
PC2CpctpoHM <-
likert( ~ . , data=ProfChal2, ## 8in x 3in  ## PC2CpctpoHM.pdf
       as.percent=TRUE,
       main="Is your job professionally challenging?",
       ReferenceZero=0,
       col=likertColor(5),  ## col= needed because ReferenceZero=0 sets all colors to positive
       ylab="",
       box.ratio=20)
PC2CpctpoHM

if (PrintPDF) {
pdf("PC2CpctpoHM.pdf", width=9, height=5)
print(PC2CpctpoHM)
dev.off()
}


## Figure 17
## Same content as Figure 1, as drawn in file likertMosaic-paper.r


## Figure 18
## screen shot from Tableau worksheet JSS_figures6-1.twbx
## TableauWorksheet.pdf

## Figure 19
## pdf file drawn by Tableau
## tableau_redone2.pdf
