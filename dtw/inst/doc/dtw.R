### R code from vignette source 'dtw.Rnw'

###################################################
### code chunk number 1: init-prompt
###################################################
options(prompt= "R> ", continue = "+  ")


###################################################
### code chunk number 2: plot-visual-explanation
###################################################
library("dtw")
data("aami3a")
ref <- window(aami3a,start=0,end=2)
test <- window(aami3a,start=2.7,end=5)
plot(dtw(test,ref,k=TRUE),type="two",off=1,match.lty=2,match.indices=20)


###################################################
### code chunk number 3: dtw-simple-example
###################################################
library("dtw")
data("aami3a")
ref <- window(aami3a,start=0,end=2)
test <- window(aami3a,start=2.7,end=5)
alignment <- dtw(test,ref)
alignment$distance


###################################################
### code chunk number 4: dtw-normalized-distance
###################################################
alignment <- dtw(test,ref,step.pattern=asymmetric)
alignment$distance


###################################################
### code chunk number 5: symmetric2-step-pattern
###################################################
symmetric2


###################################################
### code chunk number 6: plot-symmetric1-step-pattern
###################################################
plot(symmetric1)


###################################################
### code chunk number 7: plot-symmetric2-step-pattern
###################################################
plot(symmetric2)


###################################################
### code chunk number 8: plot-asymmetric-step-pattern
###################################################
plot(asymmetric)


###################################################
### code chunk number 9: plot-rj4c-step-pattern
###################################################
plot(rabinerJuangStepPattern(4,"c",TRUE))


###################################################
### code chunk number 10: plot-sakoechiba-window
###################################################
dtwWindow.plot(sakoeChibaWindow, window.size=2,reference=17, query=13)


###################################################
### code chunk number 11: plot-open-begin-end-dtw
###################################################

idx<-seq(0,6.28,len=100)
query<-sin(idx)+runif(100)/10
reference<-cos(idx)

alignmentOBE <-
  dtw(query[44:88],reference,
      keep=TRUE,step=asymmetric,
      open.end=TRUE,open.begin=TRUE)
plot(alignmentOBE,type="two",off=1)



###################################################
### code chunk number 12: dtw-multivariate-example
###################################################
query <- cbind(1:10,1)
ref <- cbind(11:15,2)
dtw(query,ref,dist.method="Manhattan")$distance


###################################################
### code chunk number 13: dtw-multivariate-example2
###################################################
cxdist <- proxy::dist(query,ref,method="Manhattan")
dtw(cxdist)$distance


###################################################
### code chunk number 14: rabiner-exercise-prepare
###################################################
lm <- matrix(nrow = 6, ncol = 6, byrow = TRUE, c(
  1, 1, 2, 2, 3, 3, 
  1, 1, 1, 2, 2, 2, 
  3, 1, 2, 2, 3, 3, 
  3, 1, 2, 1, 1, 2, 
  3, 2, 1, 2, 1, 2, 
  3, 3, 3, 2, 1, 2
))


###################################################
### code chunk number 15: rabiner-exercise-1
###################################################
alignment <- dtw(lm,step=asymmetric,keep=TRUE)
alignment$costMatrix
alignment$normalizedDistance


###################################################
### code chunk number 16: rabiner-exercise-2
###################################################
alignmentOE <- dtw(lm,step=asymmetric,keep=TRUE,open.end=TRUE)
alignmentOE$jmin
alignmentOE$normalizedDistance


###################################################
### code chunk number 17: costmatrix-figure-code-left (eval = FALSE)
###################################################
## lcm <- alignment$localCostMatrix
## image(x=1:nrow(lcm),y=1:ncol(lcm),lcm)
## text(row(lcm),col(lcm),label=lcm)
## lines(alignment$index1,alignment$index2)


###################################################
### code chunk number 18: costmatrix-figure-code-right (eval = FALSE)
###################################################
## ccm <- alignment$costMatrix
## image(x=1:nrow(ccm),y=1:ncol(ccm),ccm)
## text(row(ccm),col(ccm),label=ccm)
## lines(alignment$index1,alignment$index2)


###################################################
### code chunk number 19: costmatrix-figure-plot-left
###################################################
lcm <- alignment$localCostMatrix
image(x=1:nrow(lcm),y=1:ncol(lcm),lcm)
text(row(lcm),col(lcm),label=lcm)
lines(alignment$index1,alignment$index2)


###################################################
### code chunk number 20: costmatrix-figure-plot-right
###################################################
ccm <- alignment$costMatrix
image(x=1:nrow(ccm),y=1:ncol(ccm),ccm)
text(row(ccm),col(ccm),label=ccm)
lines(alignment$index1,alignment$index2)


###################################################
### code chunk number 21: two-way-plot-code
###################################################
library("dtw")
data("aami3a")
ref <- window(aami3a,start=0,end=2)
test <- window(aami3a,start=2.7,end=5)
plot(dtw(test,ref,k=TRUE),type="two",off=1,match.lty=2,match.indices=20)


###################################################
### code chunk number 22: three-way-plot
###################################################
## A noisy sine wave as query
## A cosine is for reference; sin and cos are offset by 25 samples
idx<-seq(0,6.28,len=100)
query<-sin(idx)+runif(100)/10
reference<-cos(idx)
dtw(query,reference,step=asymmetric,keep=TRUE)->alignment

## Beware of the reference's y axis, may be confusing
## Equivalent to plot(alignment,type="three")
#dtwPlotThreeWay(alignment)

## Highlight matches of chosen QUERY indices. We will do some index
## arithmetics to recover the corresponding indices along the warping
## curve
hq <- (0:8)/8
hq <- round(hq*100)      #  indices in query for  pi/4 .. 7/4 pi
hw <- (alignment$index1 %in% hq)   # where are they on the w. curve?
hi <- (1:length(alignment$index1))[hw]   # get the indices of TRUE elems
dtwPlotThreeWay(alignment,match.indices=hi)


###################################################
### code chunk number 23: density-plot
###################################################
dtw(query,reference,keep=TRUE,step=asymmetric)->alignment
dtwPlotDensity(alignment,normalize=TRUE)


