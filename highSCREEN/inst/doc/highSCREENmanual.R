## ----echo=TRUE-----------------------------------------------------------
set.seed(1234)
library(highSCREEN)
nc = 24
nr = 16
# create a 384-well plate with compounds and controls
replicate = matrix(abs(rnorm(nr*nc)), nr, nc)
head(replicate)

# create 384-well plate control map
cmap = data.frame(X1=c(rep("Control P", floor(nr/3)), rep(c("Control low", "Control med", "Control high"), 
(floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control N", floor(nr/3))), X2=c(rep("Control P", floor(nr/3)), 
rep(c("Control low", "Control med", "Control high"), (floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control N", 
floor(nr/3))), X3=c(rep("Control N", floor(nr/3)), rep(c("Control low", "Control med", "Control high"), 
(floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control P", floor(nr/3))), X4=c(rep("Control N", floor(nr/3)), 
rep(c("Control low", "Control med", "Control high"), (floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control P", 
floor(nr/3))))
cmap

# create 96-well plate control map
cmap = data.frame(X1=c(rep("Control P", floor(nr/3)), rep(c("Control low", "Control med", "Control high"), 
(floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control N", floor(nr/3))), X2=c(rep("Control N", floor(nr/3)), 
rep(c("Control low", "Control med", "Control high"), (floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control P", 
floor(nr/3))))
cmap = cmap[seq(1,nr,2),]
cmap

## ----echo=TRUE-----------------------------------------------------------
set.seed(1234)
library(highSCREEN)
nc = 24
nr = 16

# create 1st replicate of data matrix with compounds and controls
replicate1 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 2nd replicate of data matrix with compounds and controls
replicate2 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 3rd replicate of data matrix with compounds and controls
replicate3 = matrix(abs(rnorm(nr*nc)), nr, nc)

# combine all replicate to create "before" data set
replicates_before = list(replicate1, replicate2, replicate3)
names(replicates_before) = c("Replicate1", "Replicate2", "Replicate3")

# create 1st replicate of data matrix with compounds and controls
replicate1 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 2nd replicate of data matrix with compounds and controls
replicate2 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 3rd replicate of data matrix with compounds and controls
replicate3 = matrix(abs(rnorm(nr*nc)), nr, nc)

# combine all replicate to create "after" data set
replicates_after = list(replicate1, replicate2, replicate3)
names(replicates_after) = c("Replicate1", "Replicate2", "Replicate3")

# extract plate 3, replicate 2
extractplate(replicates_before, replicates_after, plate=3, replicate=2)

## ----echo=TRUE-----------------------------------------------------------
set.seed(1234)
library(highSCREEN)
nc = 24
nr = 16

# create control map
cmap = data.frame(X1=c(rep("Control P", floor(nr/3)), rep(c("Control low", "Control med", "Control high"), 
(floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control N", floor(nr/3))), X2=c(rep("Control N", floor(nr/3)), 
rep(c("Control low", "Control med", "Control high"), (floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control P", floor(nr/3))))
cmap = cmap[seq(1,nr,2),]
cmap

# create 1st replicate of data matrix with compounds and controls
replicate1 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 2nd replicate of data matrix with compounds and controls
replicate2 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 3rd replicate of data matrix with compounds and controls
replicate3 = matrix(abs(rnorm(nr*nc)), nr, nc)

# combine all replicate for the before data
replicates_before = list(replicate1, replicate2, replicate3)
names(replicates_before) = c("Replicate1", "Replicate2", "Replicate3")

# create 1st replicate of data matrix with compounds and controls
replicate1 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 2nd replicate of data matrix with compounds and controls
replicate2 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 3rd replicate of data matrix with compounds and controls
replicate3 = matrix(abs(rnorm(nr*nc)), nr, nc)

# combine all replicate for the after data
replicates_after = list(replicate1, replicate2, replicate3)
names(replicates_after) = c("Replicate1", "Replicate2", "Replicate3")

# extract plate 1, replicate 1
dat = extractplate(replicates_before, replicates_after, plate=1, replicate=1)

# normalize using c-score
head(normplate("Main Plate 1", dat[["datbefore"]], dat[["datafter"]], cmap, plate=1, triplicate=1, norm="cscore",
 poscont="Control P", negcont="Control N"))

# normalize using b-score (medpolish)
head(normplate("Main Plate 1", dat[["datbefore"]], dat[["datafter"]], cmap, plate=1, triplicate=1, norm="bscore"))

# normalize using z-score
head(normplate("Main Plate 1", dat[["datbefore"]], dat[["datafter"]], cmap, plate=1, triplicate=1, norm="zscore"))

## ----echo=TRUE-----------------------------------------------------------
set.seed(1234)
library(highSCREEN)
nc = 24
nr = 16

# create 1st triplicate of data matrix with compounds and controls
replicate1 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create control map
cmap = data.frame(X1=c(rep("Control P", floor(nr/3)), rep(c("Control low", "Control med", "Control high"), 
(floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control N", floor(nr/3))), X2=c(rep("Control N", floor(nr/3)), 
rep(c("Control low", "Control med", "Control high"), (floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control P", floor(nr/3))))
cmap = cmap[seq(1,nr,2),]

# create 2nd triplicate of data matrix with compounds and controls
replicate2 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 3rd triplicate of data matrix with compounds and controls
replicate3 = matrix(abs(rnorm(nr*nc)), nr, nc)

# combine all triplicates for the before data
replicates_before = list(replicate1, replicate2, replicate3)
names(replicates_before) = c("Replicate1", "Replicate2", "Replicate3")

# create 1st triplicate of data matrix with compounds and controls
replicate1 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 2nd triplicate of data matrix with compounds and controls
replicate2 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 3rd triplicate of data matrix with compounds and controls
replicate3 = matrix(abs(rnorm(nr*nc)), nr, nc)

# combine all triplicates for the after data
replicates_after = list(replicate1, replicate2, replicate3)
names(replicates_after) = c("Replicate1", "Replicate2", "Replicate3")

# extract plate 1, all triplicates
dat1 = extractplate(replicates_before, replicates_after, plate=1, replicate=1)
dat2 = extractplate(replicates_before, replicates_after, plate=1, replicate=2)
dat3 = extractplate(replicates_before, replicates_after, plate=1, replicate=3)

# normalize data of all triplicates
res1 = normplate("Main Plate 1", dat1[["datbefore"]], dat1[["datafter"]], cmap, plate=1, triplicate=1, norm="zscore") 
res2 = normplate("Main Plate 1", dat2[["datbefore"]], dat2[["datafter"]], cmap, plate=1, triplicate=2, norm="zscore")
res3 = normplate("Main Plate 1", dat3[["datbefore"]], dat3[["datafter"]], cmap, plate=1, triplicate=3, norm="zscore")

# reformat data of all triplicates
head(formatRESULT(rbind(res1, res2, res3), triplicate="Triplicate", score="score", t="Time"))

## ----echo=TRUE-----------------------------------------------------------
set.seed(1234)
library(highSCREEN)
nc = 24
nr = 16

# create 1st triplicate of data matrix with compounds and controls
replicate1 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create control map
cmap = data.frame(X1=c(rep("Control P", floor(nr/3)), rep(c("Control low", "Control med", "Control high"), 
(floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control N", floor(nr/3))), X2=c(rep("Control N", floor(nr/3)), 
rep(c("Control low", "Control med", "Control high"), (floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control P", floor(nr/3))))
cmap = cmap[seq(1,nr,2),]

# create 2nd triplicate of data matrix with compounds and controls
replicate2 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 3rd triplicate of data matrix with compounds and controls
replicate3 = matrix(abs(rnorm(nr*nc)), nr, nc)

# combine all triplicates for the before data
replicates_before = list(replicate1, replicate2, replicate3)

# create 1st triplicate of data matrix with compounds and controls
replicate1 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 2nd triplicate of data matrix with compounds and controls
replicate2 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 3rd triplicate of data matrix with compounds and controls
replicate3 = matrix(abs(rnorm(nr*nc)), nr, nc)

# combine all triplicates for the after data
replicates_after = list(replicate1, replicate2, replicate3)
names(replicates_after) = c("Replicate1", "Replicate2", "Replicate3")

# extract plate 1, replicate 1
dat11 = extractplate(replicates_before, replicates_after, plate=1, replicate=1)

# extract plate 1, triplicate 2
dat12 = extractplate(replicates_before, replicates_after, plate=1, replicate=2)

# extract plate 1, triplicate 3
dat13 = extractplate(replicates_before, replicates_after, plate=1, replicate=3)

# no normalizion (norm="raw")
res11 = normplate("Main Plate 1", dat11[["datbefore"]], dat11[["datafter"]], cmap, plate=1, triplicate=1, norm="raw")
res12 = normplate("Main Plate 1", dat12[["datbefore"]], dat12[["datafter"]], cmap, plate=1, triplicate=2, norm="raw")
res13 = normplate("Main Pltae 1", dat13[["datbefore"]], dat13[["datafter"]], cmap, plate=1, triplicate=3, norm="raw")

# combine 3 replicates
res1 = rbind(res11, res12, res13)

# reformat result
res1 = formatRESULT(res1, triplicate="Triplicate", score="score", t="Time")

# perform QC
qcplate(res1, poscont="Control P", negcont="Control N", qc1.val=0.225, qc2.val=2, 
addcont=c("Control low", "Control med", "Control high"), welltype="welltype")

## ----echo=TRUE-----------------------------------------------------------
set.seed(1234)
library(highSCREEN)

nc = 24
nr = 16

# create 1st triplicate of data matrix with compounds and controls
replicate1 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create control map
cmap = data.frame(X1=c(rep("Control P", floor(nr/3)), rep(c("Control low", "Control med", "Control high"), 
(floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control N", floor(nr/3))), X2=c(rep("Control N", floor(nr/3)), 
rep(c("Control low", "Control med", "Control high"), (floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control P", floor(nr/3))))
cmap = cmap[seq(1,nr,2),]


# create 2nd triplicate of data matrix with compounds and controls
replicate2 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 3rd triplicate of data matrix with compounds and controls
replicate3 = matrix(abs(rnorm(nr*nc)), nr, nc)

# combine all triplicates for the before data
replicates_before = list(replicate1, replicate2, replicate3)
names(replicates_before) = c("Replicate1", "Replicate2", "Replicate3")

# create 1st triplicate of data matrix with compounds and controls
replicate1 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 2nd triplicate of data matrix with compounds and controls
replicate2 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 3rd triplicate of data matrix with compounds and controls
replicate3 = matrix(abs(rnorm(nr*nc)), nr, nc)

# combine all triplicates for the after data
replicates_after = list(replicate1, replicate2, replicate3)
names(replicates_after) = c("Replicate1", "Replicate2", "Replicate3")

# extract plate 1, triplicate 1
dat1 = extractplate(replicates_before, replicates_after, plate=1, replicate=1)
# extract plate 1, triplicate 2
dat2 = extractplate(replicates_before, replicates_after, plate=1, replicate=2)
# extract plate 1, triplicate 3
dat3 = extractplate(replicates_before, replicates_after, plate=1, replicate=3)

# no normalizion
datraw1 = normplate("Main Plate 1", dat1[["datbefore"]], dat1[["datafter"]], cmap, plate=1, triplicate=1, norm="raw")
datraw2 = normplate("Main Plate 1", dat2[["datbefore"]], dat2[["datafter"]], cmap, plate=1, triplicate=2, norm="raw")
datraw3 = normplate("Main Pltae 1", dat3[["datbefore"]], dat3[["datafter"]], cmap, plate=1, triplicate=3, norm="raw")

# combine 3 triplicates
datraw = rbind(datraw1, datraw2, datraw3)

# reformat result
datraw = formatRESULT(datraw, triplicate="Triplicate", score="score", t="Time")

# compute z-factor and ssmd for each raw compound, triplicate 1
zfactor.ssmd(datraw, "Control P", "Control N", "Main Plate 1", 1)

## ----echo=TRUE-----------------------------------------------------------
set.seed(1234)
library(highSCREEN)

nc = 24
nr = 16

# create 1st triplicate of data matrix with compounds and controls
replicate1 = matrix(abs(rnorm(nr*nc)*0.01), nr, nc)

# create control map
cmap = data.frame(X1=c(rep("Control P", floor(nr/3)), rep(c("Control low", "Control med", "Control high"), 
(floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control N", floor(nr/3))), X2=c(rep("Control N", floor(nr/3)), 
rep(c("Control low", "Control med", "Control high"), (floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control P", floor(nr/3))))
cmap = cmap[seq(1,nr,2),]

# create 2nd triplicate of data matrix with compounds and controls
replicate2 = matrix(abs(rnorm(nr*nc)*0.01), nr, nc)

# create 3rd triplicate of data matrix with compounds and controls
replicate3 = matrix(abs(rnorm(nr*nc)*0.01), nr, nc)

# combine all triplicates for the before data
replicates_before = list(replicate1, replicate2, replicate3)
names(replicates_before) = c("Replicate1", "Replicate2", "Replicate3")

# create 1st triplicate of data matrix with compounds and controls
replicate1 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 2nd triplicate of data matrix with compounds and controls
replicate2 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 3rd triplicate of data matrix with compounds and controls
replicate3 = matrix(abs(rnorm(nr*nc)), nr, nc)

# combine all triplicates for the after data
replicates_after = list(replicate1, replicate2, replicate3)
names(replicates_after) = c("Replicate1", "Replicate2", "Replicate3")

# extract plate 1, triplicate 1
dat1 = extractplate(replicates_before, replicates_after, plate=1, replicate=1)

# extract plate 1, triplicate 2
dat2 = extractplate(replicates_before, replicates_after, plate=1, replicate=2)

# extract plate 1, triplicate 3
dat3 = extractplate(replicates_before, replicates_after, plate=1, replicate=3)

# no normalizion
datraw1 = normplate("Main Plate 1", dat1[["datbefore"]], dat1[["datafter"]], cmap, plate=1, triplicate=1, norm="raw")
datraw2 = normplate("Main Plate 1", dat2[["datbefore"]], dat2[["datafter"]], cmap, plate=1, triplicate=2, norm="raw")
datraw3 = normplate("Main Pltae 1", dat3[["datbefore"]], dat3[["datafter"]], cmap, plate=1, triplicate=3, norm="raw")

# combine 3 triplicates
datraw = rbind(datraw1, datraw2, datraw3)

# reformat result
datraw = formatRESULT(datraw, triplicate="Triplicate", score="score", t="Time")

# c-score normalization
datnorm1 = normplate("Main Plate 1", dat1[["datbefore"]], dat1[["datafter"]], cmap, plate=1, triplicate=1, norm="cscore", 
poscont="Control P", negcont="Control N")
datnorm2 = normplate("Main Plate 1", dat2[["datbefore"]], dat2[["datafter"]], cmap, plate=1, triplicate=2, norm="cscore", 
poscont="Control P", negcont="Control N")
datnorm3 = normplate("Main Pltae 1", dat3[["datbefore"]], dat3[["datafter"]], cmap, plate=1, triplicate=3, norm="cscore",
poscont="Control P", negcont="Control N")

# combine 3 triplicates
datnorm = rbind(datnorm1, datnorm2, datnorm3)

# reformat result
datnorm = formatRESULT(datnorm, triplicate="Triplicate", score="score", t="Time")

# identify hits
head(hits(datraw, datnorm, qc.mainplates="Main Plate 1", qc1.val=0.225, hit.val=3))

## ----echo=TRUE-----------------------------------------------------------
set.seed(1234)
library(highSCREEN)

nc = 24
nr = 16

# create 1st triplicate of data matrix with compounds and controls
replicate1 = matrix(abs(rnorm(nr*nc)*0.01), nr, nc)

# create control map
cmap = data.frame(X1=c(rep("Control P", floor(nr/3)), rep(c("Control low", "Control med", "Control high"), 
(floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control N", floor(nr/3))), X2=c(rep("Control N", floor(nr/3)), 
rep(c("Control low", "Control med", "Control high"), (floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control P", floor(nr/3))))
cmap = cmap[seq(1,nr,2),]

# create 2nd triplicate of data matrix with compounds and controls
replicate2 = matrix(abs(rnorm(nr*nc)*0.01), nr, nc)

# create 3rd triplicate of data matrix with compounds and controls
replicate3 = matrix(abs(rnorm(nr*nc)*0.01), nr, nc)

# combine all triplicates for the before data
replicates_before = list(replicate1, replicate2, replicate3)
names(replicates_before) = c("Replicate1", "Replicate2", "Replicate3")

# create 1st triplicate of data matrix with compounds and controls
replicate1 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 2nd triplicate of data matrix with compounds and controls
replicate2 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 3rd triplicate of data matrix with compounds and controls
replicate3 = matrix(abs(rnorm(nr*nc)), nr, nc)

# combine all triplicates for the after data
replicates_after = list(replicate1, replicate2, replicate3)
names(replicates_after) = c("Replicate1", "Replicate2", "Replicate3")

# extract plate 1, replicate 1
dat1 = extractplate(replicates_before, replicates_after, plate=1, replicate=1)
# extract plate 1, replicate 2
dat2 = extractplate(replicates_before, replicates_after, plate=1, replicate=2)
# extract plate 1, replicate 3
dat3 = extractplate(replicates_before, replicates_after, plate=1, replicate=3)

# no normalizion
datraw1 = normplate("Main Plate 1", dat1[["datbefore"]], dat1[["datafter"]], cmap, plate=1, triplicate=1, norm="raw")
datraw2 = normplate("Main Plate 1", dat2[["datbefore"]], dat2[["datafter"]], cmap, plate=1, triplicate=2, norm="raw")
datraw3 = normplate("Main Pltae 1", dat3[["datbefore"]], dat3[["datafter"]], cmap, plate=1, triplicate=3, norm="raw")

# combine 3 triplicates
datraw = rbind(datraw1, datraw2, datraw3)

# reformat result
datraw = formatRESULT(datraw, triplicate="Triplicate", score="score", t="Time")

# c-score normalization
datnorm1 = normplate("Main Plate 1", dat1[["datbefore"]], dat1[["datafter"]], cmap, plate=1, triplicate=1, norm="cscore", 
poscont="Control P", negcont="Control N")
datnorm2 = normplate("Main Plate 1", dat2[["datbefore"]], dat2[["datafter"]], cmap, plate=1, triplicate=2, norm="cscore", 
poscont="Control P", negcont="Control N")
datnorm3 = normplate("Main Pltae 1", dat3[["datbefore"]], dat3[["datafter"]], cmap, plate=1, triplicate=3, norm="cscore", 
poscont="Control P", negcont="Control N")

# combine 3 triplicates
datnorm = rbind(datnorm1, datnorm2, datnorm3)

# reformat result
datnorm = formatRESULT(datnorm, triplicate="Triplicate", score="score", t="Time")

# identify hits
h = hits(datraw, datnorm, qc.mainplates="Main Plate 1", qc1.val=0.225, hit.val=3)

# rank hits in disending order of mean of "after" replicate scores "ma"
head(rankhits(h))

## ----echo=TRUE-----------------------------------------------------------
set.seed(1234)
library(highSCREEN)
library(gplots)

nc = 24
nr = 16

# create 1st triplicate of data matrix with compounds and controls
replicate1 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create control map
cmap = data.frame(X1=c(rep("Control P", floor(nr/3)), rep(c("Control low", "Control med", "Control high"), 
(floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control N", floor(nr/3))), X2=c(rep("Control N", floor(nr/3)), 
rep(c("Control low", "Control med", "Control high"), (floor(nr/3)+nr-3*floor(nr/3))/3), rep("Control P", floor(nr/3))))
cmap = cmap[seq(1,nr,2),]

# create 2nd triplicate of data matrix with compounds and controls
replicate2 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 3rd triplicate of data matrix with compounds and controls
replicate3 = matrix(abs(rnorm(nr*nc)), nr, nc)

# combine all triplicates for the before data
replicates_before = list(replicate1, replicate2, replicate3)
names(replicates_before) = c("Replicate1", "Replicate2", "Replicate3")

# create 1st triplicate of data matrix with compounds and controls
replicate1 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 2nd triplicate of data matrix with compounds and controls
replicate2 = matrix(abs(rnorm(nr*nc)), nr, nc)

# create 3rd triplicate of data matrix with compounds and controls
replicate3 = matrix(abs(rnorm(nr*nc)), nr, nc)

# combine all triplicates for the after data
replicates_after = list(replicate1, replicate2, replicate3)
names(replicates_after) = c("Replicate1", "Replicate2", "Replicate3")

# extract plate 1, triplicate 1
dat11 = extractplate(replicates_before, replicates_after, plate=1, replicate=1)
# extract plate 1, triplicate 2
dat12 = extractplate(replicates_before, replicates_after, plate=1, replicate=2)
# extract plate 1, triplicate 3
dat13 = extractplate(replicates_before, replicates_after, plate=1, replicate=3)

# no normalizion (norm="raw")
res11 = normplate("Main Plate 1", dat11[["datbefore"]], dat11[["datafter"]], cmap, plate=1, triplicate=1, norm="raw")
res12 = normplate("Main Plate 1", dat12[["datbefore"]], dat12[["datafter"]], cmap, plate=1, triplicate=2, norm="raw")
res13 = normplate("Main Pltae 1", dat13[["datbefore"]], dat13[["datafter"]], cmap, plate=1, triplicate=3, norm="raw")

# combine 3 triplicates
res1 = rbind(res11, res12, res13)
# reformat result
res1 = formatRESULT(res1, triplicate="Triplicate", score="score", t="Time")

layout(matrix(c(1,2,3), 3, 1, byrow = TRUE))

# plot density of all positive controls
plotcont(subset(res1, welltype=="Control P"), main="Density of Positive Controls", xaxis.marks=seq(-1,5,0.025))

# plot density of all negative controls
plotcont(subset(res1, welltype=="Control N"), main="Density of Negative Controls", xaxis.marks=seq(-1,5,0.025))

# plot density of controls with low, medium and high concentrations
plotcont(subset(res1, welltype=="Control low" | welltype=="Control med" | welltype=="Control high"), main="Density of Controls with Low, 
Medium and High Concentrations", xaxis.marks=seq(-1,5,0.025))

# plot single plate activity levels
plotplate(res1, main="Single Plate Activity Levels")

sessionInfo()

