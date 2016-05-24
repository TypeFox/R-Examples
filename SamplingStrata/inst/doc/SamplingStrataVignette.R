### R code from vignette source 'SamplingStrataVignette.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: foo
###################################################
options(keep.source = TRUE, width = 60)
foo <- packageDescription('SamplingStrata')


###################################################
### code chunk number 2: frame1
###################################################
library(SamplingStrata)
require(memoise)
data(swissmunicipalities)


###################################################
### code chunk number 3: frame2
###################################################
swissframe <- NULL
swissframe$id <- swissmunicipalities$Nom


###################################################
### code chunk number 4: frame3
###################################################
swissframe$Y1 <- swissmunicipalities$Pop020
swissframe$Y2 <- swissmunicipalities$Pop2040
swissframe$Y3 <- swissmunicipalities$Pop4065
swissframe$Y4 <- swissmunicipalities$Pop65P


###################################################
### code chunk number 5: frame4
###################################################
library(SamplingStrata)
swissframe$X1 <- var.bin(swissmunicipalities$POPTOT, bins=18)
swissframe$X2 <- var.bin(swissmunicipalities$Surfacesbois, bins=3)
swissframe$X3 <- var.bin(swissmunicipalities$Surfacescult, bins=3)
swissframe$X4 <- var.bin(swissmunicipalities$Alp, bins=3)
swissframe$X5 <- var.bin(swissmunicipalities$Airbat, bins=3)
swissframe$X6 <- var.bin(swissmunicipalities$Airind, bins=3)


###################################################
### code chunk number 6: frame5
###################################################
swissframe$domainvalue <- swissmunicipalities$REG
swissframe <- data.frame(swissframe)


###################################################
### code chunk number 7: frame6
###################################################
head(swissframe)


###################################################
### code chunk number 8: frame7
###################################################
write.table (swissframe, "swissframe.txt", row.names=FALSE,col.names=TRUE, sep="\t", quote=FALSE)


###################################################
### code chunk number 9: frame8
###################################################
library(SamplingStrata)
data(swissframe)
head(swissframe)


###################################################
### code chunk number 10: frame9
###################################################
data(strata)
head(strata)


###################################################
### code chunk number 11: frame11
###################################################
swissstrata <- buildStrataDF(swissframe)


###################################################
### code chunk number 12: frame12
###################################################
head(swissstrata)


###################################################
### code chunk number 13: frame13
###################################################
data(swissstrata)


###################################################
### code chunk number 14: frame14 (eval = FALSE)
###################################################
## samp <- read.delim("samplePrev.txt")


###################################################
### code chunk number 15: frame15 (eval = FALSE)
###################################################
## strata <- buildStrataDF(samp)


###################################################
### code chunk number 16: frame16
###################################################
data(swisserrors)
swisserrors


###################################################
### code chunk number 17: frame17
###################################################
checkInput(swisserrors,swissstrata,swissframe)


###################################################
### code chunk number 18: frame18
###################################################
cv <- swisserrors[1,]
cv


###################################################
### code chunk number 19: frame19
###################################################
sum(bethel(swissstrata,cv))


###################################################
### code chunk number 20: frame20
###################################################
solution <- optimizeStrata(
	errors = swisserrors, 
	strata = swissstrata, 
	cens = NULL, 
	strcens = FALSE, 
	initialStrata = nrow(strata), 
	addStrataFactor = 0.00, 
	minnumstr = 2, 
	iter = 40, 
	pops = 10, 
	mut_chance = 0.05, 
	elitism_rate = 0.2,
	highvalue = 1e+08, 
	suggestions = NULL,
	realAllocation = TRUE,
	writeFiles = TRUE)
sum(ceiling(solution$aggr_strata$SOLUZ))


###################################################
### code chunk number 21: frame21
###################################################
newstrata <- updateStrata(swissstrata, solution, writeFiles = TRUE)


###################################################
### code chunk number 22: frame13
###################################################
strata_aggregation <- read.delim("strata_aggregation.txt")
head(strata_aggregation)


###################################################
### code chunk number 23: frame22
###################################################
framenew <- updateFrame(swissframe, newstrata, writeFiles=TRUE)


###################################################
### code chunk number 24: frame23
###################################################
sample <- selectSample(framenew, solution$aggr_strata, writeFiles=TRUE)


###################################################
### code chunk number 25: frame24
###################################################
evalSolution(framenew, solution$aggr_strata, nsampl=50, writeFiles=TRUE) 


###################################################
### code chunk number 26: frame25
###################################################
expected_cv <- read.csv("expected_cv.csv")
expected_cv


###################################################
### code chunk number 27: frame24
###################################################
data(swisserrors)
data(swissstrata)
data(swissframe)
#----Selection of units to be censused from the frame
framecens <- swissframe[ (swissframe$domainvalue == 1 |
                          swissframe$domainvalue == 4) & 
                         (swissframe$X2 == 1 &
                          swissframe$X3 == 1 &
                          swissframe$X4 == 1 &
                          swissframe$X5 == 1 &
                          swissframe$X6 == 1)  , ]
#----Selection of units to be sampled from the frame
# (complement to the previous)
framesamp <- swissframe[!((swissframe$domainvalue == 1 |
                           swissframe$domainvalue == 4) & 
                          (swissframe$X2 == 1 &
                           swissframe$X3 == 1 &
                           swissframe$X4 == 1 &
                           swissframe$X5 == 1 &
                           swissframe$X6 == 1)) , ]


###################################################
### code chunk number 28: frame24
###################################################
# Build strata to be censused and sampled
cens <- buildStrataDF(framecens)
sum(cens$N)
strata <- buildStrataDF(framesamp)
sum(strata$N)


###################################################
### code chunk number 29: frame25
###################################################
solution <- optimizeStrata(
	errors = swisserrors, 
	strata = strata, 
	cens = cens, 
	strcens = TRUE, 
	initialStrata = nrow(strata), 
	addStrataFactor = 0.00, 
	minnumstr = 2, 
	iter = 40, 
	pops = 10, 
	mut_chance = 0.05, 
	elitism_rate = 0.2,
	highvalue = 1e+08, 
	suggestions = NULL,
	realAllocation = TRUE,
	writeFiles = TRUE)


###################################################
### code chunk number 30: frame25
###################################################
newstrata <- updateStrata(strata, solution)
# updating sampling frame with new strata labels
framenew <- updateFrame(frame=framesamp,newstrata=newstrata)
# selection of sample from sampling strata
sample <- selectSample(frame=framenew,outstrata=solution$aggr_strata)


###################################################
### code chunk number 31: frame26
###################################################
# addition of necessary variables to 
colnames(framesamp) <- toupper(colnames(framesamp))
colnames(framecens) <- toupper(colnames(framecens))
framecens$WEIGHTS <- rep(1,nrow(framecens))
framecens$FPC <- rep(1,nrow(framecens))
framecens$LABEL <- rep("999999",nrow(framecens))
framecens$STRATUM <- rep("999999",nrow(framecens))
framecens$STRATO <- rep("999999",nrow(framecens))


###################################################
### code chunk number 32: frame27
###################################################
survey <- rbind(sample,framecens)


###################################################
### code chunk number 33: frame28
###################################################
survey$cens <- ifelse(survey$LABEL == "999999",1,0)
table(survey$cens)


