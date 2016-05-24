# Test cases for function 'anovaVCA' of R-package VCA
# 
# Author: André Schützenmeister
#
# Reference results were obtained by application of SAS 9.2 PROC NESTED to 
# test data, SAS PROC VARCOMP or SAS PROC MIXED method=type1.
#
###############################################################################


cat("\n\n***********************************************************************")
cat("\nVariance Component Analysis (VCA) - test cases for function 'anovaVCA'.")
cat("\n***********************************************************************\n\n")


### load all testdata

data(dataEP05A2_1)
data(dataEP05A2_2)
data(dataEP05A2_3)

data(dataEP05A3_MS_1)
data(dataEP05A3_MS_2)
data(dataEP05A3_MS_3)

data(dataRS0003_1)
data(dataRS0003_2)
data(dataRS0003_3)

data(dataRS0005_1)
data(dataRS0005_2)
data(dataRS0005_3)


### check results of function 'anovaVCA' ###


# check whether function 'anovaVCA' throws exceptions as expected

TF001.anovaVCA.exception <- function()
{
	checkException(anovaVCA())                                                               # no input at all 
	checkException(anovaVCA(Data=1))                                                                                               
	checkException(anovaVCA(Data=data.frame()))                 
	checkException(anovaVCA(Data=data.frame(y=1:10)))
	checkException(anovaVCA(z~day/run, Data=data.frame(y=1:10)))
	checkException(anovaVCA(y~day/run, Data=data.frame(y=1:10, day=1:10)))
}


## EP05-A2 20/2/2 Within-Lab Precision Experiments (Single-Site Experiment in EP05-A3) - balanced

TF002.anovaVCA.EP05_A2_intermediate_precision.balanced <- function()
{        
	res <- anovaVCA(y~day/run, Data=dataEP05A2_1)                                            # call function    
	checkEquals(round(as.numeric(res$aov.tab[,"VC"]),6), c(3.298986, 0.001096, 1.344492, 1.953397))     # variance components are derived from mean sum of squares --> checking SSQ und MSQ not necessary
	
	res <- anovaVCA(y~day/run, Data=dataEP05A2_2)                                            # call function    
	checkEquals(round(as.numeric(res$aov.tab[,"VC"]),6), c(8.400103, 1.853772, 2.826050, 3.720281))    # variance components are derived from mean sum of squares --> checking SSQ und MSQ not necessary
	
	res <- anovaVCA(y~day/run, Data=dataEP05A2_3)                                            # call function    
	checkEquals(round(as.numeric(res$aov.tab[,"VC"]),6), c( 35.546313, 12.231099, 7.031193, 16.284021)) # variance components are derived from mean sum of squares --> checking SSQ und MSQ not necessary
}    
# unbalanced

TF003.anovaVCA.EP05_A2_intermediate_precision.unbalanced <- function()
{ 
	res <- anovaVCA(y~day/run, Data=dataEP05A2_1[-c(11, 12, 17, 37, 45, 56, 57, 68),])              
	checkEquals(round(as.numeric(res$aov.tab[,"VC"]),6), c(3.272031, 0.278178, 0.875624, 2.118229))     
	
	res <- anovaVCA(y~day/run, Data=dataEP05A2_2[-c(2, 12, 22, 23, 24, 55, 56, 71),])              
	checkEquals(round(as.numeric(res$aov.tab[,"VC"]),6), c(8.81105, 1.425126, 3.478988, 3.906936 ))  
	
	res <- anovaVCA(y~day/run, Data=dataEP05A2_3[-c(1,6,7,36,61:65),])              
	checkEquals(round(as.numeric(res$aov.tab[,"VC"]),6), c( 33.701053, 10.384121, 7.937951, 15.378981 )) 
}

# check whether the reported SD-values are correctly computed by comparing the square-root values of
# VC-values to values of column "SD"

TF004.anovaVCA.SD_results <- function()
{ 
	res <- anovaVCA(y~day/run, Data=dataEP05A2_1[-c(11, 12, 17, 37, 45, 56, 57, 68),])              
	checkEquals(round(sqrt(as.numeric(res$aov.tab[,"VC"])),6), round(as.numeric(res$aov.tab[,"SD"]), 6))     
	
	res <- anovaVCA(y~day/run, Data=dataEP05A2_2[-c(2, 12, 22, 23, 24, 55, 56, 71),])              
	checkEquals(round(sqrt(as.numeric(res$aov.tab[,"VC"])),6), round(as.numeric(res$aov.tab[,"SD"]), 6)) 
	
	res <- anovaVCA(y~day/run, Data=dataEP05A2_3[-c(1,6,7,36,61:65),])              
	checkEquals(round(sqrt(as.numeric(res$aov.tab[,"VC"])),6), round(as.numeric(res$aov.tab[,"SD"]), 6)) 
}


# check whether the reported CV-values are correctly computed by comparing 100 times square-root values 
# of the VC-value devided by the mean to values of column "CV[%]"

TF005.anovaVCA.CV_perc_results <- function()
{ 
	res <- anovaVCA(y~day/run, Data=dataEP05A2_1[-c(11, 12, 17, 37, 45, 56, 57, 68),])              
	checkEquals(round(sqrt(as.numeric(res$aov.tab[,"VC"]))*100/res$Mean,6), round(as.numeric(res$aov.tab[,"CV[%]"]), 6))     
	
	res <- anovaVCA(y~day/run, Data=dataEP05A2_2[-c(2, 12, 22, 23, 24, 55, 56, 71),])              
	checkEquals(round(sqrt(as.numeric(res$aov.tab[,"VC"]))*100/res$Mean,6), round(as.numeric(res$aov.tab[,"CV[%]"]), 6)) 
	
	res <- anovaVCA(y~day/run, Data=dataEP05A2_3[-c(1,6,7,36,61:65),])              
	checkEquals(round(sqrt(as.numeric(res$aov.tab[,"VC"]))*100/res$Mean,6), round(as.numeric(res$aov.tab[,"CV[%]"]), 6)) 
}


# check whether values in column "%Total" were correctly computed

TF006.anovaVCA.Percent_Total_results <- function()
{ 
	res <- anovaVCA(y~day/run, Data=dataEP05A2_1)              
	checkEquals(round(as.numeric(res$aov.tab[,"VC"])*100/sum(as.numeric(res$aov.tab[-1, "VC"])),6), round(as.numeric(res$aov.tab[,"%Total"]), 6))           # exclude total variance in the sum  
	
	res <- anovaVCA(y~day/run, Data=dataEP05A2_2)              
	checkEquals(round(as.numeric(res$aov.tab[,"VC"])*100/sum(as.numeric(res$aov.tab[-1, "VC"])),6), round(as.numeric(res$aov.tab[,"%Total"]), 6)) 
	
	res <- anovaVCA(y~day/run, Data=dataEP05A2_3)              
	checkEquals(round(as.numeric(res$aov.tab[,"VC"])*100/sum(as.numeric(res$aov.tab[-1, "VC"])),6), round(as.numeric(res$aov.tab[,"%Total"]), 6)) 
}


## EP05-A3 3/5/5 Multi-Site Experiments - balanced

TF007.anovaVCA.EP05_A3_Reproducibility.balanced <- function()
{      
	res <- anovaVCA(y~site/day, Data=dataEP05A3_MS_1)
	checkEquals(round(as.numeric(res$aov.tab[,"VC"]), 6), c( 3.635232, 1.017401, 0.345882, 2.271949))
	
	res <- anovaVCA(y~site/day, Data=dataEP05A3_MS_2)
	checkEquals(round(as.numeric(res$aov.tab[,"VC"]), 6), c(5.765516, 1.010798, 1.028309, 3.726408))
	
	res <- anovaVCA(y~site/day, Data=dataEP05A3_MS_3)
	checkEquals(round(as.numeric(res$aov.tab[,"VC"]), 6), c(38.438903, 17.990790, 5.423654, 15.024459))
}   
# unbalanced

TF008.anovaVCA.EP05_A3_Reproducibility.unbalanced <- function()
{
	res <- anovaVCA(y~site/day, Data=dataEP05A3_MS_1[-c(11, 12, 17, 37, 45, 56, 57, 68), ])
	checkEquals(round(as.numeric(res$aov.tab[,"VC"]), 6), c( 3.244177, 0.785705, 0.396764, 2.061709 ))
	
	res <- anovaVCA(y~site/day, Data=dataEP05A3_MS_2[-c(2, 12, 22, 23, 24, 55, 56, 71), ])
	checkEquals(round(as.numeric(res$aov.tab[,"VC"]), 6), c( 5.773996, 0.729652, 0.917329, 4.127015 ))
	
	res <- anovaVCA(y~site/day, Data=dataEP05A3_MS_3[-c(1,6,7,36,61:65),])
	checkEquals(round(as.numeric(res$aov.tab[,"VC"]), 6), c(37.092294, 16.872368, 4.773611, 15.446314))
}


## RS0003 21 replications, testing only the residual error component
## WRI = within run imprecision

TF009.anovaVCA.WRI <- function()
{ 
	res <- anovaVCA(y~1, Data=dataRS0003_1)
	checkEquals( round(as.numeric(res$aov.tab[2,c("SS", "MS")]),5), c(22.44452, 1.12223) )
	
	res <- anovaVCA(y~1, Data=dataRS0003_2)
	checkEquals( round(as.numeric(res$aov.tab[2,c("SS", "MS")]),4), c(367.3480, 18.3674) )
	
	res <- anovaVCA(y~1, Data=dataRS0003_3)
	checkEquals( round(as.numeric(res$aov.tab[2,c("SS", "MS")]),2), c(2210.59, 110.53) )
}    

## RS0005 - Confirmation of internal data by external labs - balanced
## BDI = between day imprecision

TF010.anovaVCA.BDI.external_labs.balanced <- function()
{     
	res <- anovaVCA(y~day, Data=dataRS0005_1, NegVC=TRUE)
	checkEquals( round(as.numeric(res$aov.tab[,"VC"]), 6), c(5.193341, 1.818128, 3.375212))
	
	res <- anovaVCA(y~day, Data=dataRS0005_2, NegVC=TRUE)
	checkEquals( round(as.numeric(res$aov.tab[,"VC"]), 6), c(6.611960, 0.055541, 6.556419))
	
	res <- anovaVCA(y~day, Data=dataRS0005_3, NegVC=TRUE)
	checkEquals( round(as.numeric(res$aov.tab[,"VC"]), 5), c(75.69922, 22.45171, 53.24751))
}   

# unbalanced

TF011.anovaVCA.BDI.external_labs.unbalanced <- function()
{
	res <- anovaVCA(y~day, Data=dataRS0005_1[-c(1,10),], NegVC=TRUE)
	checkEquals( round(as.numeric(res$aov.tab[,"VC"]), 6), c(5.460289, 3.592088, 1.868201))
	
	res <- anovaVCA(y~day, Data=dataRS0005_2[-c(4,5,15),], NegVC=FALSE)                              # NegVC=FALSE, because SAS PROC NESTED does not count negative VCs to total VC
	checkEquals( round(as.numeric(res$aov.tab[,"VC"]), 6), c( 8.548521, 0, 8.548521))
	
	res <- anovaVCA(y~day, Data=dataRS0005_3[-c(2,7,9,14),], NegVC=FALSE)
	checkEquals( round(as.numeric(res$aov.tab[,"VC"]), 6), c( 62.404131, 0, 62.404131 ))
}

# test whether conversion of variables to factors is correct

TF012.anovaVCA.variable_conversion <- function()
{
	data(dataEP05A2_1)
	dat1 <- dataEP05A2_1
	dat1$day <- as.integer(as.character(dat1$day))
	dat1$run <- as.integer(as.character(dat1$run))
	
	res1 <- anovaVCA(y~day/run, dat1)
	res2 <- anovaVCA(y~day/run, dataEP05A2_1)
	
	checkEquals( res1$aov.tab, res2$aov.tab )
}

# test whether missing values in response and other variables are correctly handled

TF013.anovaVCA.missing_value_handling <- function()
{
	data(dataEP05A2_3)
	dat0 <- dataEP05A2_3
	
	dat0[c(5,15,25),    "y"] <- NA                      # generated missing data
	dat0[c(3,17,58),  "day"] <- NA
	dat0[c(51,70,77), "run"] <- NA
	
	datNoNA <- na.omit(dat0)
	
	res1 <- anovaVCA(y~day/run, dat0)
	res2 <- anovaVCA(y~day/run, datNoNA)
	
	checkEquals(as.matrix(res1$aov.tab), as.matrix(res2$aov.tab))
	checkEquals(res1$Nrm,   9)
	checkEquals(res1$Nobs, 71)
	checkEquals(res2$Nobs, 71)
}


# check numerical equivalence of a crossed-nested design (balanced)
#
# check vs. SAS PROC MIXED results using this SAS-code:
#
#   proc mixed data=sample1 method=type1 cl;
#       class lot device day run;
#       model y=;
#       random lot device lot*device*day lot*device*day*run;
#   run;

TF014.anovaVCA.crossed_nested.balanced <- function()
{
	data(VCAdata1)
	sample1 <- VCAdata1[which(VCAdata1$sample==1),]
	sample1$device <- gl(3,28,252)                                      # add device variable
	set.seed(505)
	sample1$y <- sample1$y + rep(rep(rnorm(3,,.25), c(28,28,28)),3)     # add error component according to 
	
	# write.table(sample1, file="sample1.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")      # export data, import to SAS and apply SAS-code given above
	
	res1 <- anovaVCA(y~lot+device+(lot:device:day)/run, sample1)
	
	checkEquals(round(res1$aov.tab[2, "VC"], 5), 0.01552)                        # round to precision of SAS PROC MIXED output
	checkEquals(round(res1$aov.tab[3, "VC"], 5), 0.06214)
	checkEquals(round(res1$aov.tab[4, "VC"], 5), 0.01256)
	checkEquals(round(res1$aov.tab[5, "VC"], 5), 0.05074)
	checkEquals(round(res1$aov.tab[6, "VC"], 6), 0.001152)
	
}

# now test numerical equivalence to SAS PROC MIXED results with unbalanced data 
# using following SAS-code:
#   
#   data sample1UB;
#       set sample1;
#       obs = _N_;
#       if obs in (5,31,55) then delete;
#   run;
#   
#   proc mixed data=sample1UB method=type1 asycov cl;
#       class lot device day run;
#       model y=;
#       random lot device lot*device*day lot*device*day*run;
#   run;

TF015.anovaVCA.crossed_nested.unbalanced <- function()
{
	data(VCAdata1)
	sample1 <- VCAdata1[which(VCAdata1$sample==1),]
	sample1$device <- gl(3,28,252)                                      # add device variable
	set.seed(505)
	sample1$y <- sample1$y + rep(rep(rnorm(3,,.25), c(28,28,28)),3)     # add error component according to 
	
	sample1UB <- sample1[-c(5,31,55),]                                  # delete some observations
	
	# write.table(sample1UB, file="sample1UB.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")      # export data, import to SAS and apply SAS-code given above
	
	res1UB <- anovaVCA(y~lot+device+(lot:device:day)/run, sample1UB)
	
	checkEquals(round(res1UB$aov.tab[2, "VC"], 5), 0.01507)                        # round to precision of SAS PROC MIXED output
	checkEquals(round(res1UB$aov.tab[3, "VC"], 5), 0.06284)
	checkEquals(round(res1UB$aov.tab[4, "VC"], 5), 0.01259)
	checkEquals(round(res1UB$aov.tab[5, "VC"], 5), 0.05143)
	checkEquals(round(res1UB$aov.tab[6, "VC"], 6), 0.001145)
	
}




# use subset "sample_8" of VCAdata1

TF016.anovaVCA.crossed_nested.balanced <- function()
{
	data(VCAdata1)
	sample8 <- VCAdata1[which(VCAdata1$sample==8),]
	sample8$device <- gl(3,28,252)                                      # add device variable
	set.seed(903211)
	sample8$y <- sample8$y + rep(rep(rnorm(3,,.75), c(28,28,28)),3)     # add error component according to 
	
	# write.table(sample8, file="sample8.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")      # export data, import to SAS and apply SAS-code given above
	
	res8 <- anovaVCA(y~lot+device+(lot:device:day)/run, sample8)
	
	checkEquals(round(res8$aov.tab[2, "VC"], 4), 7.3463)                        # round to precision of SAS PROC MIXED output
	checkEquals(round(res8$aov.tab[3, "VC"], 4), 2.2716)
	checkEquals(round(res8$aov.tab[4, "VC"], 4), 3.4588)
	checkEquals(round(res8$aov.tab[5, "VC"], 4), 2.0302)
	checkEquals(round(res8$aov.tab[6, "VC"], 4), 1.2219)
	
}


# use subset "sample_8" of VCAdata1 and generate unbalancedness

TF017.anovaVCA.crossed_nested.unbalanced <- function()
{
	data(VCAdata1)
	sample8 <- VCAdata1[which(VCAdata1$sample==8),]
	sample8$device <- gl(3,28,252)                                      # add device variable
	set.seed(903211)
	sample8$y <- sample8$y + rep(rep(rnorm(3,,.75), c(28,28,28)),3)     # add error component according to 
	sample8UB <- sample8[-c(7, 18, 19, 101, 191, 193, 202, 203, 204, 222),]
	
	# write.table(sample8UB, file="sample8UB.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")      # export data, import to SAS and apply SAS-code given above
	
	res8UB <- anovaVCA(y~lot+device+(lot:device:day)/run, sample8UB)
	
	checkEquals(round(res8UB$aov.tab[2, "VC"], 4), 7.2724)                        # round to precision of SAS PROC MIXED output
	checkEquals(round(res8UB$aov.tab[3, "VC"], 4), 2.2278)
	checkEquals(round(res8UB$aov.tab[4, "VC"], 4), 3.6385)
	checkEquals(round(res8UB$aov.tab[5, "VC"], 4), 1.9844)
	checkEquals(round(res8UB$aov.tab[6, "VC"], 4), 1.2451)
	
}


# checking the implemented strategy of obtaining the ANOVA Type-I DFs

TF018.anovaDF.balanced <- function()
{
	aov.fit <- anova(lm(y~day/run, dataEP05A2_1))
	rm.fit  <- anovaVCA(y~day/run, dataEP05A2_1)
	
	checkEquals(as.numeric(rm.fit$aov.tab[-1, "DF"]), as.numeric(aov.fit[,"Df"]))
}

TF019.anovaDF.unbalanced <- function()
{
	aov.fit <- anova(lm(y~day/run, dataEP05A2_1[-c(1,7,11,41,50:55),]))
	rm.fit  <- anovaVCA(y~day/run, dataEP05A2_1[-c(1,7,11,41,50:55),])
	
	checkEquals(as.numeric(rm.fit$aov.tab[-1, "DF"]), as.numeric(aov.fit[,"Df"]))
}

TF020.anovaDF.balanced <- function()
{
	aov.fit <- anova(lm(y~day/run-1, dataEP05A2_1))
	rm.fit  <- anovaVCA(y~day/run-1, dataEP05A2_1)
	
	checkEquals(as.numeric(rm.fit$aov.tab[-1, "DF"]), as.numeric(aov.fit[,"Df"]))
}

TF021.anovaDF.unbalanced <- function()
{
	aov.fit <- anova(lm(y~day/run-1, dataEP05A2_1[-c(1,7,11,41,50:55),]))
	rm.fit  <- anovaVCA(y~day/run-1, dataEP05A2_1[-c(1,7,11,41,50:55),])
	
	checkEquals(as.numeric(rm.fit$aov.tab[-1, "DF"]), as.numeric(aov.fit[,"Df"]))
}


# only interactions no main factors (actually violating marginality principle)

data(VCAdata1)
datS1to2 <- VCAdata1[VCAdata1$sample %in% 1:2,]
datS4to5 <- VCAdata1[VCAdata1$sample %in% 4:5,]
datS1to2.ub <- datS1to2[-c(8, 72, 85, 140, 152, 174, 219, 265, 284, 288, 294, 316, 324, 328, 338, 349, 432, 450, 456, 493),]
datS4to5.ub <- datS1to2[-c(56, 116, 170, 184, 211, 219, 221, 256, 257, 261, 309, 359, 376, 403, 432, 440, 459, 460, 474, 475),]


TF022.anovaDF.balanced <- function()
{
	aov.fit <- anova(lm(y~lot:sample+sample:day+sample:day:run, datS1to2))
	
	rm.fit  <- anovaVCA(y~lot:sample+sample:day+sample:day:run, datS1to2)
	
	checkEquals(as.numeric(rm.fit$aov.tab[-1, "DF"]), as.numeric(aov.fit[,"Df"]))
}

TF023.anovaDF.unbalanced <- function()
{
	aov.fit <- anova(lm(y~lot:sample+sample:day+sample:day:run-1, datS1to2))
	rm.fit  <- anovaVCA(y~lot:sample+sample:day+sample:day:run-1, datS1to2)
	
	checkEquals(as.numeric(rm.fit$aov.tab[-1, "DF"]), as.numeric(aov.fit[,"Df"]))
}

TF024.anovaDF.balanced <- function()
{
	aov.fit <- anova(lm(y~lot:sample+sample:day+sample:day:run, datS1to2.ub))
	rm.fit  <- anovaVCA(y~lot:sample+sample:day+sample:day:run, datS1to2.ub)
	
	checkEquals(as.numeric(rm.fit$aov.tab[-1, "DF"]), as.numeric(aov.fit[,"Df"]))
}

TF025.anovaDF.unbalanced <- function()
{
	aov.fit <- anova(lm(y~lot:sample+sample:day+sample:day:run-1, datS1to2.ub))
	rm.fit  <- anovaVCA(y~lot:sample+sample:day+sample:day:run-1, datS1to2.ub)
	
	checkEquals(as.numeric(rm.fit$aov.tab[-1, "DF"]), as.numeric(aov.fit[,"Df"]))
}

# checks whether the order of the model terms is kept

TF026.anovaDF.balanced.ordering <- function()
{
	aov.fit <- anova(lm(terms(y~lot+lot:day+sample+sample:day, keep.order=TRUE), datS1to2.ub))
	rm.fit  <- anovaVCA(y~lot+lot:day+sample+sample:day, datS1to2.ub)
	
	checkEquals(as.numeric(rm.fit$aov.tab[-1, "DF"]), as.numeric(aov.fit[,"Df"]))	
	
	aov.fit <- anova(lm(terms(y~lot+lot:day+sample+sample:day, keep.order=FALSE), datS1to2.ub))			# ordering not kept --> differences
	rm.fit  <- anovaVCA(y~lot+lot:day+sample+sample:day, datS1to2.ub)
	
	checkTrue(any(as.numeric(rm.fit$aov.tab[-1, "DF"]) != as.numeric(aov.fit[,"Df"])))
}


# check numerical equivalence of sweep-based and quadratic form-based SSQ computation to reference results

TF027.check.sweep_vs_qf_ssq_computation <- function()
{
	SAS.res <- c(1.431413, 0.036431, 0.036123, 0.025228, 0.061443)
	data(realData)
	dat1   <- realData[realData$PID==1,]
	fit.qf <- anovaVCA(y~lot/calibration/day/run, dat1, SSQ.method="qf")
	fit.ws <- anovaVCA(y~lot/calibration/day/run, dat1, SSQ.method="sweep")
	R.res.qf <- round(as.numeric(fit.qf$aov.tab[-1, "VC"]), 6)
	R.res.sw <- round(as.numeric(fit.qf$aov.tab[-1, "VC"]), 6)
	checkEquals(R.res.qf, SAS.res)
	checkEquals(R.res.sw, SAS.res)
	checkEquals(R.res.qf, R.res.sw)
}


# check numerical equivalence of sweep-based and quadratic form-based SSQ computation to each other

TF028.check.sweep_vs_qf_ssq_computation.balanced <- function()
{
	fit.qf <- anovaVCA(y~(sample+device+lot)/day/run, datS1to2, SSQ.method="qf")
	fit.sw <- anovaVCA(y~(sample+device+lot)/day/run, datS1to2, SSQ.method="sweep")
	
	checkEquals(c(fit.qf$aov.tab), c(fit.sw$aov.tab))
}

# check numerical equivalence of sweep-based and quadratic form-based SSQ computation to each other

TF029.check.sweep_vs_qf_ssq_computation.unbalanced <- function()
{
	fit.qf <- anovaVCA(y~(sample+device+lot)/day/run, datS1to2.ub, SSQ.method="qf")
	fit.sw <- anovaVCA(y~(sample+device+lot)/day/run, datS1to2.ub, SSQ.method="sweep")
	
	checkEquals(c(fit.qf$aov.tab), c(fit.sw$aov.tab))
}

# check numerical equivalence of sweep-based and quadratic form-based SSQ computation to each other

TF030.check.sweep_vs_qf_ssq_computation.balanced <- function()
{
	fit.qf <- anovaVCA(y~(sample+device+lot)/day/run, datS4to5, SSQ.method="qf")
	fit.sw <- anovaVCA(y~(sample+device+lot)/day/run, datS4to5, SSQ.method="sweep")
	
	checkEquals(c(fit.qf$aov.tab), c(fit.sw$aov.tab))
}

# check numerical equivalence of sweep-based and quadratic form-based SSQ computation to each other

TF031.check.sweep_vs_qf_ssq_computation.unbalanced <- function()
{
	fit.qf <- anovaVCA(y~(sample+device+lot)/day/run, datS4to5.ub, SSQ.method="qf")
	fit.sw <- anovaVCA(y~(sample+device+lot)/day/run, datS4to5.ub, SSQ.method="sweep")
	
	checkEquals(c(fit.qf$aov.tab), c(fit.sw$aov.tab))
}


# check validity of C-implementation "getAmatBat" in R-function getCmatrix for dense
# model matrices. Reference results from SAS PROC MIXED method=type1.

TF032.check.getAmatBmat.dense.C <- function()
{
	tmp.dat <- data.frame(	y=c(1.91935483870968, 1.88709677419355, 1.9741935483871, 1.88387096774194),
			v=c(1,1,2,3))
	
	fit <- anovaVCA(y~v, tmp.dat)
	
	checkEquals(round(as.numeric(fit$aov.tab[-1, "VC"]), 6), c(0.001482, 0.000520))
	checkEquals(round(as.numeric(fit$aov.tab[-1, "MS"]), 6), c(0.002373, 0.000520))
}


TF033.check.equality.HugeData <- function()
{
	memory.limit(7500)
	data(HugeData)
	try(fit <- anovaVCA(y~VC1/VC2, HugeData), silent=TRUE)
	if(class(fit) == "try-error")
	{
		cat("\n\n########################################################################\n")
		cat(  "\n'TF033.check.equality.HugeData' cannot be run due to memory limitations!")
		cat(  "\n########################################################################\n\n")	
	}
	else
		checkEquals(round(as.numeric(fit$aov.tab[-1, "VC"]), 4), c(5173.4411, 12268.7307, 13197.6746))
}






