# TODO: unit-test functions for function 'anovaMM'
# 
# Author: schueta6
###############################################################################


cat("\n\n**************************************************************************")
cat("\nVariance Component Analysis (VCA) - test cases defined in runit.anovaMM.R.")
cat("\n**************************************************************************\n\n")

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

data(VCAdata1)

datS2 <- VCAdata1[VCAdata1$sample==2, ]

# test covariance parameter estimates, fixed effects, and covariance matrix of VCs

TF001.anovaMM.balanced1 <- function()
{
	
	fit   <- anovaMM(y~(lot+device)/(day)/(run), datS2)
	
	checkEquals(as.numeric(round(fit$aov.tab[-1, "VC"], c(4,5,5))), c(0.3147, 0.05382, 0.04409))	# VCs
	
	checkEquals(as.numeric(round(fixef(fit)[,1], 4)), c(23.2393, 0.3655, -0.5570, 0, 1.1730, 0.5925, 0))
	
}


# checks whether both function yield identical results when fitting a random model, 
# once setting negative variance estimates equal to zero, and once allowing negative 
# variance estimates

TF002.anovaMM.NegVC.balanced <- function()
{
	data(dataEP05A2_1)
	
	test.dat <- dataEP05A2_1
	set.seed(123)
	test.dat$y <- test.dat$y + rnorm(40,,2.5)						# add something which yields negative estimates
	
	resRM1 <- anovaVCA(y~day/run, test.dat, NegVC=FALSE)			# constrain VC to 0
	resMM1 <- anovaMM(y~(day)/(run), test.dat, NegVC=FALSE)
	
	checkEquals(resRM1$aov.tab, resMM1$aov.tab)
	
	resRM2 <- anovaVCA(y~day/run, test.dat, NegVC=TRUE)			# allowing negative VCs
	resMM2 <- anovaMM(y~(day)/(run), test.dat, NegVC=TRUE)
	
	checkEquals(resRM1$aov.tab, resMM1$aov.tab)
}


# checks whether both function yield identical results when fitting a random model, 
# once setting negative variance estimates equal to zero, and once allowing negative 
# variance estimates

TF003.anovaMM.NegVC.unbalanced <- function()
{
	data(dataEP05A2_1)
	
	test.dat <- dataEP05A2_1
	test.dat <- test.dat[-c(4, 12:14, 31, 55, 67:70),]
	set.seed(123)
	test.dat$y <- test.dat$y + rnorm(70,,2.5)
	
	resRM1 <- anovaVCA(y~day/run, test.dat, NegVC=FALSE)			# constrain VC to 0
	resMM1 <- anovaMM(y~(day)/(run), test.dat, NegVC=FALSE)
	
	checkEquals(resRM1$aov.tab, resMM1$aov.tab)
	
	resRM2 <- anovaVCA(y~day/run, test.dat, NegVC=TRUE)				# allowing negative VCs
	resMM2 <- anovaMM(y~(day)/(run), test.dat, NegVC=TRUE)
	
	checkEquals(resRM1$aov.tab, resMM1$aov.tab)
}


# test againts SAS PROC MIXED:
#
# proc mixed data=ep5_2 method=type1 asycov cl=wald;
#   class day run;
#   model y = day;
#   random day*run/solution;
# run;

TF004.anovaMM.ProcMixed.balanced1 <- function()
{
	data(dataEP05A2_2)
	res <- anovaMM(y~day/(run), dataEP05A2_2)
	
	checkEquals(round(as.numeric(res$aov.tab[-1, "VC" ]),4), c(2.8261, 3.7203))			# VCs
	checkEquals(round(as.numeric(res$FixedEffects[c("day10", "day15", "day20"),]), 4), c(1.8259, 4.4052, 0))		# some fixed effets
}


# test againts SAS PROC MIXED:
#
# proc mixed data=ep5_2 method=type1 asycov cl=wald;
#   class day run;
#   model y = day;
#   random day*run/solution;
# run;

TF005.anovaMM.ProcMixed.unbalanced1 <- function()
{
	data(dataEP05A2_2)
	res <- anovaMM(y~day/(run), dataEP05A2_2[-c(11,12,23,32,40,41,42),])
	
	checkEquals(round(as.numeric(res$aov.tab[-1, "VC" ]),4), c(2.2966, 3.7960))			# VCs
	checkEquals(round(as.numeric(res$FixedEffects[c("day8", "day11", "day14"),]), 4), c(-0.6839, 4.2581, 3.0765))		# some fixed effets
}



TF006.anovaDF.balanced <- function()
{		
	data(Orthodont)
	Ortho <- Orthodont
	Ortho$age2 <- Ortho$age - 11
	
	aov.fit <- anova(lm(distance~Sex*age2+(Subject)*age2, Ortho))
	mm.fit  <- anovaMM( distance~Sex*age2+(Subject)*age2, Ortho)
	
	rn <- rownames(aov.fit)
	rn[length(rn)] <- "error"
	checkEquals(as.numeric(mm.fit$aov.org[rn, "DF"]), as.numeric(aov.fit[,"Df"]))
}


TF007.anovaDF.unbalanced <- function()
{
	data(Orthodont)
	Ortho <- Orthodont
	Ortho$age2 <- Ortho$age - 11
	Ortho.ub <- Ortho[-c(3, 5, 8, 25, 29, 64, 79, 82, 102), ]				# introduce unbalancedness
	Ortho.ub$Subject <- factor(as.character(Ortho.ub$Subject))
	
	aov.fit <- anova(lm(distance~Sex*age2+(Subject)*age2, Ortho.ub))
	mm.fit  <- anovaMM( distance~Sex*age2+(Subject)*age2, Ortho.ub)
	
	rn <- rownames(aov.fit)
	rn[length(rn)] <- "error"
	
	checkEquals(as.numeric(mm.fit$aov.org[rn, "DF"]), as.numeric(aov.fit[,"Df"]))
}

TF008.anovaDF.balanced <- function()
{		
	data(Orthodont)
	Ortho <- Orthodont
	Ortho$age2 <- Ortho$age - 11
	
	aov.fit <- anova(lm(distance~Sex*age2+(Subject)*age2-1, Ortho))
	mm.fit  <- anovaMM( distance~Sex*age2+(Subject)*age2-1, Ortho)
	
	rn <- rownames(aov.fit)
	rn[length(rn)] <- "error"
	
	checkEquals(as.numeric(mm.fit$aov.org[rn, "DF"]), as.numeric(aov.fit[,"Df"]))
}

TF009.anovaDF.unbalanced <- function()
{
	data(Orthodont)
	Ortho <- Orthodont
	Ortho$age2 <- Ortho$age - 11
	Ortho.ub <- Ortho[-c(3, 5, 8, 25, 29, 64, 79, 82, 102), ]				# introduce unbalancedness
	Ortho.ub$Subject <- factor(as.character(Ortho.ub$Subject))
	
	aov.fit <- anova(lm(distance~Sex*age2+(Subject)*age2-1, Ortho.ub))
	mm.fit  <- anovaMM( distance~Sex*age2+(Subject)*age2-1, Ortho.ub)
	
	rn <- rownames(aov.fit)
	rn[length(rn)] <- "error"
	
	checkEquals(as.numeric(mm.fit$aov.org[rn, "DF"]), as.numeric(aov.fit[,"Df"]))
}

TF010.anovaMM.exception_handling <- function()
{
	checkException(anovaMM())                                                               # no input at all 
	checkException(anovaMM(Data=1))                                                                                               
	checkException(anovaMM(Data=data.frame()))                 
	checkException(anovaMM(Data=data.frame(y=1:10)))
	checkException(anovaMM(z~day/run, Data=data.frame(y=1:10)))
	checkException(anovaMM(y~day/run, Data=data.frame(y=1:10, day=1:10)))
}

TF011.anovaMM.SD_results <- function()
{ 
	res <- anovaMM(y~(day)/run, Data=dataEP05A2_1[-c(11, 12, 17, 37, 45, 56, 57, 68),])              
	checkEquals(round(sqrt(as.numeric(res$aov.tab[,"VC"])),6), round(as.numeric(res$aov.tab[,"SD"]), 6))     
	
	res <- anovaMM(y~(day)/run, Data=dataEP05A2_2[-c(2, 12, 22, 23, 24, 55, 56, 71),])              
	checkEquals(round(sqrt(as.numeric(res$aov.tab[,"VC"])),6), round(as.numeric(res$aov.tab[,"SD"]), 6)) 
	
	res <- anovaMM(y~(day)/run, Data=dataEP05A2_3[-c(1,6,7,36,61:65),])              
	checkEquals(round(sqrt(as.numeric(res$aov.tab[,"VC"])),6), round(as.numeric(res$aov.tab[,"SD"]), 6)) 
}


# check whether the reported CV-values are correctly computed by comparing 100 times square-root values 
# of the VC-value devided by the mean to values of column "CV[%]"

TF012.anovaMM.CV_perc_results <- function()
{ 
	res <- anovaMM(y~(day)/run, Data=dataEP05A2_1[-c(11, 12, 17, 37, 45, 56, 57, 68),])              
	checkEquals(round(sqrt(as.numeric(res$aov.tab[,"VC"]))*100/res$Mean,6), round(as.numeric(res$aov.tab[,"CV[%]"]), 6))     
	
	res <- anovaMM(y~(day)/run, Data=dataEP05A2_2[-c(2, 12, 22, 23, 24, 55, 56, 71),])              
	checkEquals(round(sqrt(as.numeric(res$aov.tab[,"VC"]))*100/res$Mean,6), round(as.numeric(res$aov.tab[,"CV[%]"]), 6)) 
	
	res <- anovaMM(y~(day)/run, Data=dataEP05A2_3[-c(1,6,7,36,61:65),])              
	checkEquals(round(sqrt(as.numeric(res$aov.tab[,"VC"]))*100/res$Mean,6), round(as.numeric(res$aov.tab[,"CV[%]"]), 6)) 
}

TF013.anovaMM.Percent_Total_results <- function()
{ 
	res <- anovaMM(y~(day)/run, Data=dataEP05A2_1)              
	checkEquals(round(as.numeric(res$aov.tab[,"VC"])*100/sum(as.numeric(res$aov.tab[-1, "VC"])),6), round(as.numeric(res$aov.tab[,"%Total"]), 6))           # exclude total variance in the sum  
	
	res <- anovaMM(y~(day)/run, Data=dataEP05A2_2)              
	checkEquals(round(as.numeric(res$aov.tab[,"VC"])*100/sum(as.numeric(res$aov.tab[-1, "VC"])),6), round(as.numeric(res$aov.tab[,"%Total"]), 6)) 
	
	res <- anovaMM(y~(day)/run, Data=dataEP05A2_3)              
	checkEquals(round(as.numeric(res$aov.tab[,"VC"])*100/sum(as.numeric(res$aov.tab[-1, "VC"])),6), round(as.numeric(res$aov.tab[,"%Total"]), 6)) 
}


TF014.anovaMM.missing_value_handling <- function()
{
	dat0 <- dataEP05A2_3
	
	dat0[c(5,15,25),    "y"] <- NA                      # generated missing data
	dat0[c(3,17,58),  "day"] <- NA
	dat0[c(51,70,77), "run"] <- NA
	
	datNoNA <- na.omit(dat0)
	
	res1 <- anovaMM(y~day/(run), dat0)
	res2 <- anovaMM(y~day/(run), datNoNA)
	
	checkEquals(as.matrix(res1$aov.tab), as.matrix(res2$aov.tab))
	checkEquals(res1$Nrm,   9)
	checkEquals(res1$Nobs, 71)
	checkEquals(res2$Nobs, 71)
}

# testcase checks results against SAS PROC MIXED results with ANOVA Type-1 estimation
# of random effects. Dataset is taken from R-package lme4

TF015.anovaMM.sleepstudy <- function()
{
	data(sleepstudy)
	fit.mm <- anovaMM(Reaction~Days*(Subject), sleepstudy)
	
	checkEquals(round(as.numeric(fit.mm$aov.tab[-1,"VC"]), 5), c(698.52894, 35.07166, 654.94103))
}


# testcase which should generate a warning due to numerical instabilities
# using function 'chol2inv' for obtaining a matrix inverse.
# This unit-test ensures that an exceptions is thrown and that it is handled
# correctly yielding correct results. Reference results are generated with 
# SAS PROC MIXED mehtod=type1.

TF016.anovaMM.error.chol2inv <- function()
{
	data(chol2invData)
	fit <- anovaMM(value~ID+(Site), chol2invData)
	checkEquals(fit$VCoriginal, c(-0.0237369, 1.35363935))
}
