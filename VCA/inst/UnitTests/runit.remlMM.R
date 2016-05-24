# TODO: unit-test functions for function 'remlMM'
# 
# Author: schueta6
###############################################################################


cat("\n\n**************************************************************************")
cat("\nVariance Component Analysis (VCA) - test cases defined in runit.remlMM.R.")
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

TF010.anovaMM.exception_handling <- function()
{
	checkException(anovaMM())                                                               # no input at all 
	checkException(anovaMM(Data=1))                                                                                               
	checkException(anovaMM(Data=data.frame()))                 
	checkException(anovaMM(Data=data.frame(y=1:10)))
	checkException(anovaMM(z~day/run, Data=data.frame(y=1:10)))
	checkException(anovaMM(y~day/run, Data=data.frame(y=1:10, day=1:10)))
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


# testcase checks results against SAS PROC MIXED results with ANOVA Type-1 estimation
# of random effects. Dataset is taken from R-package lme4

TF015.anovaMM.sleepstudy <- function()
{
	data(sleepstudy)
	fit.mm <- anovaMM(Reaction~Days*(Subject), sleepstudy)
	
	checkEquals(round(as.numeric(fit.mm$aov.tab[-1,"VC"]), 5), c(698.52894, 35.07166, 654.94103))
}


