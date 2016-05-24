# Test cases for function 'VCAinference'
#
# Author: André Schützenmeister
#
# Note:  Confidence Intervals for all variance components but total are compared 
#        to results of SAS 9.2 PROC MIXED (method=Type1, cl options, 4 decimal digits).
#        CI-limits for total VC are compared to Roche Diagnostics SAS Intranet Module VCA
#        results (if included in the testcase).used within the Satterthwaite approximation
#        For balanced data the Fisher-Information matrix used as approximation of the
#        covariance-matrix of VCs is identical to the implemented covariance matrix of
#        VCs (see 'anovaVCA', 'anovaVCAnub', 'VCAinference' for details).
#
########################################################################################################################

cat("\n\n***************************************************************************")
cat("\nVariance Component Analysis (VCA) - test cases for function 'VCAinference'.")
cat("\n***************************************************************************\n\n")

# load data and generate datasets where negative VC estimates occur

data(dataEP05A2_1)
data1 <- dataEP05A2_1
set.seed(1)
data1$y <- data1$y + rnorm(80,,5)

data(dataEP05A2_2)
data2 <- dataEP05A2_2
set.seed(2)
data2$y <- data2$y + rnorm(80,,8)

data(dataEP05A2_3)
data3 <- dataEP05A2_3
set.seed(3)
data3$y <- data3$y + rnorm(80,,10)

data(dataEP05A3_MS_1)
data(dataEP05A3_MS_2)
data(dataEP05A3_MS_3)


# check EP05-A2 20/2/2 output for balanced data (total variance, error variance)

TF001.VCAinference.constrained_CIs.balanced <- function()
{
	fit1 <- anovaVCA(y~day/run, Data=dataEP05A2_1, NegVC=TRUE, MME=TRUE)
	INF1 <- VCAinference(fit1, VarVC=TRUE, excludeNeg=FALSE, constrainCI=TRUE) 
	
	checkEquals(round(as.numeric(INF1$ConfInt$VC$TwoSided[, "LCL"]), 4), c(2.4156,        0,      0, 1.3167))            # CI VC two-sided lower limits
	checkEquals(round(as.numeric(INF1$ConfInt$VC$TwoSided[, "UCL"]), 4), c(4.7767,  1.0322,  2.8455, 3.1980))            # CI VC two-sided upper limits
	
	fit2 <- anovaVCA(y~day/run, Data=dataEP05A2_2, NegVC=TRUE, MME=TRUE)
	INF2 <- VCAinference(fit2, VarVC=TRUE, excludeNeg=FALSE, constrainCI=TRUE)
	
	checkEquals(round(as.numeric(INF2$ConfInt$VC$TwoSided[, "LCL"]), 4), c(5.9669,       0,       0, 2.5077))            # CI VC two-sided lower limits
	checkEquals(round(as.numeric(INF2$ConfInt$VC$TwoSided[, "UCL"]), 4), c(12.7046, 4.8921,  5.8428, 6.0906))            # CI VC two-sided upper limits
	
	fit3 <- anovaVCA(y~day/run, Data=dataEP05A2_3, NegVC=TRUE, MME=TRUE)
	INF3 <- VCAinference(fit3, VarVC=TRUE, excludeNeg=FALSE, constrainCI=TRUE)
	
	checkEquals(round(as.numeric(INF3$ConfInt$VC$TwoSided[, "LCL"]), 4), c(24.8957,  0,       0,       10.9764))          # CI VC two-sided lower limits
	checkEquals(round(as.numeric(INF3$ConfInt$VC$TwoSided[, "UCL"]), 4), c(54.8935,  25.6818, 17.0897, 26.6590))         # CI VC two-sided upper limits
}    

TF002.VCAinference.unconstrained_CIs.balanced <- function()
{
	fit1 <- anovaVCA(y~day/run, Data=dataEP05A2_1, NegVC=TRUE, MME=TRUE)
	INF1 <- VCAinference(fit1, VarVC=TRUE, excludeNeg=FALSE, constrainCI=FALSE) 
	
	checkEquals(round(as.numeric(INF1$ConfInt$VC$TwoSided[, "LCL"]), 4), c(2.4156, -1.0300, -0.1565, 1.3167))            # CI VC two-sided lower limits
	checkEquals(round(as.numeric(INF1$ConfInt$VC$TwoSided[, "UCL"]), 4), c(4.7767,  1.0322,  2.8455, 3.1980))            # CI VC two-sided upper limits
	
	fit2 <- anovaVCA(y~day/run, Data=dataEP05A2_2, NegVC=TRUE, MME=TRUE)
	INF2 <- VCAinference(fit2, VarVC=TRUE, excludeNeg=FALSE, constrainCI=FALSE)
	
	checkEquals(round(as.numeric(INF2$ConfInt$VC$TwoSided[, "LCL"]), 4), c(5.9669, -1.1845, -0.1907, 2.5077))            # CI VC two-sided lower limits
	checkEquals(round(as.numeric(INF2$ConfInt$VC$TwoSided[, "UCL"]), 4), c(12.7046, 4.8921,  5.8428, 6.0906))            # CI VC two-sided upper limits
	
	fit3 <- anovaVCA(y~day/run, Data=dataEP05A2_3, NegVC=TRUE, MME=TRUE)
	INF3 <- VCAinference(fit3, VarVC=TRUE, excludeNeg=FALSE, constrainCI=FALSE)
	
	checkEquals(round(as.numeric(INF3$ConfInt$VC$TwoSided[, "LCL"]), 4), c(24.8957, -1.2196, -3.0273, 10.9764))          # CI VC two-sided lower limits
	checkEquals(round(as.numeric(INF3$ConfInt$VC$TwoSided[, "UCL"]), 4), c(54.8935,  25.6818, 17.0897, 26.6590))         # CI VC two-sided upper limits
}


TF003.VCAinference.negative_VC.balanced <- function()
{
	
	# constrainCI has to be TRUE whenever any VC-estimates were set to 0
	
	fit1 <- anovaVCA(y~day/run, data1, MME=TRUE)
	INF1 <- VCAinference(fit1, VarVC=TRUE, constrainCI=FALSE)                                                 
	checkEquals(round(as.numeric(INF1$ConfInt$VC$TwoSided[, "LCL"]), 4), c( 24.2184, NA, 0.0000, 14.4422))
	
	# if CIs for originally negative VC estimates should not be excluded, always constrain them to zero if VC estimates were constrained 
	
	fit2 <- anovaVCA(y~day/run, data1, MME=TRUE)
	INF2 <- VCAinference(fit2, VarVC=TRUE, excludeNeg=FALSE)                                                  
	checkEquals(round(as.numeric(INF2$ConfInt$VC$TwoSided[, "LCL"]), 4), c(24.2184, 0.0000, 0.0000, 14.4422))
	
	# CIs of results with negative VC estimates, which are reported unconstrained (NegVC=TRUE in anovaVCA), cannot be constrained
	
	# they can only be excluded (Default)
	
	fit3 <- anovaVCA(y~day/run, data1, NegVC=TRUE, MME=TRUE)
	INF3 <- VCAinference(fit3, VarVC=TRUE, excludeNeg=TRUE)                                                  
	checkEquals(round(as.numeric(INF3$ConfInt$VC$TwoSided[, "LCL"]), 4), c(19.2231, NA, -3.0623, 14.4422))
	
	# or not 
	
	fit4 <- anovaVCA(y~day/run, data1, NegVC=TRUE, MME=TRUE)
	INF4 <- VCAinference(fit4, VarVC=TRUE, excludeNeg=FALSE)                                                  
	checkEquals(round(as.numeric(INF4$ConfInt$VC$TwoSided[, "LCL"]), 4), c(19.2231, -14.1191, -3.0623, 14.4422))
}

# check EP05-A3 3/5/5 Multi-Site output for balanced data (total variance, error variance)

TF004.VCAinference.VC_CI.balanced <- function()
{
	fit1 <- anovaVCA(y~site/day, Data=dataEP05A3_MS_1, NegVC=TRUE, MME=TRUE)
	INF1 <- VCAinference(fit1, VarVC=TRUE, excludeNeg=FALSE, constrainCI=FALSE)
	
	checkEquals(round(as.numeric(INF1$ConfInt$VC$TwoSided[, "LCL"]), 4), c(2.0437, -1.2939, -0.3148, 1.6365))            # CI VC two-sided lower limits
	checkEquals(round(as.numeric(INF1$ConfInt$VC$TwoSided[, "UCL"]), 4), c(8.1959,  3.3287,  1.0065, 3.3674))            # CI VC two-sided upper limits
	
	fit2 <- anovaVCA(y~site/day, Data=dataEP05A3_MS_2, NegVC=TRUE, MME=TRUE)
	INF2 <- VCAinference(fit2, VarVC=TRUE, excludeNeg=FALSE, constrainCI=FALSE)                                         # total not tested
	
	checkEquals(round(as.numeric(INF2$ConfInt$VC$TwoSided[-1, "LCL"]), 4), c(-1.6806, -0.4157, 2.6842))                  # CI VC two-sided lower limits
	checkEquals(round(as.numeric(INF2$ConfInt$VC$TwoSided[-1, "UCL"]), 4), c( 3.7022,  2.4723, 5.5231))                  # CI VC two-sided upper limits
	
	fit3 <- anovaVCA(y~site/day, Data=dataEP05A3_MS_3, NegVC=TRUE, MME=TRUE)
	INF3 <- VCAinference(fit3, VarVC=TRUE, excludeNeg=FALSE, constrainCI=FALSE)                                         # total not tested
	
	checkEquals(round(as.numeric(INF3$ConfInt$VC$TwoSided[-1, "LCL"]), 4), c(-20.5980, -1.4056, 10.8222))                # CI VC two-sided lower limits
	checkEquals(round(as.numeric(INF3$ConfInt$VC$TwoSided[-1, "UCL"]), 4), c( 56.5796, 12.2530, 22.2685))                # CI VC two-sided upper limits
}



# test one-sided CIs against values obtained via Excel-based PrecPerf Reports

TF005.VCAinference.WinCAEv_PrecPerf.SD_CI <- function()
{
	data_1 <- data.frame(day=c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 20, 20, 20, 20, 21, 21, 21, 21),
			run=c(1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2),
			y=c(106.2, 106.6, 106.2, 106.1, 106.8, 106.9, 107.1, 107, 107.5, 107.5, 107.2, 107, 107.8, 107.5, 106.9, 106.8, 107.1, 106.9, 106.1, 105.7, 106.5, 106.3, 105.8, 105.8, 106, 106, 105.4, 105.4, 105.4, 105.3, 104.6, 104.5, 105.2, 105.3, 104.9, 104.8, 105.1, 105.1, 104.3, 104.3, 104.5, 104.7, 104.1, 104, 104.1, 104.1, 103.4, 103.3, 103.3, 103.2, 102.6, 102.6, 103, 103.1, 102.5, 102.6, 102.9, 102.9, 102.3, 101.9, 102.5, 102.6, 102.2, 102.2, 102.5, 102.5, 102.3, 102.1, 105.5, 105.6, 105.3, 105.5, 106.5, 106.5, 106.5, 106.5, 108.1, 108.1, 108.3, 108.5, 104.9, 105.7, 105.2, 105.7))
	
	res1 <- anovaVCA(y~day/run, data_1)
	inf1 <- VCAinference(res1)
	
	checkEquals(round(as.numeric(inf1$ConfInt$SD$OneSided["total", 2:3]), 8), c(1.45561309, 2.43889695))
	checkEquals(round(as.numeric(inf1$ConfInt$SD$OneSided["error", 2:3]), 8), c(0.12750817, 0.18324098))
}


TF006.VCAinference.WinCAEv_PrecPerf.SD_CI <- function()
{
	data_2 <-data.frame(day=c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 20, 20, 20, 20, 21, 21, 21, 21),
			run=c(1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2),
			y=c(107.7,107.4,107.3,107.5,105.7,105.8,105,104.9,104.5,105,104.9,104.8,105.4,105.3,105.2,105.3,105.6,105.5,104.9,104.9,105,104.7,104.1,104.1,104.6,104.7,104,104,104.3,104.5,104.1,104,104.5,104.3,103.8,103.8,104.1,104.3,103.9,103.7,104.1,104.1,103.3,103.2,103.5,103.5,102.6,102.5,102.7,102.8,102,101.8,102.4,102.2,101.7,101.7,102,102.1,101.4,101.3,104.8,105,104.2,104.1,104.1,104.2,103.4,103.6,104.5,104.7,104.5,104.2,105.2,105.3,104.6,104.6,105.2,105.1,104.6,104.4,103.9,103.8,103.9,103.7))
	res2 <- anovaVCA(y~day/run, data_2)
	inf2 <- VCAinference(res2)
	
	checkEquals(round(as.numeric(inf2$ConfInt$SD$OneSided["total", 2:3]), 8), c(1.05907281, 1.74636429))
	checkEquals(round(as.numeric(inf2$ConfInt$SD$OneSided["error", 2:3]), 8), c(0.10075071, 0.14478805))
}




### check confidence intervals of some crossed-nested designs (balanced and unbalanced)

# compare numerical equivalence of confidence intervals of VCs to SAS PROC MIXED results
# using following SAS-code:
#   
#   proc mixed data=sample method=type1 asycov cl;
#       class lot device day run;
#       model y=;
#       random lot device lot*device*day lot*device*day*run;
#   run;

TF007.VCAinference.CrossedNested.VC_CI.balanced <- function()
{
	data(VCAdata1)
	sample1 <- VCAdata1[which(VCAdata1$sample==1),]
	sample1$device <- gl(3,28,252)                                      # add device variable
	set.seed(505)
	sample1$y <- sample1$y + rep(rep(rnorm(3,,.25), c(28,28,28)),3)     # add error component for device
	
	# write.table(sample1UB, file="sample1UB.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")      # export data, import to SAS and apply SAS-code given above
	
	res1 <- anovaVCA(y~lot+device+(lot:device:day)/run, sample1)
	inf1 <- VCAinference(res1, VarVC=TRUE, constrainCI=FALSE, ci.method="sas")
	CIs  <- inf1$ConfInt$VC$TwoSided
	CIs  <- as.matrix(CIs[,-1])
	colnames(CIs) <- NULL
	
	checkEquals(round(CIs["lot",], 				 5), c(-0.01847, 0.04951))                        # round to precision of SAS PROC MIXED output
	checkEquals(round(CIs["device",],				 5), c(-0.06322, 0.1875))
	checkEquals(round(CIs["lot:device:day",], 	 5), c(-0.00399, 0.02910))
	checkEquals(round(CIs["lot:device:day:run",], 5), c( 0.03282, 0.06865))
	checkEquals(round(CIs["error",], 				 6), c( 0.000913, 0.001500))
	
}


# unbalanced data is obtained using following strategy:
#
#   data sample1UB;
#       set sample1;
#       obs = _N_;
#   	if obs in (11,25,26,41,50) then delete;
#   run;
#
# and then running the PROC MIXED code shown above.

TF008.VCAinference.CrossedNested.VC_CI.unbalanced <- function()
{
	data(VCAdata1)
	sample1 <- VCAdata1[which(VCAdata1$sample==1),]
	sample1$device <- gl(3,28,252)                                      # add device variable
	set.seed(505)
	sample1$y <- sample1$y + rep(rep(rnorm(3,,.25), c(28,28,28)),3)     # add error component for device
	sample1   <- sample1[-c(11,25,26,41,50),]
	
	# write.table(sample1UB, file="sample1UB.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")      # export data, import to SAS and apply SAS-code given above
	
	res1 <- anovaVCA(y~lot+device+(lot:device:day)/run, sample1)
	inf1 <- VCAinference(res1, VarVC=TRUE, constrainCI=FALSE, ci.method="sas")
	CIs  <- inf1$ConfInt$VC$TwoSided
	CIs  <- as.matrix(CIs[,-1])
	colnames(CIs) <- NULL
	
	cat("\n\nSAS PROC MIXED uses the inverse of the Fisher-Information Matrix as approximation of Covariance-Matrix of Variance Components.")
	cat("\nThis is equal to the one obtained via ANOVA-approach stated in \"Variance Components\" (Searle et al. 1991) in case")
	cat("\nof balanced designs, otherwise Var(VC) of both results may differ:\n\n")
	
	CImat <- matrix(NA, ncol=9, nrow=5)
	colnames(CImat) <- c("SAS_Est", "R_Est", "SAS_LCL", "R_LCL", "SAS_UCL", "R_UCL", "Est_Diff", "LCL_Diff", "UCL_Diff")
	rownames(CImat) <- rownames(CIs)[2:6]
	CImat[,1] <- c(0.01377, 0.06190, 0.01235, 0.05125, 0.001180)
	CImat[,2] <- round(res1$aov.tab[-1, "VC"], c(5,5,5,5,6))
	CImat[,3] <- c(-0.01375, -0.06234, -0.00418, 0.03311, 0.000932)
	CImat[,4] <- round(CIs[2:6, 1], c(5,5,5,5,6))
	CImat[,5] <- c(0.04128, 0.1861, 0.02887, 0.06940, 0.001543)
	CImat[,6] <- round(CIs[2:6, 2], c(5,4,5,5,6))
	CImat[,7] <- CImat[,"SAS_Est"] - CImat[,"R_Est"]
	CImat[,8] <- CImat[,"SAS_LCL"] - CImat[,"R_LCL"]
	CImat[,9] <- CImat[,"SAS_UCL"] - CImat[,"R_UCL"]
	
	checkEquals(round(CIs["error",], 6), c( 0.000932, 0.001543))
	
	print(CImat)
}



TF009.VCAinference.CrossedNested.VC_CI.balanced <- function()
{
	data(VCAdata1)
	sample2 <- VCAdata1[which(VCAdata1$sample==2),]
	sample2$device <- gl(3,28,252)                                      # add device variable
	set.seed(511)
	sample2$y <- sample2$y + rep(rep(rnorm(3,,.25), c(28,28,28)),3)     # add error component for device
	
	# write.table(sample1UB, file="sample1UB.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")      # export data, import to SAS and apply SAS-code given above
	
	res2 <- anovaVCA(y~lot+device+(lot:device:day)/run, sample2)
	inf2 <- VCAinference(res2, VarVC=TRUE, constrainCI=FALSE, ci.method="sas")
	CIs  <- inf2$ConfInt$VC$TwoSided
	CIs  <- as.matrix(CIs[,-1])
	colnames(CIs) <- NULL
	
	checkEquals(round(CIs["lot",], 				 4), c(-0.2240,  0.6220))                        # round to precision of SAS PROC MIXED output
	checkEquals(round(CIs["device",],				 4), c(-0.4456,  1.3055))
	checkEquals(round(CIs["lot:device:day",], 	 4), c( 0.1857,  0.4437))
	checkEquals(round(CIs["lot:device:day:run",], 5), c( 0.02677, 0.08087))
	checkEquals(round(CIs["error",], 				 5), c( 0.03495, 0.05738))	
}



TF010.VCAinference.CrossedNested.VC_CI.unbalanced <- function()
{
	data(VCAdata1)
	sample2 <- VCAdata1[which(VCAdata1$sample==2),]
	sample2$device <- gl(3,28,252)                                      # add device variable
	set.seed(511)
	sample2$y <- sample2$y + rep(rep(rnorm(3,,.25), c(28,28,28)),3)     # add error component for device
	sample2   <- sample2[-c(1,5,36,39,51),]
	
	# write.table(sample1UB, file="sample1UB.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")      # export data, import to SAS and apply SAS-code given above
	
	res2 <- anovaVCA(y~lot+device+(lot:device:day)/run, sample2)
	inf2 <- VCAinference(res2, VarVC=TRUE, constrainCI=FALSE, ci.method="sas")
	CIs  <- inf2$ConfInt$VC$TwoSided
	CIs  <- as.matrix(CIs[,-1])
	colnames(CIs) <- NULL
	
	cat("\n\nSAS PROC MIXED uses the inverse of the Fisher-Information Matrix as approximation of Covariance-Matrix of Variance Components.")
	cat("\nThis is equal to the one obtained via ANOVA-approach stated in \"Variance Components\" (Searle et al. 1991) in case")
	cat("\nof balanced designs, otherwise Var(VC) of both results may differ:\n\n")
	
	CImat <- matrix(NA, ncol=9, nrow=5)
	colnames(CImat) <- c("SAS_Est", "R_Est", "SAS_LCL", "R_LCL", "SAS_UCL", "R_UCL", "Est_Diff", "LCL_Diff", "UCL_Diff")
	rownames(CImat) <- rownames(CIs)[2:6]
	CImat[,1] <- c(0.1892, 0.4359, 0.3178, 0.05724, 0.04108)
	CImat[,2] <- round(res2$aov.tab[-1, "VC"], c(4,4,4,5,5))
	CImat[,3] <- c(-0.1909, -0.4608, 0.1860, 0.02931, 0.03241)
	CImat[,4] <- round(CIs[2:6, 1],c(4,4,4,5,5))
	CImat[,5] <- c(0.5692, 1.3326, 0.4496, 0.08517, 0.05376)
	CImat[,6] <- round(CIs[2:6, 2], c(4,4,4,5,5))
	CImat[,7] <- CImat[,"SAS_Est"] - CImat[,"R_Est"]
	CImat[,8] <- CImat[,"SAS_LCL"] - CImat[,"R_LCL"]
	CImat[,9] <- CImat[,"SAS_UCL"] - CImat[,"R_UCL"]
	
	checkEquals(round(CIs["error",], 5), c( 0.03241, 0.05376))
	
	print(CImat)
}



TF011.VCAinference.CrossedNested.VC_CI.balanced <- function()
{
	data(VCAdata1)
	sample3 <- VCAdata1[which(VCAdata1$sample==3),]
	sample3$device <- gl(3,28,252)                                      # add device variable
	set.seed(519)
	sample3$y <- sample3$y + rep(rep(rnorm(3,,.25), c(28,28,28)),3)     # add error component for device
	
	# write.table(sample1UB, file="sample1UB.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")      # export data, import to SAS and apply SAS-code given above
	
	res3 <- anovaVCA(y~lot+device+(lot:device:day)/run, sample3, NegVC=TRUE)
	inf3 <- VCAinference(res3, VarVC=TRUE, constrainCI=FALSE, excludeNeg=FALSE, ci.method="sas")
	CIs  <- inf3$ConfInt$VC$TwoSided
	CIs  <- as.matrix(CIs[,-1])
	colnames(CIs) <- NULL
	
	checkEquals(round(CIs["lot",], 				 c(5,6)), c(-0.00168,  0.000049))                        # round to precision of SAS PROC MIXED output
	checkEquals(round(CIs["device",],				 c(5,4)), c(-0.05672,  0.1700))
	checkEquals(round(CIs["lot:device:day",], 	 c(5,6)), c(-0.03297,  0.001022))
	checkEquals(round(CIs["lot:device:day:run",], c(5,4)), c( 0.05327,  0.1106))
	checkEquals(round(CIs["error",], 				 c(6,6)), c( 0.000263, 0.000432))	
}



TF012.VCAinference.CrossedNested.VC_CI.unbalanced <- function()
{
	data(VCAdata1)
	sample3 <- VCAdata1[which(VCAdata1$sample==3),]
	sample3$device <- gl(3,28,252)                                      # add device variable
	set.seed(519)
	sample3$y <- sample3$y + rep(rep(rnorm(3,,.25), c(28,28,28)),3)     # add error component for device
	sample3   <- sample3[-c(9,10,11,29,30),]
	
	# write.table(sample1UB, file="sample1UB.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")      # export data, import to SAS and apply SAS-code given above
	
	res3 <- anovaVCA(y~lot+device+(lot:device:day)/run, sample3, NegVC=TRUE)
	inf3 <- VCAinference(res3, VarVC=TRUE, constrainCI=FALSE, excludeNeg=FALSE, ci.method="sas")
	CIs  <- inf3$ConfInt$VC$TwoSided
	CIs  <- as.matrix(CIs[,-1])
	colnames(CIs) <- NULL
	
	cat("\n\nSAS PROC MIXED uses the inverse of the Fisher-Information Matrix as approximation of Covariance-Matrix of Variance Components.")
	cat("\nThis is equal to the one obtained via ANOVA-approach stated in \"Variance Components\" (Searle et al. 1991) in case")
	cat("\nof balanced designs, otherwise Var(VC) of both results may differ:\n\n")
	
	CImat <- matrix(NA, ncol=9, nrow=5)
	colnames(CImat) <- c("SAS_Est", "R_Est", "SAS_LCL", "R_LCL", "SAS_UCL", "R_UCL", "Est_Diff", "LCL_Diff", "UCL_Diff")
	rownames(CImat) <- rownames(CIs)[2:6]
	CImat[,1] <- c(-0.00085, 0.05804, -0.01741, 0.08400, 0.000337)
	CImat[,2] <- round(res3$aov.tab[-1, "VC"], c(5,5,5,5,6))
	CImat[,3] <- c(-0.00154, -0.05762, -0.03482, 0.05394, 0.000266)
	CImat[,4] <- round(CIs[2:6, 1],c(5,5,5,5,6))
	CImat[,5] <- c(-0.00017, 0.1737, -0.00000325, 0.1141, 0.000440)
	CImat[,6] <- round(CIs[2:6, 2], c(5,4,8,4,6))
	CImat[,7] <- CImat[,"SAS_Est"] - CImat[,"R_Est"]
	CImat[,8] <- CImat[,"SAS_LCL"] - CImat[,"R_LCL"]
	CImat[,9] <- CImat[,"SAS_UCL"] - CImat[,"R_UCL"]
	
	checkEquals(round(CIs["error",], 6), c( 0.000266, 0.000440))
	
	print(CImat)
}


# Test Satterthwaite methodology for confidence intervals of all variance components. This is different from the SAS PROC MIXED
# methodology, whenever Type1 estimation method is chosen, which is tested in all other test functions above. Here, all VC CIs
# are constructed using the Chi-Squared distribution and DFs are approximated according to Satterthwaite.
#
# For balanced designs, CIs are equal of those constructed by SAS PROC MIXED with method=REML. In this setting, CIs are constructed
# using the "Satterthwaite methodology".

TF013.anovaVCA.Satt_methodology_for_CI.balanced <- function()
{
	data(dataEP05A2_2)														
	res <- anovaVCA(y~day/run, dataEP05A2_2)
	inf <- VCAinference(res, VarVC=TRUE, ci.method="satterthwaite")
	
	checkEquals(round(as.numeric(inf$ConfInt$VC$TwoSided[-1, "LCL"]), 4), c(0.5835, 1.2203, 2.5077))
	checkEquals(round(as.numeric(inf$ConfInt$VC$TwoSided[-1, "UCL"]), 4), c(28.5302, 12.1410, 6.0906))
}

TF014.anovaVCA.Satt_methodology_for_CI.balanced <- function()
{
	data(dataEP05A2_3)
	res <- anovaVCA(y~day/run, dataEP05A2_3)
	inf <- VCAinference(res, VarVC=TRUE, ci.method="satterthwaite")
	
	checkEquals(round(as.numeric(inf$ConfInt$VC$TwoSided[-1, "LCL"]), 4), c(5.1779, 2.4638, 10.9764))
	checkEquals(round(as.numeric(inf$ConfInt$VC$TwoSided[-1, "UCL"]), 4), c(55.8104, 64.3402, 26.6590))
}

# now checking mixed models, where day is fixed

TF015.anovaMM.Satt_methodology.balanced <- function()
{
	data(dataEP05A2_2)
	res <- anovaMM(y~day/(run), dataEP05A2_2)
	inf <- VCAinference(res, VarVC=TRUE, ci.method="satterthwaite")
	
	checkEquals(round(as.numeric(inf$ConfInt$VC$TwoSided[-1, "LCL"]), 4), c(1.2203, 2.5077))
	checkEquals(round(as.numeric(inf$ConfInt$VC$TwoSided[-1, "UCL"]), 4), c(12.1410, 6.0906))
}

TF016.anovaMM.Satt_methodology.balanced <- function()
{
	data(dataEP05A2_3)
	res <- anovaMM(y~day/(run), dataEP05A2_3)
	inf <- VCAinference(res, VarVC=TRUE, ci.method="satterthwaite")
	
	checkEquals(round(as.numeric(inf$ConfInt$VC$TwoSided[-1, "LCL"]), 4), c(2.4638, 10.9764))
	checkEquals(round(as.numeric(inf$ConfInt$VC$TwoSided[-1, "UCL"]), 4), c(64.3402, 26.6590))
}


# This test checks whether point estimates exceeding the bounds of a confidence interval
# are correctly set to NA.

TF017.anovaVCA.Est_outsided_CI <- function()
{
	data(ReproData1)										# input data
	fit <- anovaVCA(value~Site/Day/Run, ReproData1)
	inf <- tryCatch(VCAinference(fit, ci.method="satterthwaite"),		# warning should be issued
			warning=function(w) w )
	checkTrue(is(inf, "warning"))
	inf 	<- VCAinference(fit, ci.method="satterthwaite")
	infVC	<- print(inf, what="VC")
	checkTrue(is.na(infVC["Site:Day:Run", "CI LCL"]))
	checkTrue(is.na(infVC["Site:Day:Run", "CI UCL"]))
	checkTrue(is.na(infVC["Site:Day:Run", "One-Sided LCL"]))
}


# check EP05-A2 20/2/2 output for balanced data (total variance, error variance)

TF018.VCAinference.constrained_CIs.balanced <- function()
{
	fit1 <- remlVCA(y~day/run, Data=dataEP05A2_1)
	INF1 <- VCAinference(fit1, VarVC=TRUE, excludeNeg=FALSE, constrainCI=TRUE) 
	
	checkEquals(round(as.numeric(INF1$ConfInt$VC$TwoSided[, "LCL"]), 4), c(2.4156,        0,      0, 1.3167))            # CI VC two-sided lower limits
	checkEquals(round(as.numeric(INF1$ConfInt$VC$TwoSided[, "UCL"]), 4), c(4.7767,  1.0322,  2.8455, 3.1980))            # CI VC two-sided upper limits
	
	fit2 <- remlVCA(y~day/run, Data=dataEP05A2_2)
	INF2 <- VCAinference(fit2, VarVC=TRUE, excludeNeg=FALSE, constrainCI=TRUE)
	
	checkEquals(round(as.numeric(INF2$ConfInt$VC$TwoSided[, "LCL"]), 4), c(5.9669,       0,       0, 2.5077))            # CI VC two-sided lower limits
	checkEquals(round(as.numeric(INF2$ConfInt$VC$TwoSided[, "UCL"]), 4), c(12.7046, 4.8921,  5.8428, 6.0906))            # CI VC two-sided upper limits
	
	fit3 <- remlVCA(y~day/run, Data=dataEP05A2_3)
	INF3 <- VCAinference(fit3, VarVC=TRUE, excludeNeg=FALSE, constrainCI=TRUE)
	
	checkEquals(round(as.numeric(INF3$ConfInt$VC$TwoSided[, "LCL"]), 4), c(24.8957,  0,       0,       10.9764))          # CI VC two-sided lower limits
	checkEquals(round(as.numeric(INF3$ConfInt$VC$TwoSided[, "UCL"]), 4), c(54.8935,  25.6818, 17.0897, 26.6590))         # CI VC two-sided upper limits
}    

TF019.VCAinference.unconstrained_CIs.balanced <- function()
{
	fit1 <- remlVCA(y~day/run, Data=dataEP05A2_1)
	INF1 <- VCAinference(fit1, excludeNeg=FALSE, constrainCI=FALSE) 
	
	checkEquals(round(as.numeric(INF1$ConfInt$VC$TwoSided[, "LCL"]), 4), c(2.4156, -1.0300, -0.1565, 1.3167))            # CI VC two-sided lower limits
	checkEquals(round(as.numeric(INF1$ConfInt$VC$TwoSided[, "UCL"]), 4), c(4.7767,  1.0322,  2.8455, 3.1980))            # CI VC two-sided upper limits
	
	fit3 <- remlVCA(y~day/run, Data=dataEP05A2_3)
	INF3 <- VCAinference(fit3, VarVC=TRUE, excludeNeg=FALSE, constrainCI=FALSE)
	
	checkEquals(round(as.numeric(INF3$ConfInt$VC$TwoSided[, "LCL"]), 4), c(24.8957, -1.2196, -3.0273, 10.9764))          # CI VC two-sided lower limits
	checkEquals(round(as.numeric(INF3$ConfInt$VC$TwoSided[, "UCL"]), 4), c(54.8935,  25.6818, 17.0897, 26.6590))         # CI VC two-sided upper limits
}


TF020.VCAinference.balanced.REML.satterthwaite <- function()
{
	
	# constrainCI has to be TRUE whenever any VC-estimates were set to 0
	
	fit1 <- remlVCA(y~day/run, data1)
	INF1 <- VCAinference(fit1, ci.method="satt") 		# need Satterthwaite here to have SAS PROC MIXED reference results                                                
	checkEquals(round(as.numeric(INF1$ConfInt$VC$TwoSided[, "LCL"]), 2), c( 19.64, NA, 1.51, 14.44))
	checkEquals(round(as.numeric(INF1$ConfInt$VC$TwoSided[, "UCL"]), 2), c( 37.24, NA, 89.63, 35.08))
	
	fit2 <- remlVCA(y~day/run, data2)
	INF2 <- VCAinference(fit2, ci.method="satt")
	checkEquals(round(as.numeric(INF2$ConfInt$VC$TwoSided[, "LCL"]), 4), c( 63.1044 , NA, NA, 63.1044 ))
	checkEquals(round(as.numeric(INF2$ConfInt$VC$TwoSided[, "UCL"]), 4), c( 118.2015, NA, NA, 118.2015))
}

# check EP05-A3 3/5/5 Multi-Site output for balanced data (total variance, error variance)

TF021.VCAinference.VC_CI.balanced.REML <- function()
{
	fit1 <- remlVCA(y~site/day, Data=dataEP05A3_MS_1)
	INF1 <- VCAinference(fit1, VarVC=TRUE, excludeNeg=FALSE, constrainCI=FALSE)
	
	checkEquals(round(as.numeric(INF1$ConfInt$VC$TwoSided[, "LCL"]), 4), c(2.0437, -1.2939, -0.3148, 1.6365))            # CI VC two-sided lower limits
	checkEquals(round(as.numeric(INF1$ConfInt$VC$TwoSided[, "UCL"]), 4), c(8.1959,  3.3287,  1.0065, 3.3674))            # CI VC two-sided upper limits
	
	fit2 <- remlVCA(y~site/day, Data=dataEP05A3_MS_2)
	INF2 <- VCAinference(fit2, VarVC=TRUE, excludeNeg=FALSE, constrainCI=FALSE)                                         # total not tested
	
	checkEquals(round(as.numeric(INF2$ConfInt$VC$TwoSided[-1, "LCL"]), 4), c(-1.6806, -0.4157, 2.6842))                  # CI VC two-sided lower limits
	checkEquals(round(as.numeric(INF2$ConfInt$VC$TwoSided[-1, "UCL"]), 4), c( 3.7022,  2.4723, 5.5231))                  # CI VC two-sided upper limits
	
	fit3 <- remlVCA(y~site/day, Data=dataEP05A3_MS_3)
	INF3 <- VCAinference(fit3, VarVC=TRUE, excludeNeg=FALSE, constrainCI=FALSE)                                         # total not tested
	
	checkEquals(round(as.numeric(INF3$ConfInt$VC$TwoSided[-1, "LCL"]), 4), c(-20.5980, -1.4056, 10.8222))                # CI VC two-sided lower limits
	checkEquals(round(as.numeric(INF3$ConfInt$VC$TwoSided[-1, "UCL"]), 4), c( 56.5796, 12.2530, 22.2685))                # CI VC two-sided upper limits
}



# check test case collection against SAS PROC MIXED method=type1 results 

noTF <- 22			# number of next test function




########## create unit-test functions for SAS-reference data of the real world dataset

# All leading and trailing zeros are removed as well as all non-digits (decimal point, "e-" in exponential form).
# Remaining digits are counted and returned as number of signficant digits.

Nsignif <- function(x){ return( nchar(gsub("0+$", "", gsub("^0+", "", sub("e-..$", "", sub("\\.", "", x))))) )}     # determine the number of significant digits

TestPrecision <- 1e-12

data(realData)

SASrefModel1 <- read.delim("SASresultsRealData_Model1_by_PID_Lot.txt", header=TRUE)

by_PID_Lot <- realData[,c("PID", "lot")]
by_PID_Lot <- unique(by_PID_Lot)
by_PID_Lot$PID <- as.character(by_PID_Lot$PID)
by_PID_Lot$lot <- as.character(by_PID_Lot$lot)


# Function performs tests for each combination of PID and lot in the real-world dataset
# applying model1:

model1 <- y~calibration/day/run

generic.test.function <- function(Data, prec=NULL)
{
	PID <- Data$PID[1]
	Lot <- Data$lot[1]
	
	cat("\n>>> REAL-DATA MODEL 1: Check against SAS PROC MIXED VCA reference results for PID =",PID, " and Lot =", Lot )
	
	cat("\n\nEquivalence of analysis results tested for:\n")
	cat("\n\t-all variance component (VC) estimates")
	cat("\n\t-confidence interval of error VC (LCL, UCL)")
	
	if(PID == 5 && Lot == 3 || PID == 3 && Lot == 1)
	{
		cat("\n\nDue to numerical reasons of representing floating point numbers as approximations,")
		cat("\ni.e. 0.50125 represented as 0.512499999..., this test-case has larger precision value (e-03).")
		
		prec <- 1e-03
	}
	
	cat("\n\nnumerical tolerance:", prec, "(R-results element-wise rounded to the number of significant digits of SAS PROC MIXED results)\n")
	
	SASresult <- SASrefModel1[which(SASrefModel1$PID == PID & SASrefModel1$lot == Lot),]
	
	Rresult <- anovaVCA(model1, Data, NegVC=TRUE)
	RinfRes <- VCAinference(Rresult, VarVC=TRUE, constrainCI=FALSE)$ConfInt$VC$TwoSided
	
	if(as.integer(as.character(PID)) == 2 && as.integer(as.character(Lot)) == 3)
	{
		cat("\n\nUse numeric precision of 1e-03 now!\n\n")
		prec <- 1e-03
		
		cat("R-Results:\n")
		print(Rresult)
		cat("\n\nSAS-Results:\n")
		print(SASresult)
		cat("\n\n")
	}
	
	nr <- nrow(Rresult$aov.tab)
	
	# actually check SAS- and R-results
	
	for(i in 2:nr)
	{
		checkEquals(signif(Rresult$aov.tab[i, "VC"], Nsignif(SASresult[i-1, "Estimate"])), SASresult[i-1, "Estimate"],  # check equality of all VC-estimates
				tolerance=prec)
		
		if(i == 5)																								# also check CI-limits of error-VC			
		{
			checkEquals(signif(RinfRes[i, "LCL"], Nsignif(SASresult[i-1, "Lower"])), SASresult[i-1, "Lower"],
					tolerance=prec)
			checkEquals(signif(RinfRes[i, "UCL"], Nsignif(SASresult[i-1, "Upper"])), SASresult[i-1, "Upper"],
					tolerance=prec)
		}
	}
	
	cat("\n\nSAS PROC MIXED uses the inverse of the Fisher-Information Matrix as approximation of Covariance-Matrix of Variance Components.")
	cat("\nThis is equal to the one obtained via ANOVA-approach stated in \"Variance Components\" (Searle et al. 1991) in case")
	cat("\nof balanced designs, otherwise Var(VC) of both results may differ:\n\n")
	
	CImat <- matrix(NA, ncol=6, nrow=4)
	colnames(CImat) <- c("SAS_LCL", "R_LCL", "SAS_UCL", "R_UCL", "LCL_Diff", "UCL_Diff")
	rownames(CImat) <- rownames(Rresult$aov.tab)[2:5]
	CImat[,1] <- SASresult[,"Lower"]
	CImat[,3] <- SASresult[,"Upper"]
	CImat[,2] <- RinfRes[2:5,"LCL"]
	CImat[,4] <- RinfRes[2:5, "UCL"]
	CImat[,5] <- CImat[,1] - CImat[,2]
	CImat[,6] <- CImat[,3] - CImat[,4]
	
	print(CImat)
	
	cat("\n\n")
}

formatTFno <- function(x)
{
	if(x < 10)
		return(paste("00", x, sep=""))
	else if(x < 100)
		return(paste("0", x, sep=""))
	else
		return(x)
}


# Function serves as helper function to define input-parameters of 'generic.test.function' in
# separate environments.

get.test.function <- function(Data, prec=NULL)
{
	inputData <- Data
	precision <- prec
	
	f <- function(){generic.test.function(Data=inputData, prec=precision)}
	return(f)
}

if(realWorldModel1)
{
	for(i in 1:nrow(by_PID_Lot))
	{
		fname <- paste("TF", formatTFno(noTF+i), ".PID_", by_PID_Lot[i, "PID"], "_Lot_", by_PID_Lot[i, "lot"],  sep="")
		
		tmp.data <- realData[which(realData$PID == by_PID_Lot[i, "PID"] &  realData$lot == by_PID_Lot[i, "lot"]),]
		
		assign(fname, get.test.function(Data=tmp.data, prec=TestPrecision))
	}
	
	noTF <- noTF + nrow(by_PID_Lot)
}





######### Model 2 for the real-world dataset


model2 <- y~lot/calibration/day/run

SASrefModel2 <- read.delim("SASresultsRealData_Model2_by_PID.txt", header=TRUE)

PIDs <- unique(as.character(realData$PID))

generic.test.function2 <- function(Data, prec=NULL)
{
	PID <- Data$PID[1]
	
	cat("\n>>> REAL-DATA MODEL 2: Check R VCA results v.s SAS PROC MIXED reference results for PID =",PID)
	
	cat(" (tolerance-value:", prec, ")\n")
	
	SASresult <- SASrefModel2[which(SASrefModel2$PID == PID),]
	
	Rresult <- anovaVCA(model2, Data, NegVC=TRUE)
	RinfRes <- VCAinference(Rresult, VarVC=TRUE, constrainCI=FALSE)$ConfInt$VC$TwoSided
	
	nr <- nrow(Rresult$aov.tab)
	
	# actually check SAS- and R-results
	
	for(i in 2:nr)
	{
		checkEquals(signif(Rresult$aov.tab[i, "VC"], Nsignif(SASresult[i-1, "Estimate"])), SASresult[i-1, "Estimate"],  # check equality of all VC-estimates
				tolerance=prec)
		
		if(i == 6)																								# also check CI of error-VC			
		{
			checkEquals(signif(RinfRes[i, "LCL"], Nsignif(SASresult[i-1, "Lower"])), SASresult[i-1, "Lower"],
					tolerance=prec)
			checkEquals(signif(RinfRes[i, "UCL"], Nsignif(SASresult[i-1, "Upper"])), SASresult[i-1, "Upper"],
					tolerance=prec)
		}
	}
	
	cat("\n\nSAS PROC MIXED uses the inverse of the Fisher-Information Matrix as approximation of Covariance-Matrix of Variance Components.")
	cat("\nThis is equal to the one obtained via ANOVA-approach stated in \"Variance Components\" (Searle et al. 1991) in case")
	cat("\nof balanced designs, otherwise Var(VC) of both results may differ:\n\n")
	
	CImat <- matrix(NA, ncol=6, nrow=5)
	colnames(CImat) <- c("SAS_LCL", "R_LCL", "SAS_UCL", "R_UCL", "LCL_Diff", "UCL_Diff")
	rownames(CImat) <- rownames(Rresult$aov.tab)[2:6]
	CImat[,1] <- SASresult[,"Lower"]
	CImat[,3] <- SASresult[,"Upper"]
	CImat[,2] <- RinfRes[2:6,"LCL"]
	CImat[,4] <- RinfRes[2:6, "UCL"]
	CImat[,5] <- CImat[,1] - CImat[,2]
	CImat[,6] <- CImat[,3] - CImat[,4]
	
	print(CImat)
	
	cat("\n\n")
}


# Function serves as helper function to define input-parameters of 'generic.test.function' in
# separate environments.

get.test.function2 <- function(Data, prec=NULL)
{
	inputData <- Data
	precision <- prec
	
	f <- function(){generic.test.function2(Data=inputData, prec=precision)}
	return(f)
}


if(realWorldModel2)
{
	for(i in 1:length(PIDs))
	{
		fname <- paste("TF", formatTFno(noTF+i), ".Model2_PID_", PIDs[i], sep="")
		
		tmp.data <- realData[which(realData$PID == PIDs[i]),]
		
		assign(fname, get.test.function2(Data=tmp.data, prec=TestPrecision) )
	}
	
}