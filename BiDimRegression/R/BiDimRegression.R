BiDimRegression <-
function (coord) 
{
	# function BiDimRegression - 
	#    BiDimensional Regression (Version 1.0.6)

	# PURPOSE: calculates the bidimensional regression
	#          between two 2D configurations

	# -------------------------------------------------------- 
	# USAGE: resultingMeasures <- BiDimRegression(coord); 
	# e.g.USAGE: resultingMeasures <- BiDimRegression(NakayaData); 
	#  where: coord is a data.frame consisting of A, B, X, Y:
	#	    depV1, depV2 --> vectors containing the dependent coordinates (A=dimension1, B=dimension2)
	#         indepV1, indepV2 --> vectors containing the independent coordinates (X=dimension1, Y=dimension2)
	#
	# -------------------------------------------------------- 
	# RETURN: a data.frame "resultingMeasures" consisting of essential measures of the bidimensional regression
	#   for a) the euclidean solution (data-frame "euclidean") and 
	#	b) the affine solution (data.frame "affine"): 
	#   	r: gives the regression coefficient, analogous defined to Pearson’s r. 
	#   	rsqr: gives the squared regression coefficient.
	#       diABSqr: gives the squared distortion index di for AB; di: Following Waterman and Gordon’s (1984) extension of the bidimensional 		#		regression diprovides a measure of comparison of distortions (di=distortion index), but the range of values is 0 to 1 following 	#		Friedman and Kohler (2003)
	#   	dMaxABSqr: maximal distortion factor, needed for calculating di for AB.
	#     	diXYSqr: gives the squared distortion index di for XY.
	#   	dMaxXYSqr: maximal distortion factor, needed for calculating di for XY.
	#   	scaleFactorX: returns the scaling factor on the first dimension (1.0 means no scaling; values below 1.0 indicate a contraction; values 		#		above 1.0 indicate an expansion).
	#   	scaleFactorY: returns the scaling factor on the second dimension. 
	#   	angleDEG: returns the rotation angle (in degrees). 
	#   	shear: this return value is only active, if 'affine' is set as regressType (regression type) and gives information on the shearing of 		#		the transformed configuration. 
	#	ttestDF: degrees of freedom (DF) for the t-tests regarding the model parameters (alphas and betas)
	#   	alpha1 / alpha2: following the notation of Friedman and Kohler (2003), alpha reflects the intercept vectors of the bidimensional 		#		regression  equation shown below.
	#   	beta1 / beta2: following the notation of Friedman and Kohler (2003), beta reflects the slope vectors of the bidimensional regression 		#		equation  shown below.
	#   	[only available for affine solution] beta3 / beta4: following the notation of Friedman and Kohler (2003), beta reflects the slope 		#		vectors of the bidimensional regression equation shown below.
	#   	fValue: the bidimensional regression will be tested for significance on basis of F-distributions following the advice of Nakaya (1997).
	#   	df1: degrees of freedom of the nominator used for the F-statistics propagated by Nakaya (1997); df1 = p-2, with p is the number of 		#		elements  needed to calculate the referring model: p=4 for the Euclidean and p=6 for the affine geometry {Nakaya, 1997 #5659, 		#		Table 1}. 
	#   	df2: degrees of freedom of the denominator used for the F-statistics propagated by Nakaya (1997); df2 = 2n-p, with p is the number of 		#		elements needed to calculate the referring model (see df1) and n is the number of coordinate pairs.
	#   	pValue: the significance level based on the preceding F-statistics.
	#	dAICso: the AIC difference between the regarding bidimensional regression model and the bidimensional null model (S0) according to 		#	Nakaya (1997), formula 56
	#  and c) comparative statistics (data-frame "eucVSaff") comparing the Euclidean and the Affine solution:
	#     dAIC: Difference AIC = Akaike Information Criterion.
	#	fValue: the F-value of the regarding test on basis of a F-distribution.
	#	pValue: the significance level based on the preceding F-statistics.
	#	df1: degrees of freedom of the nominator used for the F-statistics.
	#	df2: degrees of freedom of the denominator used for the F-statistics.
	#
	# -------------------------------------------------------- 
	# REFERENCES: 
	#   Carbon, C. C. (2013). BiDimRegression: Bidimensional Regression Modeling Using R. Journal of Statistical
 	#     Software, Code Snippets, 52(1), 1-11 (URL http://www.jstatsoft.org/v52/c01/)\\
	#   Tobler, W. [R.](1965). Computation of the corresponding of geographical patterns. Papers of the Regional Science Association, 15, 131-139.
	#   Tobler, W. R. (1966). Medieval distortions: Projections of ancient maps. Annals of the Association of American Geographers, 56(2), 351-360.
	#   Tobler, W. R. (1994). Bidimensional regression. Geographical Analysis, 26(3), 187-212.
	#   Friedman, A., & Kohler, B. (2003). Bidimensional regression: Assessing the configural similarity and accuracy of cognitive maps and other 		#	two-dimensional data sets. Psychological Methods, 8(4), 468-491.
	#   Nakaya, T. (1997). Statistical inferences in bidimensional regression models. Geographical Analysis, 29(2), 169-186.
	#   Waterman, S., & Gordon, D. (1984). A quantitative-comparative approach to analysis of distortion in mental maps. Professional Geographer, 		#	36(3),326-337.
	#
	#--------------------------------------------------------  
	# AUTHORSHIP:
	# main routine written in 2005 & 2011-2013 by: 
	#       Claus-Christian Carbon, PhD
	#       Department of General Psychology and Methodology
	#       University of Bamberg
	#       Markusplatz 3
	#       D-96047 Bamberg, Bavaria, Germany
	#       ccc@experimental-psychology.com
	#
	#***********************************************************
	#
	# package:
	#     library("BiDimRegression")
	# usage:
	#     resultsBiDimRegr <- BiDimRegression(coord)
	#	print(resultsBiDimRegr)
	#	summary(resultsBiDimRegr)
	# exemplary data:
	# 	data("NakayaData", package = "BiDimRegression") # uses data from Nakaya (1997)
	#	resultsBiDimRegr <- BiDimRegression(NakayaData)
	# 	print(resultsBiDimRegr)
	#
	#	data("FriedmanKohlerData1", package = "BiDimRegression") # uses data from Friedman & Kohler (2003) ex.1
	# 	data("FriedmanKohlerData2", package = "BiDimRegression") # uses data from Friedman & Kohler (2003) ex.2
	# 	data("CarbonExample1Data", package = "BiDimRegression")  # uses data from ex.1 (Mona Lisa) by Carbon (in press)
	# 	data("NakayaData", package = "BiDimRegression")  # uses data from ex.1 (Mona Lisa) by Carbon (in press)
	# 	data("CarbonExample2Data", package = "BiDimRegression")  # uses data from ex.2 (Paris) by Carbon (in press)
	# 	data("CarbonExample3Data", package = "BiDimRegression")  # uses data from ex.3 (cognitive map) by Carbon (in press)
	#
	#***********************************************************

	# calculates the bidimensional regression according to Tobler's initial ideas
	# and on basis of Nakaya's (1997) additional inferential statistics

      	# start timer 
	timerStart <- proc.time()

	# set standard variables
	timerStart <- proc.time()
	n <- dim(coord)[1]   # number of coordinates
	vecZero <- c(rep(0, n))
      	vecOne <- c(rep(1, n))
	A <- coord$depV1
	B <- coord$depV2
	X <- coord$indepV1
	Y <- coord$indepV2

	# calculating means
	Am <- mean(A)
	Bm <- mean(B)
	Xm <- mean(X)
	Ym <- mean(Y)

	# calculating (co)variances
	X2 <- sum(X^2)
	Y2 <- sum(Y^2)
	sumX2Y2 <- sum(X^2+Y^2)
	A2 <- sum(A^2)
	B2 <- sum(B^2)
	sumA2B2 <- sum(A^2+B^2)
	varA <- (sum((A-Am)*(A-Am)))/(n)
	varB <- (sum((B-Bm)*(B-Bm)))/(n)
	varX <- (sum((X-Xm)*(X-Xm)))/(n)
	varY <- (sum((Y-Ym)*(Y-Ym)))/(n)
	covABXY <- sum(((A-Am)+(B-Bm))*((X-Xm)+(Y-Ym)))/n

	# ----------- Calculating the Euclidean regression model
      	euc_par <- data.frame(
		ps = 4L, 
		name = "Euclidean"
		)
	euc_dataMatrix <- matrix(c(vecOne, vecZero, vecZero, vecOne, X, Y, -Y, X), ncol=4)
	euc_target <- matrix(c(A, B),ncol=1)
      	euc_data <- data.frame(
			y = (euc_target),
			x1 = (euc_dataMatrix[,1]),
			x2 = (euc_dataMatrix[,2]),
			x3 = (euc_dataMatrix[,3]),
			x4 = (euc_dataMatrix[,4])
		)
	euc_regression <- lm(y ~ 0 + x1 + x2 + x3 +x4, data=euc_data)
	euc_alpha <- data.frame(
			coeff = c(NA, NA),
			SE = c(NA, NA),
			tValue = c(NA, NA),
			pValue = c(NA, NA)
		)
	euc_beta <- data.frame(
			coeff = c(NA, NA),
			SE = c(NA, NA),
			tValue = c(NA, NA),
			pValue = c(NA, NA)
		)
	
	# retrieving alphas
	for(iAlpha in 1:2) {
		for(iPar in 1:4) {
			euc_alpha[iAlpha, iPar] <- summary(euc_regression)$coeff[iAlpha,iPar]
		}
	}
	# retrieving betas
	for(iBeta in 1:2) {
		for(iPar in 1:4) {
			euc_beta[iBeta , iPar] <- summary(euc_regression)$coeff[iBeta+2,iPar]
		}
	}
	
	# calculating the scale factors (=sigma)
	euc_scaleFactorX <- sqrt(euc_beta$coeff[1]*euc_beta$coeff[1] + euc_beta$coeff[2]*euc_beta$coeff[2])
	euc_scaleFactorY <- euc_scaleFactorX  # no shear (by definition) for the Euclidean solution
                        # ==> same scaling factors for dimension 1 and dimension 2

	# calculating the rotation (angle) (=theta)
	euc_angleRAD <- atan(euc_beta$coeff[2]/euc_beta$coeff[1])
	euc_angleDEG <- euc_angleRAD*180/pi

	# calculating the shear (gamma)
	euc_shear = 0L # per definition shear must be ZERO within an Euclidean geometry

	# calculating the predicted values
	euc_Apred <- euc_alpha$coeff[1]+euc_beta$coeff[1]*X-euc_beta$coeff[2]*Y
	euc_Bpred <- euc_alpha$coeff[2]+euc_beta$coeff[2]*X+euc_beta$coeff[1]*Y
	euc_Xpred <- (euc_alpha$coeff[2]*euc_beta$coeff[2]-B*euc_beta$coeff[2]-
		euc_alpha$coeff[1]*euc_beta$coeff[1]+A*euc_beta$coeff[1])/
		(euc_beta$coeff[1]*euc_beta$coeff[1]-euc_beta$coeff[2]*euc_beta$coeff[2])
	euc_Ypred <- (euc_beta$coeff[1]*B-euc_alpha$coeff[2]*euc_beta$coeff[1]+
		euc_alpha$coeff[1]*euc_beta$coeff[2]-A*euc_beta$coeff[2])/
		(euc_beta$coeff[1]*euc_beta$coeff[1]+euc_beta$coeff[2]*euc_beta$coeff[2])

	# calculating the bidimensional correlation coefficient
	euc_r <- sqrt(sum(((euc_Apred-Am)*(euc_Apred-Am)) + ((euc_Bpred-Bm)*(euc_Bpred-Bm)))/sum(((A-Am)*(A-Am)) + ((B-Bm)*(B-Bm))))
	euc_rsqr <- euc_r*euc_r

	# conducting the inference statistics following Nakaya (1997)
	euc_F <- ((2*n - euc_par$ps)/2)*(euc_rsqr/(1-euc_rsqr))
	euc_df1 <- 2L
	euc_df2 <- 2*n - euc_par$ps 	# set the degrees of freedom: df1/df2
	euc_p <- pf(euc_F, euc_df1, euc_df2, lower.tail = FALSE, log.p = FALSE)

	# ------- Calculating the distortion index following Waterman and Gordon (1984),
	# adjusted by Friedman and Kohler (2003)
	# --- first: calculating distortion index for original configuration
	euc_dDistanceXY <- sqrt(sum((X-euc_Xpred)*(X-euc_Xpred))+sum((Y-euc_Ypred)*(Y-euc_Ypred)))
	euc_dDistanceXYSqr <- euc_dDistanceXY*euc_dDistanceXY
	euc_dMaxXYSqr <- sum((X-Xm)*(X-Xm)+((Y-Ym)*(Y-Ym))) 
	euc_dMaxXY <- sqrt(euc_dMaxXYSqr)
	# --- second: calculating distortion index for target configuration
	euc_dDistanceAB <- sqrt(sum((A-euc_Apred)*(A-euc_Apred))+sum((B-euc_Bpred)*(B-euc_Bpred)))
	euc_dDistanceABSqr <- euc_dDistanceAB*euc_dDistanceAB
	euc_dMaxABSqr <- sum((A-Am)*(A-Am)+((B-Bm)*(B-Bm)))  # referring to target configuration
	euc_dMaxAB <- sqrt(euc_dMaxABSqr)
	euc_diABSqr <- euc_dDistanceABSqr/euc_dMaxABSqr
	euc_diAB <- sqrt(euc_diABSqr)

	euc_diXYSqr <- euc_dDistanceABSqr/euc_dMaxXYSqr
	euc_diXY <- sqrt(euc_diXYSqr)

	# ------- Calculation of DAIC (Difference AIC = Akaike Information Criterion)
	#      DAICso: AIC difference DAICso between a bidimensional regression model and
	#              the bidimensional null model
	#              --> if DAICso < 0, the bidimensional regression model is 
	#                  better than the bidimensional null model.
	#      [calculation according to Nakaya (1997), formula 56]
	euc_dAICso <- 2L*n*log(1-euc_rsqr)+2L*(euc_par$ps-2L)


	# +++++++ end of Euclidean regression model


	# ----------- Calculating the Affine regression model	
      	aff_par <- data.frame(
		ps = 6L, 
		name = "Affine"
		)
	aff_dataMatrix <- matrix(c(
		vecOne, vecZero, 
		vecZero, vecOne, 
		X, vecZero, 
		Y, vecZero, 
		vecZero, X, 
		vecZero, Y), ncol=6)
	aff_target <- matrix(c(A, B),ncol=1)
      	aff_data <- data.frame(
			y = (aff_target),
			x1 = (aff_dataMatrix[,1]),
			x2 = (aff_dataMatrix[,2]),
			x3 = (aff_dataMatrix[,3]),
			x4 = (aff_dataMatrix[,4]),
			x5 = (aff_dataMatrix[,5]),
			x6 = (aff_dataMatrix[,6])
		)
	aff_regression <- lm(y ~ 0 + x1 + x2 + x3 + x4 + x5 + x6, data=aff_data)
	aff_alpha <- data.frame(
			coeff = c(NA, NA),
			SE = c(NA, NA),
			tValue = c(NA, NA),
			pValue = c(NA, NA)
		)
	aff_beta <- data.frame(
			coeff = c(NA, NA),
			SE = c(NA, NA),
			tValue = c(NA, NA),
			pValue = c(NA, NA)
		)
	
	# retrieving alphas
	for(iAlpha in 1:2) {
		for(iPar in 1:4) {
			aff_alpha[iAlpha, iPar] <- summary(aff_regression)$coeff[iAlpha,iPar]
		}
	}
	# retrieving betas
	for(iBeta in 1:4) {
		for(iPar in 1:4) {
			aff_beta[iBeta , iPar] <- summary(aff_regression)$coeff[iBeta+2,iPar]
		}
	}

	# calculating the rotation (angle) (=theta)
	aff_angleRAD <- atan(aff_beta$coeff[3]/aff_beta$coeff[1])
	aff_angleDEG <- aff_angleRAD*180/pi
	if (aff_beta$coeff[1] < 0)
		{
    		aff_angleDEG <- aff_angleDEG+180
		}
	
	# calculating the shear (gamma)
	aff_shear <- (((aff_beta$coeff[4]/aff_beta$coeff[2])*sin(aff_angleRAD))+cos(aff_angleRAD))/(((aff_beta$coeff[4]/aff_beta$coeff[2])*cos		(aff_angleRAD))-sin(aff_angleRAD))

	# calculating the scale factors (=sigma)
	aff_scaleFactorX <- sqrt(aff_beta$coeff[1]*aff_beta$coeff[1]+aff_beta$coeff[3]*aff_beta$coeff[3])
	if (is.nan(aff_shear)) 
    		{
		aff_shear <- (aff_beta$coeff[1]-cos(aff_angleRAD)*aff_scaleFactorX)/aff_beta$coeff[3]
		} 
	if (is.nan(aff_shear)) 
    		{
    		aff_shear <- (sin(aff_angleRAD)*aff_scaleFactorX+aff_beta$coeff[2])/aff_beta$coeff[4]
		}
  	aff_scaleFactorY <- aff_beta$coeff[2]/(aff_shear*cos(aff_angleRAD)-sin(aff_angleRAD))
	if (is.nan(aff_scaleFactorY))
    		{
	    	aff_scaleFactorY <- aff_scaleFactorX
		}

	# calculating the predicted values
	aff_Apred <- aff_alpha$coeff[1]+aff_beta$coeff[1]*X+aff_beta$coeff[2]*Y
	aff_Bpred <- aff_alpha$coeff[2]+aff_beta$coeff[3]*X+aff_beta$coeff[4]*Y
	aff_Xpred <- -(aff_Bpred*aff_beta$coeff[2]-aff_Apred*aff_beta$coeff[4]-
		aff_alpha$coeff[2]*aff_beta$coeff[2]+aff_alpha$coeff[1]*aff_beta$coeff[4])/
		(aff_beta$coeff[2]*aff_beta$coeff[3]+aff_beta$coeff[1]*aff_beta$coeff[4])
	aff_Ypred <- -(aff_Apred*aff_beta$coeff[3]-aff_Bpred*aff_beta$coeff[1]-
		aff_alpha$coeff[1]*aff_beta$coeff[3]+aff_alpha$coeff[2]*aff_beta$coeff[1])/
		(aff_beta$coeff[2]*aff_beta$coeff[3]+aff_beta$coeff[1]*aff_beta$coeff[4])

	# calculating the bidimensional correlation coefficient
	aff_r <- sqrt(sum(((aff_Apred-Am)*(aff_Apred-Am))+((aff_Bpred-Bm)*(aff_Bpred-Bm)))/sum(((A-Am)*(A-Am))+((B-Bm)*(B-Bm))))
	aff_rsqr <- aff_r*aff_r

	# conducting the inference statistics according to 
	# Nakaya (1997), formula 50
	aff_F <- ((2*n-aff_par$ps)/4) * (aff_rsqr/(1-aff_rsqr))
	aff_df1 <- 4L
	aff_df2 <- 2*n-aff_par$ps;   # set the degrees of freedom: df1/df2
	aff_p <-pf(aff_F, aff_df1, aff_df2, lower.tail = FALSE, log.p = FALSE)

	# ------- Calculating the distortion index following Waterman and Gordon (1984),
	# adjusted by Friedman and Kohler (2003)
	# --- first: calculating distortion index for original configuration
	aff_dDistanceXY <- sqrt(sum((X-aff_Xpred)*(X-aff_Xpred))+sum((Y-aff_Ypred)*(Y-aff_Ypred)))
	aff_dDistanceXYSqr <- aff_dDistanceXY*aff_dDistanceXY
	aff_dMaxXYSqr <- sum((X-Xm)*(X-Xm)+((Y-Ym)*(Y-Ym))) 
	aff_dMaxXY <- sqrt(aff_dMaxXYSqr)
	# --- second: calculating distortion index for target configuration
	aff_dDistanceAB <- sqrt(sum((A-aff_Apred)*(A-aff_Apred))+sum((B-aff_Bpred)*(B-aff_Bpred)))
	aff_dDistanceABSqr <- aff_dDistanceAB*aff_dDistanceAB
	aff_dMaxABSqr <- sum((A-Am)*(A-Am)+((B-Bm)*(B-Bm)))  # referring to target configuration
	aff_dMaxAB <- sqrt(aff_dMaxABSqr)
	aff_diABSqr <- aff_dDistanceABSqr/aff_dMaxABSqr
	aff_diAB <- sqrt(aff_diABSqr)

	aff_diXYSqr <- aff_dDistanceABSqr/aff_dMaxXYSqr
	aff_diXY <- sqrt(aff_diXYSqr)

	# ------- Calculation of DAIC (Difference AIC = Akaike Information Criterion)
	#      DAICso: AIC difference DAICso between a bidimensional regression model and
	#              the bidimensional null model
	#              --> if DAICso < 0, the bidimensional regression model is 
	#                  better than the bidimensional null model.
	#      [calculation according to Nakaya (1997), formula 56]
	aff_dAICso <- 2*n*log(1-aff_rsqr)+2*(aff_par$ps-2)

	#+++++++++++ end of affine solution


	# ---- Calculation of DAICs (Difference AIC = Akaike Information Criterion)
	#      between the different fitted bidimensional regression models
	#      [see Nakaya (1997), table 4]
	# -- comparative test between Euclidean and Affine regression model
	dAICea <- 2*n*log((1-aff_rsqr)/(1-euc_rsqr))+2*(aff_par$ps-euc_par$ps)   
	f_ea <- ((2*n-aff_par$ps)/(aff_par$ps-euc_par$ps))*
		((aff_rsqr-euc_rsqr)/(1-aff_rsqr))
	df1_ea <- as.integer(aff_par$ps-euc_par$ps)
	df2_ea <- as.integer(2*n-aff_par$ps)
	p_ea <- pf(f_ea, aff_df1, aff_df2, lower.tail = FALSE, log.p = FALSE)

	# df of parameter t-tests
	tTestDF <- aff_regression$df[1]

	# stop timer for calculating processing time
	timerStop <- proc.time()
      	timeDeltaMS <- timerStop[3] - timerStart[3]


      	# return all the results to a data.frame
	res_euc <- data.frame(r=euc_r, rsqr=euc_rsqr, diABSqr=euc_diABSqr, dMaxABSqr=euc_dMaxABSqr, diXYSqr=euc_diXYSqr, dMaxXYSqr=euc_dMaxXYSqr, 	scaleFactorX=euc_scaleFactorX, scaleFactorY=euc_scaleFactorY, angleDEG=euc_angleDEG, shear=euc_shear, ttestDF=tTestDF, alpha1=euc_alpha[1,], 	alpha2=euc_alpha[2,], beta1=euc_beta[1,], beta2=euc_beta[2,], beta3=euc_beta[3,], beta4=euc_beta[4,], fValue=euc_F, df1=euc_df1, df2=euc_df2, 	pValue=euc_p, dAICso=euc_dAICso)
	res_aff <- data.frame(r=aff_r, rsqr=aff_rsqr, diABSqr=aff_diABSqr, dMaxABSqr=aff_dMaxABSqr, diXYSqr=aff_diXYSqr, dMaxXYSqr=aff_dMaxXYSqr, 	scaleFactorX=aff_scaleFactorX, scaleFactorY=aff_scaleFactorY, angleDEG=aff_angleDEG, shear=aff_shear, ttestDF=tTestDF, alpha1=aff_alpha[1,], 	alpha2=aff_alpha[2,], beta1=aff_beta[1,], beta2=aff_beta[2,], beta3=aff_beta[3,], beta4=aff_beta[4,], fValue=aff_F, df1=aff_df1, df2=aff_df2, 	pValue=aff_p, dAICso=aff_dAICso)
	euclideanVSaffine <- data.frame(dAIC=dAICea, fValue=f_ea, pValue=p_ea, df1=df1_ea, df2=df2_ea)
	
	results_sum <- data.frame(euclidean=res_euc, affine=res_aff, eucVSaff=euclideanVSaffine)
	class(results_sum) <- "BiDimRegression"

	invisible(results_sum)   # returns the measures of the bidimensional regression
	
	############################################################################################################
	# Copyright ©2012-2013. Claus-Christian Carbon (CCC). All Rights Reserved. 
	# Permission to use, copy, modify, and distribute this software and its documentation for 
	# educational, research, and not-for-profit purposes, without fee and without a signed licensing agreement,
	# is hereby granted, provided that the above copyright notice, this paragraph and the following 
	# two paragraphs appear in all copies, modifications, and distributions. 
	# Contact the author for commercial licensing opportunities.
	#
	# Created by Claus-Christian Carbon, Department of General Psychology and Methodology, 
	# University of Bamberg, Bavaria, Germany
	#
	# IN NO EVENT SHALL THE AUTHOR BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR 
	# CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS 
	# DOCUMENTATION, EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	#
	# THE AUTHOR SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
	# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION,
	# IF ANY, PROVIDED HEREUNDER IS PROVIDED 'AS IS'. THE AUTHOR HAS NO OBLIGATION TO PROVIDE MAINTENANCE,
	# SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
	############################################################################################################
}
