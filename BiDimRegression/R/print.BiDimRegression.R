print.BiDimRegression <-
function(x, ...)
{
	# version of the routine
	verRoutine <- "1.0.6"

      # ---- OUTPUT the data ---- 
	for(iOutputLoop in 0:1) {
		if (iOutputLoop == 1) {
			# declaration of the output file 	
			sink(file="BiDimRegrOutput.txt", type="output", append=TRUE)
		}

      # print the results of the BiDimensional Regression Analysis
	cat(sprintf('\n-------------------\n'))
	cat(sprintf('BiDimRegression %s \n', verRoutine))
	cat(sprintf('Date-Time: %s\n', date()))
	cat(sprintf('----------------------------------------------------------------------\n'))

	euc_alpha=matrix(c(x$euclidean.alpha1.coeff, x$euclidean.alpha1.SE, x$euclidean.alpha1.tValue, 
		x$euclidean.ttestDF, x$euclidean.alpha1.pValue, 
            x$euclidean.alpha2.coeff, x$euclidean.alpha2.SE, x$euclidean.alpha2.tValue, 
		x$euclidean.ttestDF, x$euclidean.alpha2.pValue),
	      nrow=2, byrow = TRUE)
	euc_beta=matrix(c(x$euclidean.beta1.coeff, x$euclidean.beta1.SE, 
		x$euclidean.beta1.tValue, x$euclidean.ttestDF, x$euclidean.beta1.pValue,
 		x$euclidean.beta2.coeff, x$euclidean.beta2.SE, x$euclidean.beta2.tValue, 
		x$euclidean.ttestDF, x$euclidean.beta2.pValue),
		nrow=2, byrow = TRUE)
	aff_alpha=matrix(c(x$affine.alpha1.coeff, x$affine.alpha1.SE, x$affine.alpha1.tValue, 
		x$affine.ttestDF, x$affine.alpha1.pValue,
	      x$affine.alpha2.coeff, x$affine.alpha2.SE, x$affine.alpha2.tValue, 
		x$affine.ttestDF, x$affine.alpha2.pValue),
		nrow=2, byrow = TRUE)
	aff_beta=matrix(c(x$affine.beta1.coeff, x$affine.beta1.SE, x$affine.beta1.tValue, 
		x$affine.ttestDF, x$affine.beta1.pValue,
 		x$affine.beta2.coeff, x$affine.beta2.SE, x$affine.beta2.tValue, x$affine.ttestDF, x$affine.beta2.pValue,
		x$affine.beta3.coeff, x$affine.beta3.SE, x$affine.beta3.tValue, x$affine.ttestDF, x$affine.beta3.pValue,
		x$affine.beta4.coeff, x$affine.beta4.SE, x$affine.beta4.tValue, x$affine.ttestDF, x$affine.beta4.pValue),
            nrow=4, byrow = TRUE)
	
	# -- Overall analysis
	overallStats=matrix(c(
            x$euclidean.r, x$euclidean.rsqr, x$euclidean.fValue, x$euclidean.df1, x$euclidean.df2, x$euclidean.pValue, 
		x$affine.r, x$affine.rsqr, x$affine.fValue, x$affine.df1, x$affine.df2, x$affine.pValue),
		nrow=2, byrow = TRUE)
      colnames(overallStats) <- c("r", "r-sqr", "F-value", "df1", "df2", "p-value")
	rownames(overallStats) <- c("Euclidean", "Affine")

	cat(sprintf('\n--- overall statistics ---\n'))
	printCoefmat(overallStats, P.values=TRUE, has.Pvalue=TRUE, digits=3, tst.ind=4, signif.stars=TRUE)

	# show warning for unidentified model
	if (x$euclidean.df2 < x$euclidean.df1)	{
		cat(sprintf('WARNING: Euclidean model is not defined'))
		}
	if (x$affine.df2 < x$affine.df1) {	
		cat(sprintf('\t\t\t\t\tWARNING: Affine model is not defined\n'))
		}

	cat(sprintf('----------------------------------------------------------------------\n'))
	cat(sprintf('\n--- parameters ---\n'))
	cat(sprintf('\n- Euclidean -\n'))

	# -- Parameter analyses
	colNames <- c("Parameter", "Std.Err", "t-value", "df", "p-value")
	rowNamesAlphas <- c("alpha1", "alpha2")
	rowNames2Betas <- c("beta1 ", "beta2 ")
	rowNames4Betas <- c("beta1 ", "beta2 ", "beta3 ", "beta4 ")
	colnames(euc_alpha) <- colNames 
	colnames(euc_beta) <- colNames 
	colnames(aff_alpha) <- colNames 
	colnames(aff_beta) <- colNames 
	rownames(euc_alpha) <- rowNamesAlphas
	rownames(euc_beta) <- rowNames2Betas
	rownames(aff_alpha) <- rowNamesAlphas
	rownames(aff_beta) <- rowNames4Betas

	printCoefmat(euc_alpha, P.values = TRUE, has.Pvalue = TRUE, digits=3, tst.ind=4)
	cat(sprintf('\n'))
	printCoefmat(euc_beta, P.values = TRUE, has.Pvalue = TRUE, digits=3, tst.ind=4)
	cat(sprintf('\n'))

	cat(sprintf('\n- Affine -\n'))
	printCoefmat(aff_alpha, P.values = TRUE, has.Pvalue = TRUE, digits=3, tst.ind=4)
	cat(sprintf('\n'))
	printCoefmat(aff_beta, P.values = TRUE, has.Pvalue = TRUE, digits=3, tst.ind=4)
	cat(sprintf('\n'))

	cat(sprintf('----------------------------------------------------------------------\n'))
	cat(sprintf('\n--- details ---\n'))
	cat(sprintf('\n- Euclidean -\t\t\t- Affine -\n'))
	cat(sprintf('scaleX\t= scaleY = %4.3f\tscaleX\t= %4.3f, scaleY = %4.3f\n', 
		x$euclidean.scaleFactorX, x$affine.scaleFactorX, x$affine.scaleFactorY))
	cat(sprintf('shear\t= %4.3f\t\t\tshear\t= %4.3f\n', 
		x$euclidean.shear, x$affine.shear))
	cat(sprintf('angle\t= %4.3f DEG\t\tangle\t= %4.3f DEG\n', 
		x$euclidean.angleDEG, x$affine.angleDEG))

	cat(sprintf('---\t\t\t\t---\n'))
	cat(sprintf('DAIC (agst.0)\t= %4.2f\tDAIC (agst.0)\t= %4.2f\n', 
		x$euclidean.dAICso, x$affine.dAICso))

	cat(sprintf('---\t\t\t\t---\n'))
	cat(sprintf('dMaxABSqr\t= %4.3f\tdMaxABSqr\t= %4.3f\n', 
		x$euclidean.dMaxABSqr, x$affine.dMaxABSqr))
	cat(sprintf('diABSqr  \t= %4.3f\t\tdiABSqr  \t= %4.3f\n', 
		x$euclidean.diABSqr, x$affine.diABSqr))
	cat(sprintf('dMaxXYSqr\t= %4.3f\tdMaxXYSqr\t= %4.3f\n', 
		x$euclidean.dMaxXYSqr, x$affine.dMaxXYSqr))
	cat(sprintf('diXYSqr  \t= %4.3f\t\tdiXYSqr  \t= %4.3f\n', 
		x$euclidean.diXYSqr, x$affine.diXYSqr))
	cat(sprintf('----------------------------------------------------------------------\n'))
	cat(sprintf('\n--- comparative statistics of fitted bidimensional regxsion models\n\n'))

 	comparativeStats=matrix(c(
            x$eucVSaff.fValue, as.integer(x$eucVSaff.df1), as.integer(x$eucVSaff.df2), x$eucVSaff.pValue),
		nrow=1, byrow = TRUE)
      colnames(comparativeStats) <- c("F-value", "df1", "df2", "p-value")
	rownames(comparativeStats) <- c("Euclidean vs. Affine")
	printCoefmat(comparativeStats, P.values=TRUE, has.Pvalue=TRUE, digits=3, cs.ind=1, signif.stars=TRUE)

	cat(sprintf('\n'))

	if (x$eucVSaff.df2 < x$eucVSaff.df1) {	
		cat(sprintf('WARNING: model is not defined\n'))
		}
	if (x$eucVSaff.pValue <= .05) {
		if (x$eucVSaff.dAIC<0) {
        		superiorSolution <- '(significantly better: Affine solution)'
    		} else {
        		superiorSolution <- '(significantly better: Euclidean solution)'
		}
    	} else { 
    		superiorSolution = '(not significantly different solutions)'
		}

   	cat(sprintf('\nDAICea = %4.3f %s\n\n', x$eucVSaff.dAIC, superiorSolution))
	cat(sprintf('**********************************************************************\n\n\n'))
	}

      sink()   # close the output result file
}
