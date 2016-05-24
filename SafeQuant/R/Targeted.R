
#' S3 class object describing a calibration curve and storing some figures of merit
#' @param eset ExpressionSet
#' @param method to calculate Limit of Detection / Limit of Quantification. c("blank","low") 
#' @return calibrationCurve object
#' @export
#' @import Biobase
#' @importFrom stats qt
#' @note  No note
#' @details No details
#' @references Statistical characterization of multiple-reaction monitoring mass spectrometry (MRM-MS) assays for quantitative proteomics, Mani et al. (2012), \url{http://www.ncbi.nlm.nih.gov/pubmed/23176545} 
#' @examples print("No examples")
calibrationCurve <- function(eset,method="blank"){
	
	out <- list()
	class(out) <- "calibrationCurve"
	
	### store the curve
	out$curve <- data.frame(concentration=rep(fData(eset)$concentration[fData(eset)$concentration > 0],ncol(eset))
			, intensity= as.vector(unlist(exprs(eset)[fData(eset)$concentration > 0,]))	)
	logCurve <- log10(out$curve)
	out$fit <- lm(intensity ~ concentration, data=logCurve)
	
	out$lod <- NA
	out$loq <- NA
	### store ExpressionSet matching assay
	out$eset <- eset
	#print(names( fData(out$eset)))
	#[1] "Peptide.Sequence" "Protein.Name"     "Precursor.Mz"     "Precursor.Charge"
	#[5] "Retention.Time"   "concentration"    "dilutionCurveId" 

	out$label <- paste(fData(out$eset)$dilutionCurveId[1],fData(out$eset)$Protein.Name[1])
	
	if("blank" %in% method ){
		
		# LIMIT OF DETECTION & LIMIT OF QUANTIFICATION
		
		# Blank Sample Method, Mani et al. http://www.biomedcentral.com/1471-2105/13/S16/S9
		# Assuming that random measurement errors are normally distributed, and with 5% 
		# risk of incorrectly claiming detection in the absence of analyte (alpha) or missing the detection of analyte (beta),
		# LOD = 3.29 sdB and LOQ = 3 x LOD = 10 sigmaB where sigmaB is the standard deviation of the blank sample.
		#
		# OR
		#Currie LA: Limits for qualitative detection and quantitative determination. Application to radiochemistry. Analytical Chemistry 1968, 40(3):586-593.
		
		# In Mani et al. the LOD is calculated as the sd of "blank (light r.t. region) to heavy peptide ratio". This given the lod as a fraction of the heavy peptide concentration.
		# We calculate the LOD as a fraction of the most intense peptide concentration
		
		maxConc <-  max(fData(eset)$concentration,na.rm=T)
		blankInt <- exprs(eset)[fData(eset)$concentration == 0,] 
		maxIntMedian <-  median( exprs(eset)[fData(eset)$concentration %in% maxConc,] ,na.rm=T) 
		blankMeasuredConc <- (blankInt/maxIntMedian)*maxConc
		
		out$lod <- 3.29 * sd(blankMeasuredConc,na.rm=T)
		out$loq <- 3 * out$lod
	}else if("low" %in% method ){
		
		# Blank and low Concentration Sample Method, Mani et al. http://www.biomedcentral.com/1471-2105/13/S16/S9 and http://www.nature.com/nbt/journal/v27/n7/full/nbt.1546.html
		# LOD = meanB +t(1-b) (sdB + sdS)/sqrt(n)
		
		maxConc <-  max(fData(eset)$concentration,na.rm=T)
		blankInt <- exprs(eset)[fData(eset)$concentration == 0,] 
		lowConc <- min(fData(eset)$concentration[fData(eset)$concentration > 0],na.rm=T)
		lowInt <- exprs(eset)[fData(eset)$concentration %in% lowConc,]
		maxIntMedian <-  median( exprs(eset)[fData(eset)$concentration %in% maxConc,] ,na.rm=T) 
		blankMeasuredConc <- (blankInt/maxIntMedian)*maxConc
		lowMeasuredConc <- (lowInt/maxIntMedian)*maxConc
		n <- length(blankMeasuredConc)+length(lowMeasuredConc)
		out$lod <- mean(blankMeasuredConc,na.rm=T) +  (qt(0.95,n-1) * (sd(blankMeasuredConc,na.rm=T)+sd(lowMeasuredConc,na.rm=T))/sqrt(n))
		out$loq <- 3 * out$lod
			
	}
	
	return(out)
	
}

#' @export
#' @importFrom graphics plot.default
plot.calibrationCurve <- function(x, pch=19, lwd=2,cex.axis=1.5,cex.lab=1.5,cex=1.5,cex.main=1.5,xlab="Concentration",ylab="Area Under Curve",...){
	
	### include lod if out of range
	xlim=c(min(x$lod,min(x$curve[,1],na.rm=T),na.rm=T),max(x$curve[,1],na.rm=T))
	
	plot.default(x$curve
			, log="xy"
			, main=x$label
 			, pch=pch
			, cex.axis=cex.axis 
			, cex.lab=cex.lab
			, cex=cex
			, cex.main=1.5
			, xlab=xlab
			, ylab=ylab
			, xlim=xlim
			, ...
		
	)
	
	### linear model for points above loq
	slope <- NA
	selData <- log10(x$curve[x$curve[,1] > x$loq  ,])
	if(length(unique(selData[,1])) > 1){
		fit <- lm(intensity ~ concentration, data=selData)
		slope <- coef(fit)[2]
		abline(fit, col="darkgrey",lwd=lwd)
	}
	
	abline(v=x$lod, col="red",lwd=lwd, lty=2)
	abline(v=x$loq, col="blue",lwd=lwd, lty=2)
	
	### bottom right 
	digits <- 3
	legend("bottomright"
			,legend=c(as.expression(bquote(R^2*"" == .(round(summary(x$fit)$r.squared,digits)))
					)
					,paste("K=",round(slope,digits)) 
					,paste("LOD=",round(x$lod,digits)) 
					,paste("LOQ=",round(x$loq,digits)) 
			)
			,text.col=c("black","darkgrey","red","blue"), box.lwd=NA, cex=1.5 )
	
}
