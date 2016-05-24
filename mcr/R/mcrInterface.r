###############################################################################
##
## mcrInterface.r
##
## Interface function for computing method comparisons.
##
## Copyright (C) 2011 Roche Diagnostics GmbH
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

#' Comparison of Two Measurement Methods Using Regression Analysis
#' 
#' \code{mcreg} is used to compare two measurement methods by means of regression analysis.
#' Available methods comprise ordinary and weighted linear regression,
#' Deming and weighted Deming regression and Passing-Bablok regression. Point estimates of regression
#' parameters are computed together with their standard errors and confidence intervals.
#'
#' The regression analysis yields regression coefficients 'Inercept' and 'Slope' of the regression
#' \eqn{Testmethod = Intercept + Slope * Referencemethod}. There are methods for computing the systematical
#' bias between reference and test method at a decision point Xc, \eqn{Bias(Xc) = Intercept + (Slope-1) * Xc},
#' accompanied by its corresponding standard error and confidence interval. One can use plotting
#' method \code{plotBias} for a comprehensive view of the systematical bias.
#'   
#' Weighted regression for heteroscedastic data is available for linear and Deming regression and
#' implemented as a data point weighting with the inverted squared value of the reference method.
#' Therefore calculation of weighted regression (linear and Deming) is available only for positive values (>0).
#' Passing-Bablok regression is only available for non-negative values (>=0).
#'
#' Confidence intervals for regression parameters and bias estimates are calculated
#' either by using analytical methods or by means of resampling methods ("jackknife", "bootstrap", "nested bootstrap").
#' An analytical method is available for all types of regression except for weighted Deming.
#' For Passing-Bablok regression the option "analytical" calculates confidence intervals for the regression parameters
#' according to the non-parametric approach given in the original reference.
#' 
#' The "jackknife" (or leave one out resampling) method was suggested by Linnet for
#' calculating confidence intervals of regression parameters
#' of Deming and weighted Deming regression. It is possible to calculate jackknife confidence
#' intervals for all types of regression. Note that we do not recommend this method for Passing-Bablok
#' since it has a tendency of underestimating the variability
#' (jackknife is known to yield incorrect estimates for errors of quantiles).
#' 
#' The bootstrap method requires additionally choosing a value for \code{method.bootstrap.ci}.
#' If bootstrap is the method of choice, "BCa", t-bootstrap ("tBoot") and simple "quantile" confidence intervals are recommended
#' (See Efron B. and Tibshirani R.J.(1993),Carpenter J., Bithell J. (2000)).
#' The "nestedbootstrap" method can be very time-consuming but is necessary for calculating t-bootstrap
#' confidence intervals for weighted Deming or Passing-Bablok regression. For these regression methods there are no analytical
#' solutions for computing standard errors, which therefore have to be obtained by nested bootstrapping.
#' 
#' Note that estimating resampling based confidence intervals for Passing-Bablok regressions can take very long
#' for larger data sets due to the high computational complexity of the algorithm. To mitigate this drawback
#' an adaption of the Passing-Bablok algorithm has been implemented (\code{"PaBaLarge"}), which yields approximative results. This approach
#' does not build the complete upper triangular matrix of all 'n*(n-1)/2' slopes. It subdivides the range of slopes into
#' 'NBins' classes, and sorts each slope into one of these bins. The remaining steps are the same as for the exact \code{"PaBa"}
#' algorithm, except that these are performed on the binned slopes instead of operating on the matrix of slopes.
#' 
#' Our implementation of exact Passing-Bablok regression (\code{"PaBa"}) provides two alternative metrics for regression slopes which
#' can result in different regression estimates.
#' As a robust regression method PaBa is essentially invariant to the parameterization of regression slopes,
#' however in the case of an even number of all pairwise slopes the two central slopes are averaged to estimate the final regression slope.
#' In this situation using an angle based metric (\code{slope.measure="radian"}) will result in
#' a regression estimate that is geometrically centered between the two central slopes, whereas the tangent measure (\code{slope.measure="tangent"})
#' proposed in Passing and Bablok (1983) will be geometrically biased towards a higher slope. See below for a pathological example.
#' Note that the difference between the two measures is neglectable for data sets with reasonable sample size (N>20) and correlation.   
#'
#'
#' @param x measurement values of reference method, or two column matrix.
#' @param y measurement values of test method.
#' @param mref.name name of reference method (Default "Method1").
#' @param mtest.name name of test Method (Default "Method2").
#' @param sample.names names of cases (Default "S##").
#' @param error.ratio ratio between squared measurement errors of reference and 
#'        test method, necessary for Deming regression (Default 1).
#' @param alpha value specifying the 100(1-alpha)\% confidence level for confidence intervals (Default is 0.05).
#' @param method.reg regression method.  It is possible to choose between five regression methods: 
#'          \code{"LinReg"} - ordinary least square regression.\cr 
#'          \code{"WLinReg"} - weighted ordinary least square regression.\cr
#'          \code{"Deming"} - Deming regression.\cr 
#'          \code{"WDeming"} - weighted Deming regression.\cr
#'          \code{"PaBa"} - Passing-Bablok regression.\cr  
#'          \code{"PaBaLarge"} - approximative Passing-Bablok regression for large datasets, operating on \code{NBins} classes of constant slope angle which each slope
#'                               is classified to instead of building the complete triangular matrix of all N*N/2 slopes.
#' @param method.ci method of confidence interval calculation. The function 
#'         contains four basic methods for calculation of confidence intervals for regression coefficients.
#'          \code{"analytical"} - with parametric method.\cr
#'          \code{"jackknife"} - with leave one out resampling.\cr
#'          \code{"bootstrap"} - with ordinary non-parametric bootstrap resampling.\cr
#'          \code{"nested bootstrap"} - with ordinary non-parametric bootstrap resampling.\cr 
#' @param method.bootstrap.ci bootstrap based confidence interval estimation method. 
#' @param nsamples number of bootstrap samples.
#' @param nnested number of nested bootstrap samples.
#' @param rng.seed integer number that sets the random number generator seed for bootstrap sampling. If set to NULL currently in the R session used RNG setting will be used. 
#' @param rng.kind type of random number generator for bootstrap sampling. Only used when rng.seed is specified, see set.seed for details.
#' @param iter.max maximum number of iterations for weighted Deming iterative algorithm.
#' @param threshold numerical tolerance for weighted Deming iterative algorithm convergence.
#' @param na.rm remove measurement pairs that contain missing values (Default is FALSE).
#' @param NBins number of bins used when 'reg.method="PaBaLarge"' to classify each slope in one of 'NBins' bins covering the range of all slopes
#' @param slope.measure angular measure of pairwise slopes used for exact PaBa regression (see below for details).\cr   
#'          \code{"radian"} - for data sets with even sample numbers median slope is calculated as average of two central slope angles.\cr
#'          \code{"tangent"} - for data sets with even sample numbers median slope is calculated as average of two central slopes (tan(angle)).\cr
#' @return "MCResult" object containing regression results. The function \code{\link{getCoefficients}} or
#'          \code{\link{printSummary}} can be used to obtain or print a summary of the results.
#'          The function \code{\link{getData}} allows to see the original data.
#'          An S4 object of class "MCResult" containing at least the following slots: 
#'  \item{data}{measurement data in wide format, one pair of observations per sample. Includes samples ID, 
#'              reference measurement, test measurement.}
#'  \item{para}{numeric matrix with estimates for slope and intercept, corresponding standard deviations and confidence intervals.}
#'  \item{mnames}{character vector of length two containing names of analytical methods.}
#'  \item{regmeth}{type of regression type used for parameter estimation.}
#'  \item{cimeth}{method used for calculation of confidence intervals.}
#'  \item{error.ratio}{ratio between squared measurement errors of reference  and test method, necessary for Deming regression.}
#'  \item{alpha}{confidence level using for calculation of confidence intervals.}
#' @references  Bland, J. M., Altman, D. G. (1986) 
#'              Statistical methods for assessing agreement between two methods of clinical measurement. 
#'              \emph{Lancet}, \bold{i:} 307--310.
#'  
#'              Linnet, K. (1993)
#'              Evaluation of Regression Procedures for Methods Comparison Studies.
#'              \emph{CLIN. CHEM.} \bold{39/3}, 424--432.
#'  
#'              Linnet, K. (1990)
#'              Estimation of the Linear Relationship between the Measurements of two Methods with Proportional Errors.
#'              \emph{STATISTICS IN MEDICINE}, Vol. \bold{9}, 1463--1473.
#' 
#'              Neter, J., Wassermann, W., Kunter, M. (1985)
#'              \emph{Applied Statistical Models.} 
#'              Richard D. Irwing, INC.
#' 
#'              Looney, S. W. (2010) 
#'              Statistical Methods for Assessing Biomarkers.
#'              \emph{Methods in Molecular Biology}, vol. \bold{184}: \emph{Biostatistical Methods}. Human Press INC.
#' 
#'              Passing, H., Bablok, W. (1983) 
#'              A new biometrical procedure for testing the equality of measurements from two different analytical methods. 
#'              Application of linear regression procedures for method comparison studies in clinical chemistry, Part I. 
#'              \emph{J Clin Chem Clin Biochem}.  Nov; \bold{21(11)}:709--20.
#' 
#'              Efron, B., Tibshirani, R.J. (1993)
#'              \emph{An Introduction to the Bootstrap}. Chapman and Hall.
#' 
#'              Carpenter, J., Bithell, J. (2000)
#'              Bootstrap confidence intervals: when, which, what? A practical guide for medical statisticians.
#'              \emph{Stat Med}, \bold{19 (9)}, 1141--1164.
#' 
#'              \emph{CLSI EP9-A2}. Method Comparison and Bias Estimation Using Patient Samples; Approved Guideline.
#' @seealso \code{\link{plotDifference}}, \code{\link{plot.mcr}}, \code{\link{getResiduals}}, \code{\link{plotResiduals}}, \code{\link{calcResponse}}, \code{\link{calcBias}}, \code{\link{plotBias}}, \code{\link{compareFit}}
#' @author Ekaterina Manuilova \email{ekaterina.manuilova@@roche.com}, Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}, Fabian Model \email{fabian.model@@roche.com}
#' @examples
#' library("mcr")
#' data(creatinine,package="mcr")
#' x <- creatinine$serum.crea
#' y <- creatinine$plasma.crea
#' # Deming regression fit.
#' # The confidence intercals for regression coefficients
#' # are calculated with analytical method
#' model1<- mcreg(x,y,error.ratio=1,method.reg="Deming", method.ci="analytical",
#'                mref.name = "serum.crea", mtest.name = "plasma.crea", na.rm=TRUE)
#' # Results
#' printSummary(model1)
#' getCoefficients(model1)
#' plot(model1)
#' # Deming regression fit.
#' # The confidence intervals for regression coefficients
#' # are calculated with bootstrap (BCa) method
#' model2<- mcreg(x,y,error.ratio=1,method.reg="Deming",
#'                method.ci="bootstrap", method.bootstrap.ci = "BCa",
#'                mref.name = "serum.crea", mtest.name = "plasma.crea", na.rm=TRUE)
#' compareFit(model1, model2) 
#' 
#' ## Pathological example of Passing-Bablok regression where measure for slope angle matters  
#' x1 <- 1:10; y1 <- 0.5*x1; x <- c(x1,y1); y <- c(y1,x1) 
#' m1 <- mcreg(x,y,method.reg="PaBa",method.ci="analytical",slope.measure="radian",
#'             mref.name="X",mtest.name="Y")
#' m2 <- mcreg(x,y,method.reg="PaBa",method.ci="analytical",slope.measure="tangent",
#'             mref.name="X",mtest.name="Y")
#' plot(m1, add.legend=FALSE,identity=FALSE,
#'      main="Radian vs. tangent slope measures in Passing-Bablok regression\n(pathological example)",
#'      ci.area=FALSE,add.cor=FALSE)
#' plot(m2, ci.area=FALSE,reg.col="darkgreen",reg.lty=2,identity=FALSE,add.legend=FALSE,
#'      draw.points=FALSE,add=TRUE,add.cor=FALSE)
#' includeLegend(place="topleft",models=list(m1,m2),model.names=c("PaBa Radian","PaBa Tangent"),
#'               colors=c("darkblue","darkgreen"),lty=c(1,2),design="1",digits=2)

mcreg <- function(x,y=NULL,error.ratio=1,alpha=0.05,
              mref.name = NULL,
              mtest.name = NULL, 
              sample.names = NULL,
              method.reg=c("PaBa","LinReg","WLinReg","Deming","WDeming", "PaBaLarge"),
              method.ci=c("bootstrap","jackknife","analytical","nestedbootstrap"),
              method.bootstrap.ci=c("quantile","Student","BCa","tBoot"),
              nsamples=999,nnested=25,rng.seed=NULL,rng.kind="Mersenne-Twister",
              iter.max=30,threshold=0.000001,na.rm=FALSE, NBins=1000000,slope.measure=c("radian","tangent"))
{
	## Match choice parameters
    method.reg <- match.arg(method.reg)
    method.ci <- match.arg(method.ci)
    method.bootstrap.ci <- match.arg(method.bootstrap.ci)
    slope.measure <- match.arg(slope.measure)
  
	## Check input data
    if( method.reg %in% c("Deming", "WDeming") )
    {
        stopifnot(!is.na(error.ratio))
        stopifnot(is.numeric(error.ratio))
        stopifnot(length(error.ratio) > 0)
        stopifnot(error.ratio >= 0) 
        
        if( method.reg == "WDeming")
        {
            stopifnot(!is.na(iter.max))
            stopifnot(is.numeric(iter.max))
            stopifnot(length(iter.max) > 0)
            stopifnot(round(iter.max) == iter.max)  
            stopifnot(iter.max > 0)
            stopifnot(!is.na(threshold))
            stopifnot(is.numeric(threshold))
            stopifnot(length(threshold) > 0)
            stopifnot(threshold >= 0)
        }
    }
    
    if( method.reg %in% c("bootstrap","nestedbootstrap") )
    {
        stopifnot(!is.na(nsamples))
        stopifnot(is.numeric(nsamples))
        stopifnot(length(nsamples) > 0)
        stopifnot(nsamples > 0)
        
        if( method.reg == "nestedbootstrap")  
        {
            stopifnot(!is.na(nnested))
            stopifnot(is.numeric(nnested))
            stopifnot(length(nnested) > 0)
            stopifnot(nnested > 0)    
        }
    }

    stopifnot(!is.na(alpha))
    stopifnot(is.numeric(alpha))
    stopifnot(length(alpha) > 0)
    stopifnot(alpha>0 & alpha<1)

    if(is.null(y)) 
	{
        stopifnot(is.matrix(x))
        stopifnot(ncol(x)==2)
        stopifnot(nrow(x)>2)
		data <- x
	} 
	else 
	{
        stopifnot(is.numeric(x) & is.numeric(y))
        stopifnot(length(x)==length(y))
        stopifnot(length(x) > 2)	
        stopifnot(!all(is.na(x)) & !all(is.na(y)))
		data <- cbind(x,y)
	}
	
    colnames(data) <- c("x","y")
    npoints.old <- dim(data)[1]
    if(!is.null(sample.names)) 									
        stopifnot(length(sample.names) == nrow(data))
    
	## Remove NAs
    if(na.rm) 
	{
        cc <- complete.cases(data) 
        data <- data[cc,]
        sample.names <- sample.names[cc]
        npoints.new <- dim(data)[1]
        if(npoints.old != npoints.new)
        {
            cat(paste(  "Please note: \n",round(npoints.old-npoints.new), " of ",npoints.old,
                        " observations contain missing values and have been removed.\nNumber of data points in analysis is ",
                        npoints.new,".\n", sep=""))
        }
    }
    stopifnot(!any(is.na(data)))
    stopifnot((sd(data[,"x"])>0)&(sd(data[,"y"])>0))
    
    #----- Limits for data points
    
    if(is.null(mref.name)) 
        mref.name <- "Method1"
    if(is.null(mtest.name)) 
        mtest.name <- "Method2"
    if(mref.name == mtest.name)
        cat("Names of test method and reference method are equal!")
  
    mnames <- c(mref.name, mtest.name)   
    
    ## For PaBa set error ratio 1 for correct weighted residual calculation
    if(method.reg %in% c("PaBa", "PaBaLarge")) error.ratio <- 1

    if(method.reg=="WDeming" & min(data)<=0) 
        stop("Weighted Deming regression for non-positive values is not available.")
    if(method.reg %in% c("PaBa", "PaBaLarge") & min(data)<0) 
        stop("Passing-Bablok regression for negative values is not available.")	
    if(method.reg=="WLinReg" & any(data==0)) 
        stop("Weighted linear regression for zero values is not available.")
    if(method.reg %in% c("WDeming","PaBa") & method.ci=="bootstrap" & method.bootstrap.ci=="tBoot") 
        stop(paste("The analytical estimation of coefficients' SE for", method.reg ,"is not available. \nFor tBoot please choose nested bootstrap."))	

    ## Analytical confidence intervals
    if(method.ci=="analytical") 
    {
        if(method.reg %in% c("PaBa", "PaBaLarge")) 
        {
            ## Determine if slope 1 or -1 is expected
            posCor <- cor(data[,"x"], data[,"y"], method="kendall") >= 0
            if(method.reg == "PaBa")
            {                
                ## Compute slope matrix
                angM <- mc.calcAngleMat(data[,"x"], data[,"y"], posCor=posCor)
                ## Regression
                mc.res <- mc.paba(angM, data[,"x"], data[,"y"], posCor=posCor, alpha=alpha, slope.measure=slope.measure)		
            }
            else    # PaBaLarge: regression using the approximative PaBa-algorithm
            {
                mc.res <- mc.paba.LargeData( data[,"x"], data[,"y"], posCor=posCor, NBins=NBins, alpha=alpha )
            }
            xmean <- as.numeric(NA)
            weight <- rep(1, nrow(data))
        }
        else if(method.reg %in% c("LinReg","WLinReg","Deming")) 
        {
            mcr <- switch(  method.reg,
                            LinReg=mc.linreg(data[,"x"], data[,"y"]),
                            WLinReg=mc.wlinreg( data[,"x"], data[,"y"]),
                            Deming=mc.deming(data[,"x"], data[,"y"], error.ratio))
                    
            mc.res <- mc.analytical.ci(mcr$b0, mcr$b1, mcr$se.b0, mcr$se.b1, nrow(data), alpha)
            xmean <- mcr$xw
            weight <- mcr$weight
        } 
        else 
            stop(paste("The analytical estimation of confidence intervals for", method.reg, "is not available."))

        result <- newMCResultAnalytical(method.names = mnames, sample.names = sample.names, 
                                        wdata=data,para=mc.res, xmean=xmean, regmeth=method.reg,
                                        cimeth=method.ci, error.ratio=error.ratio, alpha=alpha, weight=weight)
    }
    ## Jackknife confidence intervals
    else if(method.ci == "jackknife") 
    {
        cat("Jackknife based calculation of standard error and confidence intervals according to Linnet's method.\n")
        res <- mc.bootstrap(method.reg=method.reg, jackknife=TRUE,bootstrap="none", X=data[,"x"], Y=data[,"y"],
                            error.ratio=error.ratio, nsamples=nsamples, NBins=NBins, slope.measure=slope.measure,
                            threshold=threshold, iter.max=iter.max, nnested=nnested)
        B0jack <- res$B0jack
        B1jack <- res$B1jack
        glob.coef <- res$glob.coef
        jackB0 <- mc.calcLinnetCI(B0jack, glob.coef[1], alpha)
        jackB1 <- mc.calcLinnetCI(B1jack, glob.coef[2], alpha)
        weight <- res$weight
        mc.res <- mc.make.CIframe(b0=jackB0$est, b1=jackB1$est, se.b0=jackB0$se, se.b1=jackB1$se, CI.b0=jackB0$CI, CI.b1=jackB1$CI)
        result <- newMCResultJackknife( method.names = mnames, sample.names = sample.names, wdata=data, para=mc.res, regmeth=method.reg,
                                        cimeth=method.ci, B0jack=B0jack, B1jack=B1jack, alpha=alpha, glob.coef=glob.coef,
                                        error.ratio=error.ratio, weight=weight)
    }
    ## Bootstrap confidence intervals
    else if(method.ci %in% c("bootstrap", "nestedbootstrap")) 
    {
        ## Set random number generator
        rng.old <- NULL
        if(!is.null(rng.seed))
        {
            stopifnot(is.numeric(rng.seed))
            stopifnot(length(rng.seed)==1)
            if(exists(".Random.seed")) rng.old <- .Random.seed # store old RNG
            set.seed(seed=rng.seed,kind=rng.kind)
        }
        else
        {
            ## No defined RNG initialization => set parameters to NA for storage in object
            rng.seed <- as.numeric(NA)
            rng.kind <- as.character(NA)
        }
        
        #if (method.ci=="nestedbootstrap" & method.bootstrap.ci %in% c("Student","BCa")) warning("You need nested bootstrap only for 'tBoot' method.\n Please use option 'bootstrap' ")
        if(method.bootstrap.ci == "BCa") 
        {
            res <- mc.bootstrap(method.reg=method.reg, jackknife=TRUE, bootstrap=method.ci, X=data[,"x"], Y=data[,"y"], NBins=NBins, slope.measure=slope.measure,
                                error.ratio=error.ratio, nsamples=nsamples, threshold=threshold, iter.max=iter.max, nnested=nnested)
            B0jack <- res$B0jack
            B1jack <- res$B1jack
		} 
        else 
        {
            res <- mc.bootstrap(method.reg=method.reg, jackknife=FALSE, bootstrap=method.ci, X=data[,"x"], Y=data[,"y"], NBins=NBins, slope.measure=slope.measure,
                                error.ratio=error.ratio, nsamples=nsamples, threshold=threshold, iter.max=iter.max, nnested=nnested)
        }               
        xmean <- res$xmean
        B0 <- res$B0
        B1 <- res$B1
        MX <- res$MX
        sigmaB0 <- res$sigmaB0
        sigmaB1 <- res$sigmaB1
        glob.coef <- res$glob.coef
        glob.sigma <- res$glob.sigma
        npoints <- length(data[,"x"])
        nsamples <- res$nsamples
        nnested <- res$nnested
        weight <- res$weight
        
        if(method.bootstrap.ci == "Student") 
        {
            bootB0 <- mc.calc.Student(B0,xhat=glob.coef[1],alpha, npoints)
            bootB1 <- mc.calc.Student(B1,xhat=glob.coef[2],alpha, npoints)
        } 
        else if(method.bootstrap.ci == "BCa") 
        {
            bootB0 <- mc.calc.bca(Xboot=B0, Xjack=B0jack, xhat=glob.coef[1], alpha)
            bootB1 <- mc.calc.bca(Xboot=B1, Xjack=B1jack, xhat=glob.coef[2], alpha)
            bootB0$se <- glob.sigma[1]
            bootB1$se <- glob.sigma[2]
            bootB0$se <- NA
            bootB1$se <- NA
        } 
        else if(method.bootstrap.ci == "tBoot") 
        {
            bootB0 <- mc.calc.tboot(Xboot=B0, Sboot=sigmaB0, xhat=glob.coef[1], shat=glob.sigma[1], alpha)
            bootB1 <- mc.calc.tboot(Xboot=B1, Sboot=sigmaB1, xhat=glob.coef[2], shat=glob.sigma[2], alpha)
        }  
        else if(method.bootstrap.ci == "quantile") 
        {
            bootB0 <- mc.calc.quantile(Xboot=B0, alpha)
            bootB1 <- mc.calc.quantile(Xboot=B1, alpha)
            bootB0$se <- NA
            bootB1$se <- NA
        }
				
        if(method.reg=="WDeming") 
            cat("The global.sigma is calculated with Linnet's method\n")
		
        mc.res <- mc.make.CIframe(b0=glob.coef[1], b1=glob.coef[2], se.b0=bootB0$se, se.b1=bootB1$se, CI.b0=bootB0$CI, CI.b1=bootB1$CI)
        
        if(method.bootstrap.ci == "BCa")
        {
            result <- newMCResultBCa(method.names = mnames, sample.names = sample.names, 
                                     wdata=data, para=mc.res, xmean=xmean, regmeth=method.reg,
                                     cimeth=method.ci, bootcimeth=method.bootstrap.ci, B0jack=B0jack, B1jack=B1jack,
                                     B0=B0, B1=B1, MX=MX, sigmaB0=sigmaB0, sigmaB1=sigmaB1, alpha=alpha, glob.coef=glob.coef,
                                     glob.sigma=glob.sigma, nsamples=nsamples, nnested=nnested, error.ratio=error.ratio,
                                     rng.seed=rng.seed,rng.kind=rng.kind, weight=weight)
        }  
        else 
        {
            result <- newMCResultResampling(method.names = mnames, sample.names = sample.names, 
                                            wdata=data, para=mc.res, xmean=xmean, regmeth=method.reg,
                                            cimeth=method.ci, bootcimeth=method.bootstrap.ci,
                                            B0=B0, B1=B1, MX=MX, sigmaB0=sigmaB0, sigmaB1=sigmaB1, alpha=alpha, glob.coef=glob.coef,
                                            glob.sigma=glob.sigma, nsamples=nsamples, nnested=nnested, error.ratio=error.ratio,
                                            rng.seed=rng.seed,rng.kind=rng.kind, weight=weight)
        }
        
        ## Restore old random number generator
        if(!is.null(rng.old)) .Random.seed <- rng.old
    }

    ## Return MCResult object
    return(result)
}