##################################################################################

##' Class "NpdeRes"
##' 
##' The results component of a NpdeObject object
##' 
##' @name NpdeRes-class
##' @aliases NpdeRes NpdeRes-class, show,NpdeRes-method print,NpdeRes-method 
##' showall,NpdeRes-method summary,NpdeRes-method test,NpdeRes-method
##' [,NpdeRes-method [<-,NpdeRes-method
##' @docType class
##' @section Objects from the Class: NpdeRes objects are created during a call to   \code{\link{npde}} or \code{\link{autonpde}} as the "results" slot in a NpdeObject object. An NpdeRes object contains the following slots:
##' 
##' \describe{
##' \item{res}{a dataframe containing the results. Columns include id (group), xobs (observed X), yobs (observed Y), cens (indicator for censored data), as well as the actual results: ypred (model population predictions), pd (prediction discrepancies), npde (normalised prediction distribution errors), ycomp (completed data), ydobs (decorrelated observed data).}
##' \item{N}{number of subjects}
##' \item{ntot.obs}{total number of non-missing observations}
##' \item{ploq}{a vector giving the probability that a given observation is LOQ, according to the model}
##' \item{icens}{index of (non-missing) censored observations}
##' \item{not.miss}{a vector of boolean indicating for each observation whether it is missing (FALSE) or available (TRUE)}
##' \item{pd.sim}{pd computed for a number of simulated datasets (optional, used to obtain prediction intervals on the distribution of pd)}
##' \item{npde.sim}{npde computed for a number of simulated datasets (optional, used to obtain prediction intervals on the distribution of npde)}
##' }
##' @section Methods:
##' \describe{
##'   \item{print(npde.res):}{Prints a summary of object npde.res}
##'   \item{show(npde.res):}{Prints a short summary of object npde.res}
##'   \item{showall(npde.res):}{Prints a detailed summary of object npde.res}
##'   \item{plot(npde.res):}{Plots the data in npde.res. More details can be found in \code{\link{plot.NpdeRes}}}
##'   \item{summary(npde.res):}{Returns a summary of object npde.res in list format}
##' }
##' @seealso \code{\link{npde}}, \code{\link{autonpde}}, \code{\link{plot.NpdeRes}}, \code{\link{NpdeObject}}
##' @keywords classes internal
##' @examples
##' 
##' data(theopp)
##' 
##' methods(class="NpdeRes")
##' 
##' showClass("NpdeRes")
##' 
##' @exportClass NpdeRes

setClass(
  Class="NpdeRes",
  representation=representation(
  	res="data.frame",          # a data frame containing the results: id, xobs, yobs, cens (1 for censored data, 0 for observed data), ypred (predicted responses), pd (predicted responses), npde (npde), ycomp (completed responses), ydobs (decorrelated responses)
  	icens="numeric",		# index of the censored observations (non-missing)
    not.miss="logical",		# vector of logical, TRUE if present (=not missing), FALSE for missing data
    ploq="numeric",		# probability to be below LOQ
    xerr="numeric",		# error code
  	pd.sim="matrix",          # a matrix with pd for a sample of the simulated datasets (as many lines as non-missing observations)
  	npde.sim="matrix",          # a matrix with npde for a sample of the simulated datasets (as many lines as non-missing observations)
  	ntot.obs="numeric"		# total number of observations
  ),
  validity=function(object){
#    cat ("--- Checking NpdeRes object ---\n")
  	if(length(object@res)>0) {
  		if(length(object@not.miss)>0 && length(object@not.miss)!=dim(object@res)[1]) {
  			cat("Size mismatch between indicator variable not.miss and the data.\n")
  			return(FALSE)
  		}
  		if(length(object@ntot.obs)>0 && length(object@not.miss)>0 && object@ntot.obs!=dim(object@res[object@not.miss,])[1]) {
  			cat("Size mismatch between the total number of observations and the data.\n")
  			return(FALSE)
  		}
  	}
    return(TRUE)
  }
)

setMethod(
  f="initialize",
  signature="NpdeRes",
  definition= function (.Object){
#    cat ("--- initialising NpdeSimData Object --- \n")
    return (.Object )
  }
)

##################################################################################

##' Get/set methods for NpdeRes object
##' 
##' Access slots of a NpdeRes using the object["slot"] format
##' 
##' @keywords methods

#### NpdeRes
# Getteur
setMethod(
  f ="[",
  signature = "NpdeRes" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "res"={return(x@res)},
    "icens"={return(x@icens)},
    "not.miss"={return(x@not.miss)},
    "ploq"={return(x@ploq)},
    "xerr"={return(x@xerr)},    
  	"pd.sim"={return(x@pd.sim)},
  	"npde.sim"={return(x@npde.sim)},
    "ntot.obs"={return(x@ntot.obs)},
    stop("No such attribute\n")
   )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "NpdeRes" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "res"={x@res<-value},
    "icens"={x@icens<-value},
    "not.miss"={x@not.miss<-value},
    "ploq"={x@ploq<-value},
    "xerr"={x@xerr<-value},    
  	"pd.sim"={x@pd.sim<-value},
  	"npde.sim"={x@npde.sim<-value},  				
    "ntot.obs"={x@ntot.obs<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)

##################################################################################
# print/show/showall
# alias in class documentation

setMethod("print","NpdeRes",
					function(x,nlines=10,...) {
						digits<-2;nsmall<-2
						cat("Object of class NpdeRes\n")
						cat("    resulting from a call to npde or autonpde\n")
						if(length(x@res)[1]==0) cat("    currently empty\n") else {
							cat("    containing the following elements:\n")
							if(!is.null(x@res$ypred)) {
								cat("    predictions (ypred)\n")
								print(summary(x@res$ypred[x@not.miss]))
							}
							if(!is.null(x@res$pd)) {
								cat("    prediction discrepancies (pd)\n")
								print(summary(x@res$pd[x@not.miss]))
							}
							if(!is.null(x@res$npde)) {
								cat("    normalised prediction distribution errors (npde)\n")
								print(summary(x@res$npde[x@not.miss]))
							}
							if(!is.null(x@res$ycomp)) {
								cat("    completed responses (ycomp) for censored data\n")
							}
							if(!is.null(x@res$ydobs)) {
								cat("    decorrelated responses (ydobs)\n")
							}
							if(length(x@ploq)>0)  cat("    ploq: probability of being <LOQ for each observation\n")
							if(length(x@ntot.obs)) {
								cat("  the dataframe has ",x@ntot.obs,"non-missing observations ")
								if(dim(x@res)[1]>x@ntot.obs) cat("and",dim(x@res)[1],"lines")
								cat(".\n")
							}
							# ECO TODO: FINIR
							if(nlines!=0 & length(x@res)>0) {
								tab<-x@res[x@not.miss,]
								cat("First",nlines,"lines of results, removing missing observations:\n")
								if(nlines==(-1)) {
									print(tab)
								} else {
									nrowShow <- min (nlines , nrow(tab))
									print(tab[1:nrowShow,])
								}
							}
						}}
					)

# ECO TODO: FINIR, mettre le test ici ???
setMethod("show","NpdeRes",
					function(object) {
						cat("Object of class NpdeRes\n")
						if(length(object@res)[1]==0) cat("    currently empty\n") else {
							cat("  containing the following elements:\n")
							if(!is.null(object@res$ypred)) {
								cat("    predictions (ypred)\n")
							}
							if(!is.null(object@res$pd)) {
								cat("    prediction discrepancies (pd)\n")
							}
							if(!is.null(object@res$npde)) {
								cat("    normalised prediction distribution errors (npde)\n")
							}
							if(!is.null(object@res$ycomp)) {
								cat("    completed responses (ycomp) for censored data\n")
							}
							if(!is.null(object@res$ydobs)) {
								cat("    decorrelated responses (ydobs)\n")
							}
							if(length(object@ploq)>0)  cat("    ploq: probability of being <LOQ for each observation\n")
							if(length(object@ntot.obs)) {
								cat("  the dataframe has ",object@ntot.obs,"non-missing observations ")
								if(dim(object@res)[1]>object@ntot.obs) cat("and",dim(object@res)[1],"lines")
								cat(".\n")
							}
						}}
					)

setMethod("showall","NpdeRes",
					function(object) {
						cat("Object of class NpdeRes\n")
						if(length(object@res)[1]==0) cat("    currently empty\n") else {
							cat("    containing the following elements:\n")
							if(!is.null(object@res$ypred)) {
								cat("    predictions (ypred)\n")
							}
							if(!is.null(object@res$pd)) {
								cat("    prediction discrepancies (pd)\n")
							}
							if(!is.null(object@res$npde)) {
								cat("    normalised prediction distribution errors (npde)\n")
							}
							if(!is.null(object@res$ycomp)) {
								cat("    completed responses (ycomp) for censored data\n")
							}
							if(!is.null(object@res$ydobs)) {
								cat("    decorrelated responses (ydobs)\n")
							}
							if(length(object@ploq)>0)  cat("    ploq: probability of being <LOQ for each observation\n")
							if(length(object@ntot.obs)) {
								cat("  the dataframe has ",object@ntot.obs,"non-missing observations ")
								if(dim(object@res)[1]>object@ntot.obs) cat("and",dim(object@res)[1],"lines")
								cat(".\n")
							}
							nlines<-10
							cat("First",nlines,"lines of results, removing missing observations:\n")
							tab<-object@res[object@not.miss,]
							nrowShow <- min (nlines , nrow(tab))
							print(tab[1:nrowShow,])
						}
					}
					)

setMethod("summary","NpdeRes",
	function(object, print=TRUE, ...) {
		if(length(object@res)==0) {
			cat("Object of class NpdeRes, empty.\n")
			return()
		}
		res<-list(N=object@N,data=object@res,ntot.obs=object@ntot.obs)
		res$ploq<-object@ploq
		invisible(res)
		}
)

##################################################################################
#
#' Plots a NpdeRes object
#'
#' Plots distribution and scatterplots for the npde in a NpdeRes object. Users are advised to use the plot() function on the NpdeObject object resulting from a call to npde() or autonpde() instead of trying to plot only the results element of this object.
#' 
#' @param x a NpdeRes object
#' @details Four graphs are produced:
#' \describe{
#' \item{a quantile-quantile plot}{plot of the npde versus the corresponding quantiles of a normal distribution, with the line y=x overlayed.}
#' \item{a histogram of the npde}{the shape of the normal distribution is also shown}
#' \item{two scatterplots of the npde}{a plot of the npde versus the independent variable X and a plot of the npde versus the empirical mean of the predicted distribution; for these last two graphs, we plot the lines corresponding to y=0 	and to the 5% and 95% critical value of the normal distribution delimiting a 90% prediction interval for the npde}
#' }
#' @references K. Brendel, E. Comets, C. Laffont, C. Laveille, and F.Mentre. Metrics for external model evaluation with an application to the population pharmacokinetics of gliclazide. \emph{Pharmaceutical Research}, 23:2036--49, 2006.
##' @seealso \code{\link{set.plotoptions}}
#' @keywords plot internal
#' @examples
#' 
##' data(theopp)
##' 
##' @importFrom graphics plot
##' @method plot NpdeRes
##' @export

#setMethod("plot","NpdeRes",
plot.NpdeRes<-
					function(x, y, ...) {
						if(length(x@res)==0) {
							cat("No data to plot.\n")
							return()
						} else xres<-x@res[x@not.miss,]
						if(length(xres$npde)>0) {
							nclass<-10
							npde<-xres$npde
							par(mfrow=c(2,2))
							qqnorm(sort(npde),xlab="Sample quantiles (npde)",ylab="Theoretical Quantiles",
										 cex.lab=1.5,main="Q-Q plot versus N(0,1) for npde")
							qqline(sort(npde))
							#Histogram of npde, with N(0,1) superimposed on the plot
							xh<-hist(npde,nclass=nclass,xlab="npde",main="",cex.lab=1.5)
							xpl<-min(npde)+c(0:100)/100*(max(npde)-min(npde))
							ypl<-dnorm(xpl)
							ypl<-ypl/max(ypl)*max(xh$counts)
							lines(xpl,ypl,lwd=2)
							
							#residuals
							plot(xres$xobs,npde,xlab="X",ylab="npde",cex.lab=1.5)
							abline(h=0,lty=2)
							x1<-qnorm(0.05)
							abline(h=x1,lty=3);abline(h=(-x1),lty=3)
							
							plot(xres$ypred,npde,xlab="Predicted Y",ylab="npde",cex.lab=1.5)
							abline(h=0,lty=2)
							abline(h=x1,lty=3);abline(h=(-x1),lty=3)    
						}
						if(length(xres$pd)>0) {
							pd<-xres$pd
							nclass<-10
							par(mfrow=c(2,2))
							samp<-sort(pd);ndat<-length(samp)
							theo<-c(1:ndat)/ndat
							qqplot(samp,theo,xlab="Sample quantiles (pd)",ylab="Theoretical Quantiles",
										 cex.lab=1.5,main="Q-Q plot versus U(0,1) for pd")
							segments(0,0,1,1)
							#Histogram of pd, with N(0,1) superimposed on the plot
							xh<-hist(pd,nclass=nclass,xlab="pd",main="",cex.lab=1.5)
							abline(h=ndat/nclass,lty=2,lwd=2)
							
							#residuals
							plot(xres$xobs,pd,xlab="X",ylab="pd",cex.lab=1.5)
							abline(h=0.5,lty=2)
							abline(h=0.05,lty=3);abline(h=0.95,lty=3)
							plot(xres$ypred,pd,xlab="Predicted Y",ylab="pd",cex.lab=1.5)
							abline(h=0.5,lty=2)
							abline(h=0.05,lty=3);abline(h=0.95,lty=3)
						}
					}
#)

##################################################################################

#' Kurtosis
#' 
#' Computes the kurtosis.
#' 
#' If \eqn{N = \mathrm{length}(x)}{N = length(x)}, then the kurtosis of \eqn{x}
#' is defined as: \deqn{N sum_i (x_i-\mathrm{mean}(x))^4 (sum_i
#' (x_i-\mathrm{mean}(x))^2)^(-2) - }{N sum_i (x_i-mean(x))^4 (sum_i
#' (x_i-mean(x))^2)^(-2) - 3}\deqn{3}{N sum_i (x_i-mean(x))^4 (sum_i
#' (x_i-mean(x))^2)^(-2) - 3}
#' 
#' @usage kurtosis(x)
#' @param x a numeric vector containing the values whose kurtosis is to be
#' computed. NA values are removed in the computation.
#' @return The kurtosis of \code{x}.
#' @references G. Snedecor, W. Cochran. \emph{Statistical Methods}, 
#' Wiley-Blackwell, 1989
#' @keywords univar
#' @examples
#' 
#' x <- rnorm(100)
#' kurtosis(x)
#' @export

kurtosis<-function (x) 
{
#from Snedecor and Cochran, p 80
    x<-x[!is.na(x)]
    m4<-sum((x - mean(x))^4)
    m2<-sum((x - mean(x))^2)
    kurt<-m4*length(x)/(m2**2)-3
    return(kurtosis=kurt)
}

#' Skewness
#' 
#' Computes the skewness.
#' 
#' If \eqn{N = \mathrm{length}(x)}{N = length(x)}, then the skewness of \eqn{x}
#' is defined as \deqn{N^{-1} \mathrm{sd}(x)^{-3} \sum_i (x_i -
#' \mathrm{mean}(x))^3.}{ N^(-1) sd(x)^(-3) sum_i (x_i - mean(x))^3.}
#' 
#' @usage skewness(x)
#' @param x a numeric vector containing the values whose skewness is to be
#' computed. NA values are removed in the computation.
#' @return The skewness of \code{x}.
#' @references G. Snedecor, W. Cochran. \emph{Statistical Methods}, 
#' Wiley-Blackwell, 1989
#' @keywords univar
#' @examples
#' 
#' x <- rnorm(100)
#' skewness(x)
#' 
#' @export

skewness<-function (x) 
{
#from Snedecor and Cochran, p 79
    x<-x[!is.na(x)]
    m3<-sum((x - mean(x))^3)
    m2<-sum((x - mean(x))^2)
    skew<-m3/(m2*sqrt(m2/length(x)))
    return(skewness=skew)
}


##' @S3method gof.test numeric
#' @export

gof.test.numeric<-function(object,which="npde",parametric=TRUE, ...) {
# Default is to compare the distribution in object to N(0,1)	
	args1<-match.call(expand.dots=TRUE)
	verbose<-TRUE
	i1<-match("verbose",names(args1))
	if(!is.na(i1) && !is.na(as.logical(as.character(args1[[i1]])))) verbose<-as.logical(as.character(args1[[i1]]))
	object<-object[!is.na(object)]
	sev<-var(object)*sqrt(2/(length(object)-1))
	semp<-sd(object)
	n1<-length(object)
	sem<-semp/sqrt(length(object))
	res<-list(mean=mean(object),se.mean=sem,var=var(object),se.var=sev, kurtosis=kurtosis(object),skewness=skewness(object))
	object<-object[!is.na(object)]
	myres<-rep(0,4)
	if(parametric) {
		if(which=="pd") y<-t.test(object,mu=0.5) else y<-t.test(object)
	} else {
		if(which=="pd") y<-wilcox.test(object,mu=0.5) else y<-wilcox.test(object)
	}
	myres[1]<-y$p.val
# ECO TODO: ici utiliser le test de Magalie
	if(which=="pd") y<-ks.test(object,"punif",min=min(object,na.rm=TRUE), max=max(object,na.rm=TRUE)) else y<-shapiro.test(object)
	myres[3]<-y$p.val
	
# test de variance pour 1 échantillon
# chi=s2*(n-1)/sigma0 et test de H0={s=sigma0} vs chi2 à n-1 df
#    if(parametric) {
	chi<-(semp**2)*(n1-1)
	if(which=="pd") chi<-chi*12 # X~U(0,1) => var(X)=1/12
	y<-2*min(pchisq(chi,n1-1),1-pchisq(chi,n1-1))
	myres[2]<-y
#    } else {
# ECO TODO: non-parametric equivalent of variance test for one-sample ?
#    }
	xcal<-3*min(myres[1:3])
	myres[4]<-min(1,xcal)
	if(parametric) 
		names(myres)<-c("  t-test                    ","  Fisher variance test      ","  SW test of normality      ", "Global adjusted p-value     ") else 
			names(myres)<-c("  Wilcoxon signed rank test ","  Fisher variance test      ", "  SW test of normality      ","Global adjusted p-value     ")
	if(which=="pd") names(myres)[3]<-"KS test of uniformity       "
	res$p.value<-myres
	res$nobs<-n1
#	if(verbose) print.gof.test(res, which=which)
	invisible(res)
}

##' @S3method gof.test NpdeRes
#' @export

# Performs test on the selected variable (which=one of npde, pd or npd)
### parametric=TRUE: use parametric tests for mean (t.test) and variance (Fisher)
### if na.omit (default), missing values are removed before testing

# test on npde, npd:
### mean=0
### variance=1
### normality (Shapiro-Wilks)
# test on pd:
### mean=0.5
### variance=1/12
### uniformity (Kolmogorov-Smirnov)

# ECO TODO: non-parametric equivalent of variance test for one-sample ?
### ECO TODO: test for pd pd~U(0,1) using test from Magalie

gof.test.NpdeRes<-function(object,which="npde",parametric=TRUE, ...) {
# Performs test on the selected variable (one of npde, pd or npd)
	args1<-match.call(expand.dots=TRUE)
	if(length(object@res)==0) {
		cat("No data available.\n")
		return()
	}
	if(!which%in%c("pd","npde","npd")) {
		cat("Tests can be performed on one of: npde (default), pd, npd. Please choose one using the which argument.\n")
		return()
	}
	if(which=="npde" & length(object@res$npde)==0) {
		cat("    Missing npde object to plot.\n")
		return()
	}
	if(which%in%c("pd","npd") & length(object@res$pd)==0) {
		cat("    Missing pd object to plot.\n")
		return()
	}
	npde<-switch(which,npde=object@res$npde,pd=object@res$pd, npd=qnorm(object@res$pd))
	npde<-npde[object@not.miss] # Removing values for MDV=1 (pd, npde not computed)
	verbose<-TRUE
	i1<-match("verbose",names(args1))
	if(!is.na(i1) && !is.na(as.logical(as.character(args1[[i1]])))) verbose<-as.logical(as.character(args1[[i1]]))
	args1<-match.call(expand.dots=TRUE)
	i1<-match("na.action",names(args1))
	na.action<-"na.omit"
		if(!is.na(i1) && as.character(args1[[i1]]) %in% c("na.omit","na.fail","na.exclude","na.pass")) na.action<-as.character(args1[[i1]])
		if(na.action=="na.fail" & sum(is.na(npde))>0) {
			cat("Missing values and na.action is set to na.fail.\n")
			return()
		}
		if(na.action=="na.pass" & sum(is.na(npde))>0) {
			if(verbose) cat("Warning: there are missing values and na.action is set to na.pass. Results and tests will be obtained removing the missing data (nmis=", sum(is.na(npde)),").\n")
			}
	npde<-eval(call(na.action,npde))
	res<-gof.test(npde,which=which,parametric=parametric)
	if(verbose)
		print.gof.test(res, which=which)
	invisible(res)
}

#' Prints a summary of a gof.test result

print.gof.test<-function(object, which="npde", ...) {
	cat("---------------------------------------------\n")
	cat("Distribution of",which,":\n")
	cat("      nb of obs:",object$nobs,"\n")
	cat("           mean=",format(object$mean,digits=4),"  (SE=",format(object$se.mean,digits=2),")\n")
	cat("       variance=",format(object$var,digits=4),"  (SE=",format(object$se.var,digits=2),")\n")
	cat("       skewness=",format(object$skewness,digits=4),"\n")
	cat("       kurtosis=",format(object$kurtosis,digits=4),"\n")
	cat("---------------------------------------------\n\n")
	cat("Statistical tests\n")
	for(i in 1:4) {
		cat(names(object$p.value)[i],": ")
		#if (myres[i]<1) 
		cat(format(object$p.value[i],digits=3)) 
		#else cat(myres[i])
		if(as.numeric(object$p.value[i])<0.1 & as.numeric(object$p.value[i])>=0.05) 
			cat(" .")
		if(as.numeric(object$p.value[i])<0.05) cat(" *")
		if(as.numeric(object$p.value[i])<0.01) cat("*")
		if(as.numeric(object$p.value[i])<0.001) cat("*")
		cat("\n")
	}
	cat("---\n")
	cat("Signif. codes: '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 \n")
	cat("---------------------------------------------\n")
}
##################################################################################
