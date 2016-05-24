
#' Compares the rates of two variables against one another.
#' @param first.interval An \code{rds.interval.estimate} object fit with either "Gile" or "Salganik" uncertainty.
#' @param second.interval An \code{rds.interval.estimate} object fit with either "Gile" or "Salganik" uncertainty.
#' @param M The number of bootstrap resamplings to use
#' @details This function preforms a bootstrap test comparing the 
#' the rates of two variables against one another.
#' @examples
#' \dontrun{
#' data(faux)
#' int1 <- RDS.bootstrap.intervals(faux, outcome.variable=c("X"), 
#'	weight.type="RDS-II", uncertainty="Salganik", N=1000,
#'	number.ss.samples.per.iteration=1000, 
#' 	confidence.level=0.95, number.of.bootstrap.samples=100)
#' int2 <- RDS.bootstrap.intervals(faux, outcome.variable=c("Y"), 
#' 	weight.type="RDS-II", uncertainty="Salganik", N=1000,
#'	number.ss.samples.per.iteration=1000,
#'	confidence.level=0.95, number.of.bootstrap.samples=100)
#' RDS.compare.proportions(int1,int2)
#' }
#' @export
RDS.compare.proportions <- function(first.interval,second.interval,M=10000){
	
	bsresult1 <- attr(first.interval,"bsresult") 
	boots1 <- if(is.matrix(bsresult1)) bsresult1 else bsresult1$bsests
	levs1 <- colnames(boots1)
	nlev1 <- length(levs1)
	bs1 <- boots1[sample(1:nrow(boots1),replace=TRUE,size=M),]
	tmp <- first.interval$estimate
	obs <- tmp[match(names(tmp),levs1)]
	bs1 <- sweep(bs1,2,colMeans(bs1)-obs)

	bsresult2 <- attr(second.interval,"bsresult") 
	boots2 <- if(is.matrix(bsresult2)) bsresult2 else bsresult2$bsests
	levs2 <- colnames(boots2)
	nlev2 <- length(levs2)
	bs2 <- boots2[sample(1:nrow(boots2),replace=TRUE,size=M),]
	tmp <- second.interval$estimate
	obs <- tmp[match(names(tmp),levs2)]
	bs2 <- sweep(bs2,2,colMeans(bs2)-obs)

	res <- matrix(NA,ncol=nlev2,nrow=nlev1)
	colnames(res) <- levs2
	rownames(res) <- levs1
	for(i in 1:nlev1){
		for(j in 1:nlev2){
			bs <- bs1[,i] - bs2[,j]
			pval <- 2*min(c(mean(bs<=0),mean(bs>=0)))
			res[i,j] <- pval
		}
	}
	res <- as.data.frame(res)
	class(res) <- c("pvalue.table","data.frame")
	res
}


#' Compares the rates of two variables against one another.
#' @param data An object of class \code{rds.interval.estimates.list} with attribute \cr
#' \code{variables} containing a character vector of names of objects of class\cr
#' \code{rds.interval.estimate}.
#' @param variables A character vector of column names to select from \code{data}.
#' @param confidence.level The confidence level for the confidence intervals. The default is 0.95 for 95\%.
#' @param number.of.bootstrap.samples The number of Monte Carlo draws to determine the null distribution of the likelihood
#' ratio statistic.
#' @param plot Logical, if TRUE then a plot is produces of the null distribution of the likelihood
#' ratio statistic with the observed statistics plotted as a vertical dashed line.
#' @param seed The value of the random number seed. Preset by default to allow reproducability.
#' @return An object of class \code{pvalue.table} containing the cross-tabulation of p-values
#' for comparing the two classes
#' @export
RDS.compare.two.proportions <- function(data,variables,confidence.level=0.95,
   number.of.bootstrap.samples=5000,plot=FALSE,seed=1){
    vars <- attr(data,"variables")
    vars <- vars[match(variables,colnames(data),nomatch=NULL)]
	first.interval <- vars[1][[1]]
	bsresult1 <- attr(first.interval,"bsresult")
	boots1 <- if(is.matrix(bsresult1)) bsresult1 else bsresult1$bsests
	levs1 <- colnames(boots1)
	nlev1 <- length(levs1)
	bs1 <-boots1[sample(1:nrow(boots1),replace=TRUE,size=number.of.bootstrap.samples),]
	tmp <- first.interval$estimate
	obs <- tmp[match(names(tmp),levs1)]
	bs1 <- sweep(bs1,2,colMeans(bs1)-obs)

	second.interval <- vars[2][[1]]
	bsresult2 <- attr(second.interval,"bsresult") 
	boots2 <- if(is.matrix(bsresult2)) bsresult2 else bsresult2$bsests
	levs2 <- colnames(boots2)
	nlev2 <- length(levs2)
	bs2 <- boots2[sample(1:nrow(boots2),replace=TRUE,size=number.of.bootstrap.samples),]
	tmp <- second.interval$estimate
	obs <- tmp[match(names(tmp),levs2)]
	bs2 <- sweep(bs2,2,colMeans(bs2)-obs)

	res <- matrix(NA,ncol=nlev2,nrow=nlev1)
	colnames(res) <- levs2
	rownames(res) <- levs1
	for(i in 1:nlev1){
		for(j in 1:nlev2){
			bs <- bs2[,i] - bs1[,j]
			pval <- 2*min(c(mean(bs<=0),mean(bs>=0)))
			res[i,j] <- pval
		}
	}
	res <- as.data.frame(res)
	class(res) <- c("pvalue.table","data.frame")
#
    x <- as.numeric( sapply(vars,function(x) {x$interval[2]}) )
    sigma <- as.numeric( sapply(vars,function(x) {x$interval[10]}) )
    if ('par' %in% plot) {par(mfrow = c(1,2))} else {par(mfrow = c(1,1))}

    #  obsL <- diff(x)
    obsL <- 0
    
    if("distributions" %in% plot){
    binn <- 500
    maxl <- 1.1*max(obsL,stats::quantile(bs2[,2]-bs1[,2],0.99))
    minl <- (1/1.1)*min(obsL,stats::quantile(bs2[,2]-bs1[,2],0.01))
    r <- seq(minl, maxl, length = binn + 1)[-1] - 0.5*(maxl-minl)/binn
#   yl <- locfit::locfit.raw(locfit::lp(bs2[,2]-bs1[,2],nn=0.1,h=0.8), xlim=c(minl,maxl))
#   gpdf <- predict(yl, newdata=r)
    Range=range(bs2[,2]-bs1[,2])
    a=bgk_kde(bs2[,2]-bs1[,2],n=2^(ceiling(log(maxl-minl)/log(2))),
	        MIN=1.1*Range[1]-0.1*Range[2],
	        MAX=1.1*Range[2]-0.1*Range[1])
    gpdf <- stats::spline(x=a[1,],y=a[2,],xout=r)$y
    scalef <- binn/sum(gpdf)
    gpdf <- gpdf * scalef
    #   maxl <- r[which.max(cumsum(gpdf)>binn*0.99)]
    plot(x=r,y=gpdf,xlim=c(minl,maxl), type="l", ylab="Density", xlab="Difference in Proportions",main="Distribution of the Difference in Proportions",sub="The vertical line is the hypothesis of no difference in proportions")
    abline(v=obsL,lty=2)
    }
    
    if("estimates" %in% plot){
     yminus <- x - 1.96*sigma
     yplus  <- x + 1.96*sigma
     errorbar(x=seq_along(x),y=x,yminus=yminus, yplus=yplus,xlab="sequence",ylab="estimate", main="Trend of Estimates")
    }
    

    
  cat(sprintf("The p-value of a difference between the proportions is %s\n",format.pval(res[2,2])))
  if(res[2,2] < (1-confidence.level)){
    cat(sprintf("The hypothesis of equal proportions is rejected (at the %s%% level).\n",format(100*(1-confidence.level))))
  }else{
    cat(sprintf("The hypothesis of equal proportions is not rejected (at the %s%% level).\n",format(100*(1-confidence.level))))
  }
  
   invisible(res)
}

#' Displays a pvalue.table
#' @param x a pvalue.table object
#' @param ... additional parameters passed to print.data.frame.
#' @export
#' @method print pvalue.table
print.pvalue.table <- function(x,...){
	cat("P-Value table comparing rates:\n")
	print.data.frame(x,...)
}
