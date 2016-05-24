# TODO: Add comment
# 
# Author: ecor
###############################################################################

###  http://cran.r-project.org/web/views/Distributions.html

rm(list=ls())
library(MASS)
library(RMAWGEN)
##library(RHEG)
library(fitdistrplus)
source('~/Dropbox/iasma/RMAWGENdev/RHEG/R/HEG.R', chdir = TRUE)
### additional function 
## fitting with a parametric probability distribution
qqfit <- function(x,FUN=qexp,par=list()) {

	if (!is.list(par)) {
		names <- names(par)
		par <- as.list(par)
		names(par) <- names
		
	}
	
	
	pval <- ecdf(x)(x)
	args <- par
	args[["p"]] <- as.vector(pval)
	qval <- do.call(what=FUN,args=args)
	return(qval)
	
}

## modified signature for 'ks.test'  or other test function suitably to run this script
gof.test.mod <- function(x,y="pexp",what="ks.test",par=list(),...){
	x <- sort(x)
	if (is.numeric(y)) y <- sort(y)
	if (!is.list(par)) {
		names <- names(par)
		par <- as.list(par)
		names(par) <- names
		
	}
	
	args <- list(...)
	if (length(par)>1) args[names(par)] <- par 
	if (length(par)==1) args[[names(par)]] <- par[[1]] 
	
	args[["y"]] <- y 
	args[["x"]] <- x
	
	htest <- do.call(what,args)
	return(htest)
}



data(trentino)


year_min <- 1961
year_max <- 1990


period <- PRECIPITATION$year>=year_min & PRECIPITATION$year<=year_max
season <- PRECIPITATION$month %in% c(3,4,5) ## SPRING!
station <- names(PRECIPITATION)[!(names(PRECIPITATION) %in% c("day","month","year"))]
 
prec_mes <- PRECIPITATION[period,station]   
season <- season[period]
## removing nonworking stations (e.g. time series with NA)

accepted <- array(TRUE,length(names(prec_mes)))
names(accepted) <- names(prec_mes)
for (it in names(prec_mes)) {
	accepted[it]  <- (length(which(!is.na(prec_mes[,it])))==length(prec_mes[,it]))
}


prec_mes <- prec_mes[season,accepted]		

stations <- names(prec_mes)

season <- PRECIPITATION$month %in% c(3,4,5) ## SPRING!

## List of probability distributions to test for precipitation.
distributions <- list(exp=dexp,gamma=dgamma,weibull=dweibull) ###,heg=dheg)


##distributions <- list(exp=dexp)  ##,gpd=dgpd,heg=dheg)
## Initialization Value For each tested probability distribution

distributionPar <- lapply(X=names(distributions),FUN=function(x) {
		 dfun <- get(paste("d",x,sep=""))
		 out <- formals(dfun)
		 out <- out[!(names(out) %in% c("x","log"))]
		 for (f in 1:length(out)) {
			 
			 if (class(out[[f]])!="numeric") out[[f]] <- 1
		 }
		
		 out$densfun <- dfun
		 out$namefun <- x
		 
		 return(out)
})
names(distributionPar) <- names(distributions)


## Specific correction on density function arguments 
distributionPar$gamma <- distributionPar$gamma[names(distributionPar$gamma)!="scale"]



distributionParForEachStation <- list()

for (it in stations[2:3]) {
	data <- prec_mes[,it]
	data <- data[data>0]
	
	distributionParForEachStation[[it]] <-  lapply(X=distributionPar,FUN=function(x,data){
				
			 start <- x[!(names(x) %in% c("densfun","namefun"))]
	
		     out <- fitdistr(data,x$densfun,start=start,method=method) ##,method="mle") ###,method="CG") ###,fix.arg=args.fix)
			 return(out)
				
			},data=data)				


###fitdistr(x,dexp,start=start[names(start) %in% names(formals(dexp))])
	
	
	
}







stop("stop here")

prec_exp <- prec_mes
prec_gamma <- prec_mes
prec_weibull <- prec_mes


list_exp <- list()         ## dexp
list_gamma <- list()       ## dgamma 
list_weibull <- list()     ## dweibull

for (it in names(prec_mes)) {
	
	x <- prec_mes[,it]
	x <- x[x>0]
	start <- list(shape=1,scale=1,rate=1,beta=1,xi=1)
	
	if (length(x)>0) {
		list_exp[[it]] <- list()
		list_gamma[[it]] <- list()
		list_weibull[[it]] <- list()
		list_gpd <- list()
		list_heg[[it]] <- list()
		
		
		list_exp[[it]]$par <- fitdistr(x,dexp,start=start[names(start) %in% names(formals(dexp))])
		list_gamma[[it]]$par <- fitdistr(x,dgamma,start=start[names(start) %in% names(formals(dgamma))][1:2])
		list_weibull[[it]]$par <- fitdistr(x,dweibull,start=start[names(start) %in% names(formals(dweibull))])
		list_weibull[[it]]$par <- fitdistr(x,dweibull,start=start[names(start) %in% names(formals(dweibull))])
		
		# Kolgomorov-Smirnov tests 
		
		list_exp[[it]]$htest <- gof.test.mod(x=x,y="pexp",par=list_exp[[it]]$par$estimate)
		list_gamma[[it]]$htest <- gof.test.mod(x=x,y="pgamma",par=list_gamma[[it]]$par$estimate)
		list_weibull[[it]]$htest <- gof.test.mod(x=x,y="pweibull",par=list_weibull[[it]]$par$estimate)
		
		### Theoretical Quantile Precipitation 
		
		prec_exp[,it][prec_mes[,it]>0] <- qqfit(prec_mes[,it][prec_mes[,it]>0],FUN=qexp,par=list_exp[[it]]$par$estimate)
		prec_gamma[,it][prec_mes[,it]>0] <- qqfit(prec_mes[,it][prec_mes[,it]>0],FUN=qgamma,par=list_gamma[[it]]$par$estimate)
		prec_weibull[,it][prec_mes[,it]>0] <- qqfit(prec_mes[,it][prec_mes[,it]>0],FUN=qweibull,par=list_weibull[[it]]$par$estimate)

		x <- prec_exp[,it][prec_mes[,it]>0] 
		list_exp[[it]]$htest0 <- gof.test.mod(x=x,y="pexp",par=list_exp[[it]]$par$estimate)
		x <- prec_gamma[,it][prec_mes[,it]>0] 
		list_gamma[[it]]$htest0 <- gof.test.mod(x=x,y="pgamma",par=list_gamma[[it]]$par$estimate)
		x <- prec_weibull[,it][prec_mes[,it]>0] 
		list_weibull[[it]]$htest0 <- gof.test.mod(x=x,y="pweibull",par=list_weibull[[it]]$par$estimate)
		
	} else {
	
		list_exp[it] <- NULL
		list_gamma[it] <- NULL
		list_weibull[it] <- NULL
		
	}
	
	
	
	
	
}

## resume p-Values 

pvalues_exp <- unlist(lapply(X=list_exp,FUN=function(x) {x$htest$p.value}))
pvalues_gamma <- unlist(lapply(X=list_gamma,FUN=function(x) {x$htest$p.value}))
pvalues_weibull <- unlist(lapply(X=list_weibull,FUN=function(x) {x$htest$p.value}))

pvalues0_exp <- unlist(lapply(X=list_exp,FUN=function(x) {x$htest0$p.value}))
pvalues0_gamma <- unlist(lapply(X=list_gamma,FUN=function(x) {x$htest0$p.value}))
pvalues0_weibull <- unlist(lapply(X=list_weibull,FUN=function(x) {x$htest0$p.value}))



## scratch 

#		val <- prec_mes[,it]
#		probs <- pexp(val,rate=list_exp[[it]]["rate"])
#		val[prec_mes[,it]>0] <- quantile(val[prec_mes[,it]>0],probs=probs[prec_mes[,it]>0])
#		prec_exp[,it] <- val
#		val <- prec_mes[,it]
#		pval <- ecdf(val[prec_mes[,it]>0])(val[prec_mes[,it]>0])
#		qval <- qexp(pval,rate=list_exp[[it]]["rate"])
#		val[prec_mes[,it]>0] <- qval
#		prec_exp[,it] <- val





