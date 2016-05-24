NULL
#' generate
#' 
#' Implementation of \code{generate} method for \code{YuleWalkerCoefficientBlockmatrices} or \code{CCGammaObject} S3 object.  It generates a multivarite random series according using a VAR model with coefficient obtained by \code{\link{CoeffYWeq}} (\code{YuleWalkerCoefficientBlockmatrices} S3 object) . Alternatively it generates by applying a first-order Markov Cain from \code{CCGammaObject} S3 Object.
#'
#' @param x \code{YuleWalkerCoefficientBlockmatrices} S3 object. See \code{\link{CoeffYWeq}}
#' @param FUN random function of the probability distribution used for noise random generation. Default is \code{\link{rnorm}}. See \url{http://cran.r-project.org/web/views/Distributions.html} 
#' @param n number of generations requested 
#' @param names null object or string vectors or names of the variables to be generated simultaneously. Default is \code{NULL}.
# @param K number of the variables to be generated simultaneously, i.e. the K parameters of a VAR. It is automatically detected by \code{x}, \code{names} or \code{cov}, if one of these is not \code{NULL}. 
##### ### @param cov null object or covariance matrix of the random variables to be generated simultaneously. Default is \code{NULL}, not used in case this information can be detected from \code{x}.
#' @param xprev null object or initial condition of the multivariate random process to be generated. Default is \code{NULL}. 
#' @param names_x names of the elements of a \code{YuleWalkerCoefficientBlockmatrices} S3 object. See examples. 
#' @param nearPD logical. If \code{TRUE} (Default) the function \code{\link{nearPD}} is applied to the covariance matrix of residuals. If \code{FALSE}, the estimated covariance matrix is not verified to be positive definite and error may occur during the function execution.  
#' @param precipitation.indicator logical value. Default is \code{FALSE}. If it is  \code{TRUE}, the output is transformed to a number between 0 (no precipitation) and 1 (high precipitation) taking into account the probabability of no occurence and the cumulate probability function referred to \code{FUN}.  
#' @param year_min,year_max first and last years of generation period
#' @param ... additional arguments for \code{FUN}
#' 
#' @return a matrix or a data frame object
#' 
#' 
#' 
#' @title generate
# @name generate
# @rdname generate
#' @rdname generate
#' @method generate YuleWalkerCoefficientBlockmatrices
#' @S3method generate YuleWalkerCoefficientBlockmatrices
#' @aliases generate generate.YuleWalkerCoefficientBlockmatrices 
#' @importFrom RGENERATE generate
#' @export
# @import methods
#' @seealso \code{\link{CoeffYWeq}},\code{\link{CCGammaToBlockmatrix}}
#' 
#' @examples 
#' 
#' library(RMRAINGEN)
#' 
#' set.seed(125)
#' data(trentino)
#' 
#' year_min <- 1961
#' year_max <- 1990
#' 
#' period <- PRECIPITATION$year>=year_min & PRECIPITATION$year<=year_max
#' station <- names(PRECIPITATION)[!(names(PRECIPITATION) %in% c("day","month","year"))]
#' prec_mes <- PRECIPITATION[period,station]  
#' 
#' ## removing nonworking stations (e.g. time series with NA)
#' accepted <- array(TRUE,length(names(prec_mes)))
#' names(accepted) <- names(prec_mes)
#' for (it in names(prec_mes)) {
#' 		 accepted[it]  <- (length(which(!is.na(prec_mes[,it])))==length(prec_mes[,it]))
#' }
#'
#' prec_mes <- prec_mes[,accepted]
#' ## the dateset is reduced!!! 
#' prec_mes <- prec_mes[,1:2]
#' 
#' ## Not Run in the examples, uncomment to run the following lines 
#' # coeff <- CoeffYWeq(data=prec_mes,p=1,tolerance=0.001)
#' 
#' # generation <- generate(coeff,n=10,names=names(prec_mes))
#' 
#'
#


#' ## Not Run in the examples, uncomment to run the following lines 
#' # origin <- paste(year_min,1,1,sep="-")
#' 
#'
#' # coeff_monthly <- CoeffYWeq(data=prec_mes,p=1,tolerance=0.001,sample="monthly",origin=origin)
#' 
#' 
#' # generation_monthly <- generate(coeff_monthly,year_min=year_min,year_max=year_max,
#' #					names=names(prec_mes))
#' 
#' 
#' ###  generation with CCGammaObject
#' 
#' # CCGamma <- CCGamma(data=prec_mes,lag=0,tolerance=0.001,only.matrix=FALSE)
#' 
#' # generation_CCGamma <- generate(x=CCGamma,n=100,names=names(prec_mes)) 
#' 
#' # CCGamma_monthly <- CCGamma(data=prec_mes,lag=0,tolerance=0.001,only.matrix=FALSE,
#' #                            sample="monthly",origin=origin)
#' ## generation_CCGamma <- generate(x=CCGamma_monthly,year_min=year_min,year_max=year_max,
#' ##                                names=names(prec_mes)) 
#' 




generate.YuleWalkerCoefficientBlockmatrices <- function(x,FUN=rnorm,n=100,names=NULL,xprev=NULL,names_x=c("A","Sigma_u","CCGammaInfo"),nearPD=TRUE,precipitation.indicator=FALSE,...) {
	
	out <- NULL 
	
	cov <- x[[names_x[2]]][1,1]
	
	if (nearPD) cov <- nearPD(cov)$mat 
	
	A <- x[[names_x[1]]] 
	
	p <- ncol(A)
	
	K <- nrow(cov)
	
	
	
	
	
	
	if (is.null(xprev)) {
		
		nxprev <- 0 
		
		
	} else if (!is.null(xprev)) {
		
		if (is.null(nrow(xprev)) | is.null(ncol(xprev))) {
					
			print("Warning in generate.YuleWalkerCoefficientBlockmatrices: xprev has a NULL number of rows or columns!!!")
			xprev <- NULL		
		}
		
		if (ncol(xprev)!=K) {
			
			print("Warning in generate.YuleWalkerCoefficientBlockmatrices: xprev has a number of columns different form K and it is ingnored!!!")
			xprev <- NULL
		}
		
		if (nrow(xprev)!=p) {
			
			print("Warning in generate.YuleWalkerCoefficientBlockmatrices: xprev has a number of rows different form autoregression order p  and it is ingnored!!!")
			xprev <- NULL
		}
		
		nxprev <- nrow(xprev)
		if (!is.null(nxprev)) nxprev <- 0
		
	}	
	
	out <- array(-9999,c(n+nxprev,K))
	out[,] <- NA
	
	
	
	
	
	
	
	out[1:n+nxprev,] <- as.matrix(generate(x=NULL,FUN=FUN,n=n,K=K,names=names,cov=cov,...))
	
	if (nxprev>0) {
		out[1:nxprev,] <- xprev
	}
	
	#str(out)
	for (i in 1:n+nxprev) {
		
		for (j in 1:p) {
		#	print(i)
		#	print(j)
		#if (i>j) print((out[i-j,]))
		if (i>j) out[i,] <- A[1,j] %*% as.matrix(out[i-j,])+out[i,]
		
		}
	}	
	#

	out <- as.data.frame(out)
	if (!is.null(names)) names(out) <- names 
	
	### INVERSE FUN TRANSFORMATION 
	
	
	if (precipitation.indicator & !is.null(x[[names_x[3]]])) {
		
		pfun <- as.character(substitute(FUN))
		pfun <- get(paste("p",substring(pfun,first=2),sep=""))
		
		prob <- x[[names_x[3]]]$p0_v1
		
		if (ncol(out)!=length(prob)) {
			
			warning("precipitation.indicator option does not work:ncol(out)!=length(prob), Values are NOT transformed!!! ")
			return(out)
		}
		for (c in 1:ncol(out)) {
			
			
			val <- as.vector(pfun(out[,c],...))
			val[which(val<prob[c])] <- 0
			val[which(val>=prob[c])] <- (val[which(val>=prob[c])]-prob[c])/(1-prob[c])
			
			
			out[,c] <- val
			
		}
		
	}
	
	
	
	return(out)
	
	
	
}


NULL
#' @name generate
#' @rdname generate
#' @method generate YuleWalkerCoefficientBlockmatricesPerEachMonth
#' @S3method generate YuleWalkerCoefficientBlockmatricesPerEachMonth
#' @aliases generate generate.YuleWalkerCoefficientBlockmatricesPerEachMonth
#' @export
#' 
#' 

generate.YuleWalkerCoefficientBlockmatricesPerEachMonth <- function(x,FUN=rnorm,year_min=1961,year_max=1990,names=NULL,xprev=NULL,names_x=c("A","Sigma_u","CCGammaInfo"),nearPD=TRUE,precipitation.indicator=FALSE,...) {
	
	### get p and K
	cov <- x[[1]][[names_x[2]]][1,1]
	
###	if (nearPD) cov <- nearPD(cov)$mat 
	
	A <- x[[1]][[names_x[1]]] 
	
	p <- ncol(A)
	
	K <- nrow(cov)
	
	months <- 1:12
	
	start <- paste(year_min,1,1,sep="-")
	end <- paste(year_max,12,31,sep="-")
	days <- as.character(seq(as.Date(start),as.Date(end),by="days"))
	
	out <- as.data.frame(array(NA,c(length(days),K)))
	if (!is.null(names)) names(out) <- names 
	
	out <- adddate(out,origin=start)
	
	pfun <- as.character(substitute(FUN))
	pfun <- get(paste("p",substring(pfun,first=2),sep=""))
	
	
	for (year in year_min:year_max) {
		
		for (m in months) {
			
			index <- which(out$year==year & out$month==m)
			 
			temp <- out[index,] 
			ignore.date <- !(names(temp) %in% c("day","month","year"))
		
			temp[,ignore.date] <- generate(x=x[[m]],FUN=FUN,n=length(index),xprev=xprev,names=names,names_x=names_x,nearPD=nearPD,precipitation.indicator=FALSE,...)  
			
			
			xprev <- temp[nrow(temp)-0:(p-1),ignore.date]
		
			if (precipitation.indicator & !is.null(x[[m]][[names_x[3]]])) {
				
				prob <- x[[m]][[names_x[3]]]$p0_v1
				
				for (c in 1:ncol(temp[,ignore.date])) { 
				
					val <- as.vector(pfun(temp[,ignore.date][,c],...)) ## DA MODIFICARE 
				
					val[val<prob[c]] <- 0
					val[val>=prob[c]] <- (val[val>=prob[c]]-prob[c])/(1-prob[c])
					
					
					temp[,ignore.date][,c] <- val
					
				
				}
				
			}
			## put temp into out 
			out[index,] <- temp 
			
		}
		
		
		
	}
	
###	out <- NULL
	

	out <- out[,!(names(out) %in% c("day","month","year"))]
	return(out)
	
}

NULL 
#' generate
#' 
#' @name generate
#' @rdname generate
#' @method generate CCGammaObject
#' @S3method generate CCGammaObject
#' @aliases generate generate.CCGammaObject
#' @export
#' 
#' 


generate.CCGammaObject <- function(x,n=100,names=NULL,xprev=NULL,precipitation.indicator=TRUE,...) {
	
	
	cov <- x$nooccurence_gcorrelation
	pd <- x$p0_v1  
	mt <- x$TransintionMatrixMCFirstOrder
	K <- nrow(cov)
	if (is.null(names(x))) names <- paste("V",1:nrow(cov))
	out <- generate(x=NULL,FUN=rnorm,n=n,K=K,names=names,cov=cov,...)
	
	if (!precipitation.indicator) {
		
		return(out)
	}
	
    out[,] <- pnorm(as.matrix(out))
	## Collect in a unique data.frame the conditional probability of no precipitation occurence 
	npd <- as.data.frame(array(NA,c(nrow(mt[[1]]),length(mt))))

	names(npd) <- names(out)
	for (l in 1:ncol(npd)) {
		
		npd[,l] <- as.vector(mt[[l]][,"dry"])
		
	}
	
	### r=1 first row 
	
	if (is.null(xprev)) { 
	
		pd <- x$p0_v1
	} else {
		
		precipitation <- xprev>0
		pd <- npd[1,] ## no precipitation conditioned to no precipitation
		pd[precipitation] <- npd[2,][precipitation]  ## no precipitation conditioned to precipitation
		
	}
	## Firest row 
	r=1

	
	
	# ec 20130324
	val <- out[r,] 
	val <- (val-pd)/(1-pd)
	val[val<0] <- 0

	out[r,] <- val
	## end 20130324
	for (r in 2:nrow(out)) { 
	
		## estimate pd
	
	###	print(out[r-1,])
		precipitation <- out[r-1,]>0
	####	print(precipitation)
		pd <- npd[1,] ## no precipitation conditioned to no precipitation
		pd[precipitation] <- npd[2,][precipitation]  ## no precipitation conditioned to precipitation
		

# ec 20130324
		val <- out[r,] 

		val <- (val-pd)/(1-pd)
		val[val<0] <- 0
		out[r,] <- val
##		print(val)
## end 20130324
## commented by ec 20130324
	##	preci <- as.numeric(out[r,]>pd)
	##	val[val!=0] <- (val[val!=0]-pd)/(1-pd)
	##	out[r,] <- val
	
## commented by ec 20130324		

	}
	
	

	return(out)
	
	
	
	
}


NULL 
#' @name generate
#' @rdname generate
#' @method generate CCGammaObjectListPerEachMonth
#' @S3method generate CCGammaObjectListPerEachMonth
#' @aliases generate generate.CCGammaObjectListPerEachMonth
#' @export
#' 
#' 



generate.CCGammaObjectListPerEachMonth <- function(x,year_min=1961,year_max=1990,names=NULL,xprev=NULL,precipitation.indicator=TRUE,...) {

	
	CCGamma <- x[[1]]

	if (class(CCGamma)=="CCGammaObjectList") { 
	
		CCGamma <- CCGamma[[1]]
	}
	
	K <- nrow(CCGamma$nooccurence_gcorrelation)
	
	
	months <- 1:12
	
	start <- paste(year_min,1,1,sep="-")
	end <- paste(year_max,12,31,sep="-")
	days <- as.character(seq(as.Date(start),as.Date(end),by="days"))
	
	out <- as.data.frame(array(NA,c(length(days),K)))
	
	if (!is.null(names)) names(out) <- names 
	
	out <- adddate(out,origin=start)
	
	
	for (year in year_min:year_max) { 
		for (m in months) {
			
			CCGamma <- x[[m]]
			if (class(CCGamma)=="CCGammaObjectList") { 
				
				CCGamma <- CCGamma[[1]]
			}
			
			
			index <- which(out$year==year & out$month==m)
			
			temp <- out[index,] 
			ignore.date <- !(names(temp) %in% c("day","month","year"))
			
			
			temp[,ignore.date] <- generate(x=CCGamma,n=nrow(temp),names=names,xprev=xprev,precipitation.indicator=precipitation.indicator,...) 
			
			xprev <- temp[nrow(temp),ignore.date]
			
			out[index,] <- temp 
			
			
		}
		
		
		
	}
	
	
	out <- out[,!(names(out) %in% c("day","month","year"))]
	return(out)
	
}









