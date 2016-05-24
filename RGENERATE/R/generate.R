# TODO: Add comment
# 
# Author: ecor
###############################################################################

NULL

#' generate
#' 
#' 
#' 
#' @rdname generate
#' @export

generate <- function (x=NULL,...)  {
	
	
	return(UseMethod("generate",x))
	
}


NULL
#' 
#' 

#' 
#' @rdname generate
#' @method generate default
#' @S3method generate default
#' @aliases generate 
#' @export



generate.default <- function (x,FUN=rnorm,n=100,K=3,names=NULL,cov=NULL,gap.filling=NULL,...)  {
	
	if (!is.null(names)) K <- length(names)
	if (!is.null(cov)) K <- ncol(cov)
	
	out <- FUN(n*K,...)
	
	out <- array(out,c(n,K))
	out <- as.data.frame(out)
	if (!is.null(names)) {
		if (length(names)==K) names(out) <- names
	}
# TRASH .... str	
#	if (is.null(B)) B <- t(chol(summary(var)$covres))
	
#	if (is.null(xprev)) xprev <- rnorm(ncol(B))
	
	
	if (!is.null(cov)) {
		B <- t(chol(cov))
		out[,] <- ((as.matrix(out)) %*% t(B))
	}
	
	if (!is.null(gap.filling)) {
		
		
		out[!is.na(gap.filling)] <- gap.filling[!is.na(gap.filling)]
		
	}
	
	return(out)
	
}


NULL
#' 
#' It generates a multivarite random series according to the model \code{x}
#' 
#' 
#' 
#' 
#' @param x null object or the model used for random generation , e.g. a VAR model as a \code{\link{varest-class}} or \code{\link{varest2-class}} object. Default is \code{NULL}.
#' @param FUN random function of the probability distribution used for noise random generation. Default is \code{\link{rnorm}}. See \url{http://cran.r-project.org/web/views/Distributions.html} 
#' @param n number of generations requested 
#' @param names null object or string vectors or names of the variables to be generated simultaneously. Default is \code{NULL}.
#' @param K number of the variables to be generated simultaneously, i.e. the K parameters of a VAR. It is automatically detected by \code{x}, \code{names} or \code{cov}, if one of these is not \code{NULL}. 
#' @param cov null object or covariance matrix of the random variables to be generated simultaneously. Default is \code{NULL}, not used in case this information can be detected from \code{x}.
#' @param noise null object or a generic external noise for \code{x} model residuals, e.g. standard white noise, for random generation with the model \code{x}. Default is \code{NULL}. If \code{NULL} the noise is automatically calculated. 
#' @param exogen null object or amatrix or data frame with exogeneous variables (predictors) id requested by \code{x}. Default is \code{NULL}  
#' @param xprev null object or initial condition of the multivariate random process to be generated. Default is \code{NULL}. 
#' @param gap.filling data frame with time series with gabs (\code{NA} values) to be filled. Default is \code{NULL} and not considered, otherwise the method returns this data frame with \code{NA}  row replaced with generated (e.g auto-regressed) values. 
#' @param GPCA.row.gap.filling.option logical value. Default is \code{TRUE}. In case of \code{\link{GPCAvarest2-class}} objects, If \code{gap.filling} contains both \code{NA} and finite values in the same row, 
#' this row will contains all \code{NA} values after GPCA. In this case all row values are generated through auto-regression. If \code{GPCA.row.gap.filling.option} all insterted non-NA \code{gap.filling} values   are repleced before returning the function value. 
#' Otherwise, in the rows with \code{NA}s all values are re-generated. The option \code{TRUE} is not safe in case the gaps are vary long becouse the genereted values is used for subsequent auto-regrossion.
#' @param type character string used in some method implementations. See \code{\link{inv_GPCA}}. In the matrix implementation, default is \code{"autoregression"}, i.e. the matrix is used as a vector auto-regression coefficient, if it is \code{"covariance"} the method genereted a sample with covariance matrix given by \code{x}.
#' @param ... further arguments for \code{FUN}
#' 
#' @return a matrix or a data frame object
#' 
#' @seealso \code{\link{getVARmodel}}
#' 
#' @title generate
#' @name generate
#' @rdname generate
#' @method generate varest
#' @S3method generate varest
#' @aliases generate generate.varest 
#' 
#' @export
#' 
#' @import RMAWGEN 
#' @examples 
#' 
#' 
#' 
#' 
#' library(RGENERATE)
#' 
#' set.seed(122)
#' NSTEP <- 1000
#' x <- rnorm(NSTEP)
#' y <- x+rnorm(NSTEP)
#' z <- c(rnorm(1),y[-1]+rnorm(NSTEP-1))
#' df <- data.frame(x=x,y=y,z=z)
#' var <- VAR(df,type="none")
#' gg <- generate(var,n=20) 
#'
#' cov <- cov(gg)
#' 
#' ggg <- generate(FUN=rnorm,n=NSTEP,cov=cov)
#' 
#'  
#' library(RMAWGEN)
#' exogen <- as.data.frame(x+5)
#' gpcavar <- getVARmodel(data=df,suffix=NULL,p=3,n_GPCA_iteration=5,
#' n_GPCA_iteration_residuals=5,exogen=exogen)
#' gpcagg <- generate(gpcavar,n=20,exogen=exogen) 
#' 
#' ## Generate an auto-regrassive time-series with a generic matrix 
#' 
#' A <- diag(c(1,-1,1))
#' mgg <- generate(A,n=100)
#' 
#' ### Gap Filling Examples
#' 
#'  dfobs <- df
#'  dfobs[20:30,] <- NA 
#'  n <- nrow(df)
#'  dffill <- generate(gpcavar,n=n,exogen=exogen,gap.filling=dfobs,names=names(dfobs)) 
#'  
#' qqplot(dfobs$y,dffill$y)
#' abline(0,1)
#' 
#' ### Gap filling with matrix 
#' 
#' mgg_n <- mgg
#' mgg_n[20:30,2] <- NA 
#' 
#' mgg_nfill <- generate(A,gap.filling=mgg_n)
#' 
#' print(mgg_n[1:31,])
#' print(mgg_nfill[1:31,])
#'
#' dfobs2 <- df
#' dfobs2[20:30,2] <- NA
#' n <- nrow(df)
#' dffill2 <- generate(gpcavar,n=n,exogen=exogen,gap.filling=dfobs2,names=names(dfobs2)) 
#' 
#' qqplot(dfobs$y,dffill$y)
#' abline(0,1)
#'  
#' ### generation with 'generetion.matrix' 
#' ### and matrix 'x' is a covariance matrix 
#' 
#' covariance <- array(0.5,c(3,3))
#' 
#' diag(covariance) <- 1
#' 
#' set.seed(127)
#' ngns <- 1000
#' gg1 <- generate(FUN=rnorm,n=ngns,cov=covariance)
#' set.seed(127)
#' gg2 <- generate(covariance,type="covariance",n=ngns)
#' 
#' 
#' ## generate with a list of covariance matrix 
#' ndim <- 5
#' dim <- c(ndim,ndim)
#' CS1 <- array(0.3,dim)
#' CS2 <- array(0.5,dim)
#' CS3 <- array(0.7,dim)
#' CS4 <- array(0.1,dim)
#' 
#' diag(CS1) <- 1
#' diag(CS2) <- 1
#' diag(CS3) <- 1
#' diag(CS4) <- 1
#' 
#' list <- list(CS1=CS1,CS2=CS2,CS3=CS3,CS4=CS4)
#' 
#' series <- rep(1:4,times=4,each=100)
#' series <- sprintf("CS%d",series)
#' names_A <- sprintf("A%d",1:ndim)
#' ggs <- generate(list,factor.series=series,FUN=rnorm,type="covariance",names=names_A)
#' 
#' ggs_CS1 <- ggs[series=="CS1",]
#' cov(ggs_CS1)
#' 
#' ggs_CS3 <- ggs[series=="CS3",] 
#' cov(ggs_CS3)

generate.varest <- function (x,FUN=rnorm,n=100,names=NULL,noise=NULL,exogen=NULL,xprev=NULL,gap.filling=NULL,...)  {

	# x is varest objesct
#	out <- require("vars")
#	if (out==FALSE) {
#		
#		print("The function cannot work correctly and returns FALSE")
#		return(out)
#		
#	}
	
	K <-x$K
	p <- x$p
	
	if (x$type!="none") {
	   message("Warning in generate.varest: VAR model type is different from 'none' ")
	   message("Check VAR Model, this method might not work successfully!!")
	}	
	nexogen <- ncol(x$datamat)-K*(p+1)
	
	## CHECK exogen var!!! 
	if ((is.null(exogen)) & (nexogen!=0)) {
		message("Error in generate.varest method: exogen variables (predictors) are needed for VAR multirealization") 
		message("Check VAR and exogen variables!!!")
		print(nexogen)
		stop()
	} else if (!is.null(exogen)) if (nexogen!=ncol(exogen)) {
			message("Error in generate.varest method:  corrected exogen variables (predictors) are needed for VAR multirealization") 
			message("Check VAR and exogen variables!!!")
			stop()
	}

	cvar <- coef(x)
	lcvar <- lapply(X=cvar,FUN=function(x){ x[,"Estimate"]})
	mcvar <- matrix(unlist(lcvar),nrow=length(lcvar),byrow=TRUE)
	rownames(mcvar) <- names(lcvar)
	colnames(mcvar) <- names(lcvar[[1]])
	
	
	
	
	if (!is.null(exogen)) n <- nrow(exogen)
	
	if (is.null(xprev)) xprev <- generate(FUN=FUN,n=p,K=K,names=names,cov=summary(x)$covres,...)
	nxprev <- length(xprev)

	if (is.null(noise)) { 
		noise <- generate(FUN=FUN,n=n,K=K,names=names,cov=summary(x)$covres,...)
	
	
	}
	
	
	if (!is.null(exogen)) {
		if (nrow(noise)!=nrow(exogen)) { 
			warning("Warning in generate.varest method:  exogen realizations (predictors) are not equal to noise realizations") 
			nmin <- min(c(nrow(exogen),nrow(noise)))
			exogen <- exogen[1:nmin,]
			noise <- noise[1:nmin,]
			
		}
	}	
	
	out <- noise
# noise and exogen are trasformed into matrix
	noise <- as.matrix(noise)
	if (!is.null(exogen)) exogen <- as.matrix(exogen)
	
	xprev <- as.vector(t(as.matrix(xprev[nrow(xprev):1,])))
	
	if (!is.null(gap.filling)) {
		
		if (nrow(noise)!=nrow(gap.filling)) {
			warning("Warning in generate.varest method: gap.filling is not equal to noise realizations (i.e. the number of realizations!!), then it is not considered!!!") 
			gap.filling <- NULL
			
			}
		
		
	}
	
	if (is.null(gap.filling)) { 
		
	 	count.gab.filling <-  1:nrow(noise)
		gap.filling.action=FALSE
		
	} else {
		
		count.gab.filling <- lapply(X=1:nrow(gap.filling), FUN=function(x,df) {length(df[x,])!=length(which(!is.na(df[x,])))},df=gap.filling)
		count.gab.filling <- unlist(count.gab.filling)
		count.gab.filling <- which(count.gab.filling)
		gap.filling.action=TRUE
		###print(count.gab.filling)
		
	} 
	for (i in 1:nrow(noise)) {
		
##		print(xprev)
##		print(i)
		
		if (is.null(exogen)) { 
			xprevx <- as.vector(xprev)
		} else { 
			xprevx <- as.vector(xprev)

			xprevx <- (c(as.vector(xprev),as.vector((exogen[i,]))))
		}
		
		if (i %in% count.gab.filling) {
			
			xout <- as.vector(mcvar %*% as.matrix(xprevx)) + as.vector(noise[i,])
			
			if (gap.filling.action==TRUE) {
				
				xout_val <- as.vector(unlist(gap.filling[i,]))
				
				xout[!is.na(xout_val)] <- xout_val[!is.na(xout_val)]
			}
		} else { 
			
			xout <- as.vector(unlist(gap.filling[i,]))
			
	##		xout_var <- as.vector(mcvar %*% as.matrix(xprevx)) + as.vector(noise[i,])
	##		xout[is.na(xout)] <- xout_var[is.na(xout)]
			
		}
		
	
		if (p>1) { 
			xprev <- c(xout,xprev[1:(K*(p-1))])
		} else { 
			xprev[] <- as.vector(xout)
		}
		
		out[i,] <- xout
	
	}	
	
	return(out)
	
}	

NULL

#' 
#' @rdname generate
#' @method generate varest2
#' @S3method generate varest2
#' @aliases generate 
#' @export

generate.varest2 <- function (x,FUN=rnorm,n=100,names=NULL,noise=NULL,exogen=NULL,xprev=NULL,gap.filling=NULL,...) { 

	out <- generate(x=x@VAR,FUN=FUN,n=n,names=names,noise=noise,exogen=exogen,xprev=xprev,gap.filling=gap.filling,...)
	
	return(out)
	
 
}	


NULL

#' 
#' 
#' @param extremes  see \code{\link{inv_GPCA}}
#' 
#' @rdname generate
#' @method generate GPCAvarest2
#' @S3method generate GPCAvarest2
#' @aliases generate 
#' @export

generate.GPCAvarest2 <- function (x,FUN=rnorm,n=100,names=NULL,noise=NULL,exogen=NULL,xprev=NULL,extremes=TRUE,type=3,gap.filling=NULL,GPCA.row.gap.filling.option=TRUE,...) { 
	
	# vedere codice RMAWGEN 
	
	K <-x@VAR$K
	p <- x@VAR$p
	
	if (!is.null(exogen)) n <- nrow(exogen)
	if (!is.null(noise))  n <- min(c(nrow(noise),nrow(exogen)))
	
	if (is.null(noise)) {
		
		cov <- cov(residuals(x))
		noise <- generate(FUN=FUN,n=n,K=K,names=names,cov=cov,...)
		
	}

	varresiduals <- inv_GPCA(x=noise,GPCA_param=x@GPCA_residuals,type=type,extremes=extremes)
	
	
	gap.filling.GPCA <- NULL
    if (is.null(gap.filling)) {
		
		GPCA_data <- x@GPCA_data
		GPCA.gap.filling <- NULL
		GPCA.row.gap.filling.option==FALSE
	} else { 
		GPCA_data <- GPCA(gap.filling, n = length(x@GPCA_data)-1, extremes = extremes)
		GPCA.gap.filling <- GPCA_data$final_results
	}
	
	 
	out <- generate(x=x@VAR,FUN=FUN,n=n,names=names,noise=varresiduals,exogen=exogen,xprev=xprev,gap.filling=GPCA.gap.filling ,...)

		
	
	out <- inv_GPCA(x=out,GPCA_param=GPCA_data,type=type,extremes=extremes)

	
	names(out) <- names
	
	if (GPCA.row.gap.filling.option==TRUE) {
		
		count <- !is.na(gap.filling)
		out[count] <- gap.filling[count]
		
		
	}
	
	return(out)
	
	
}	

NULL
#' 
#' 

#' @rdname generate
#' @method generate matrix
#' @S3method generate matrix
#' @aliases generate 
#' @export


#####################################################################
#####################################################################
#####################################################################

generate.matrix <- function (x,FUN=rnorm,n=100,noise=NULL,xprev=NULL,names=NULL,gap.filling=NULL,type=c("autoregression","covariance"),...) { 
	
	# vedere codice RMAWGEN 
	
	
	
	
	
	if (!is.null(noise))  n <- nrow(noise)
	type <- type[1]
	
	if (type=="covariance") {
		
		K <- nrow(x)
		out <- generate(FUN=FUN,n=n,cov=x,names=names,...)
		return(out)
	}
	
	
	if (is.null(noise)) {
		K <- nrow(x)
		noise <- generate(FUN=FUN,n=n,K=K,...)
		
		
	}
	
	if (is.null(xprev)) {
		
		xprev <- t(noise[1,])
	}
	out <- noise 
	
	if (!is.null(gap.filling)) {
		
		if (nrow(noise)!=nrow(gap.filling)) {
			warning("Warning in generate.matrix method: gap.filling is not equal to noise realizations (i.e. the number of realizations!!), then it is not considered!!!") 
			gap.filling <- NULL
			
		}
		
		
	}
	
	if (is.null(gap.filling)) { 
		
		count.gab.filling <- 1:nrow(noise)
		gap.filling.action=FALSE
		
	} else {
		
		count.gab.filling <- lapply(X=1:nrow(gap.filling), FUN=function(x,df) {length(df[x,])!=length(which(!is.na(df[x,])))},df=gap.filling)
		count.gab.filling <- unlist(count.gab.filling)
		count.gab.filling <- which(count.gab.filling)
		gap.filling.action=TRUE
		
		
	} 
	
	
	
	
	for (i in 1:nrow(noise)) { 
	       
			
			if (i %in% count.gab.filling) {
				xnext <- as.vector(x %*% as.matrix(xprev)) + as.vector(noise[i,])	
				if (gap.filling.action==TRUE) {
					xnext_val <- as.vector(gap.filling[i,])
					xnext[!is.na(xnext_val)] <- xnext_val[!is.na(xnext_val)]
				}
			} else { 
			
				xnext <- as.vector(gap.filling[i,])
		##		xnext_var <- as.vector(x %*% as.matrix(xprev)) + as.vector(noise[i,])	
		##		xnext[is.na(xnext)] <- xnext_var[is.na(xnext)]
				
			}
			out[i,] <- xnext
			xprev <- t(out[i,])
			
	
	}
	
	if (!is.null(names)) names(out) <- names
	
	return(out)
	
	
}	

NULL
#' 
#' 
#' @param  factor.series factor series used by 'factor.series' 
#' 
#' 
#' @rdname generate
#' @method generate list
#' @S3method generate list
#' @aliases generate 
#' @export


generate.list <- function(x,factor.series=names(x),n=NA,...) {
	
	out <- NULL
	it <- as.character(factor.series[1])
	
	rows <- which(as.character(factor.series)==it)
	temp <- generate(x[[it]],n=length(rows),...)
	
	out <- array(NA,c(length(factor.series),ncol(temp)))
	out <- as.data.frame(out)
	names(out) <- names(temp)
	
	out[rows,] <- temp
    ##for (it )
	###factor.series.a <- factor.series[-rows]
	for (id in names(x)[names(x)!=it]) {
		
	
		rows <- which(as.character(factor.series)==id)
		if (length(rows)>0) out[rows,] <- generate(x[[id]],n=length(rows),...)
		
		
	}
	

	return(out)
	
	
	
}

NULL
#' 
#' 
#' @param  origin start date for generation. See \code{\link{adddate}}
#' 
#' 
#' @rdname generate
#' @method generate MonthlyList
#' @S3method generate MonthlyList
#' @aliases generate 
#' @export

generate.MonthlyList <- function(x,origin,n,...) {
	
	monthly.factor  <- adddate(data=as.data.frame(1:n),origin=origin)$month
	
	class(x) <- "list"
	names(x) <- sort(unique(monthly.factor))
	out <- generate(x,factor.series=factor(monthly.factor),n=NA,...)
	
	
}


