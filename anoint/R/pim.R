setClass("pim",
	representation(
		formula = "anoint",
		coef = "list",
		exact = "logical",
		interval = "numeric",
		LRT = "numeric",
		boot.pim = "matrix",
		vcov = "matrix"
	)
)

pim <- function(object,exact=TRUE,interval=c(-3,3),n.boot=NULL,...){
	
	if(exact){
		fit <- pim.exact(object@formula, object@data, interval=interval,...)
	}
	else{
		fit <- pim.approx(object@formula, object@data)
	}
	
	
	pim <- new("pim",
		formula = object,
		coef = list(alpha=fit$alpha,beta.control=fit$beta,theta=fit$theta),
		exact = exact,
		LRT = fit$LRT,
		interval = interval,
		boot.pim = matrix(0),
		vcov = matrix(0)
	)
	
	if(!is.null(n.boot)){		
		pim@boot.pim <- boot.pim(pim,n.boot)
		pim@vcov <- vcov.pim(coef(pim),pim@boot.pim)
		}
	if(is.null(n.boot)&exact){
		pim@vcov <- fit$vcov
	}

pim
}


pim.print <- function(x,...){
				
				cat("\nBaseline/treatment effects:\n\n")
				print.mat <- cbind(x@coef[[1]])
				
				if(nrow(print.mat)==2)
					row.names(print.mat) <- c("Control","Treatment")
				else
					row.names(print.mat) <- "Treatment"
				colnames(print.mat) <- "Estimate"
				
				print(print.mat)
				
				print.mat <- x@coef[[2]]
				colnames(print.mat) <- "Estimate"
				cat("\nPrognostic effects:\n\n")
				print(print.mat)

				cat("\nResponsiveness parameter:\n\n")
						
				p <- pchisq(x@LRT,df=1,lower.tail=FALSE)
				print.mat <- matrix(c(x@coef[[3]],x@LRT,p),ncol=3)
				row.names(print.mat) <- "theta"
				colnames(print.mat) <- c("Estimate","LRT","p-value")
				
				print(print.mat)
			}

setMethod("print","pim",pim.print)

setMethod("show","pim",function(object) pim.print(object))

setMethod("coef","pim",
	function(object,...){
		mat <- unlist(object@coef)
		
		baseline <- c("Treatment")
		if(length(object@coef[[1]])==2) baseline <- c("Control",baseline)
		
		names(mat) <- c(baseline,row.names(object@coef[[2]]),"theta")
		
		mat
	}
)

setMethod("confint","pim",
	function(object, parm, level=0.95,...){
		
		if(!object@exact&all(dim(object@boot.pim)==c(1,1))){
			stop("No bootstrap resamples were specified in pim call.")
		}
		else if(object@exact&all(dim(object@boot.pim)==c(1,1))){
	
		alpha <- (1-level)/2
		cov.mat <- coef(object)
		cov.mat <- cov.mat[names(cov.mat)!="theta"]
		cov.mat <- c(cov.mat[names(cov.mat)!="Treatment"],cov.mat[names(cov.mat)=="Treatment"])
		cov.p <- length(cov.mat)
		df <- nrow(object@formula@data)-cov.p
		se <- sqrt(diag(vcov(object)))
		
		z <- qt(alpha,df=df)
		lower <- cov.mat-z*se
		upper <- cov.mat+z*se
		
			if(missing(parm)){
				index <- 1:length(lower)
			}
			else{
				index <- match(names(cov.mat),parm)
				index <- index[!is.na(index)]
			}
		
		ci.mat <- cbind(lower[index],upper[index])
		row.names(ci.mat) <- names(cov.mat)[index]
		colnames(ci.mat) <- paste(c(alpha*100,(1-alpha)*100),"%")
		}
		else{		
	
		alpha <- (1-level)/2
		
		lower <- apply(object@boot.pim,1,function(x)quantile(x,alpha))
		upper <- apply(object@boot.pim,1,function(x)quantile(x,1-alpha))
		
			if(missing(parm)){
				index <- 1:length(lower)
			}
			else{
				index <- match(row.names(object@boot.pim),parm)
				index <- index[!is.na(index)]
			}
		
		ci.mat <- cbind(lower[index],upper[index])
		row.names(ci.mat) <- row.names(object@boot.pim)[index]
		colnames(ci.mat) <- paste(c(alpha*100,(1-alpha)*100),"%")
	}
	
	ci.mat
})

setMethod("vcov","pim",function(object,...) object@vcov)

boot.pim <- function(object,n=200){
	
	#BOOTSTRAP RESAMPLES
	resamples <- lapply(1:n,function(void) sample(nrow(object@formula@data),replace=TRUE))

	obs.coef <- coef(object)
	
	fitter <- function(index){
		
		boot.object <- anoint(object@formula@formula@formula,
											object@formula@data[index,],
											family=object@formula@formula@family)
		
		#USE APPROXIMATE METHOD FOR SPEED
		est.coef <- coef(pim(boot.object,exact=FALSE))
	}
	
	fits <- sapply(resamples, fitter)
	
fits
}

vcov.pim <- function(obs.coef,boot.pim){
	
	O <- matrix(obs.coef,nrow(boot.pim),ncol(boot.pim))
	Sigma <- apply(O-boot.pim,2,function(x){outer(x,x)})
	Sigma <- matrix(rowSums(Sigma)/ncol(O),nrow(O),nrow(O))		
	row.names(Sigma) <- row.names(boot.pim)
	colnames(Sigma) <- row.names(boot.pim)
	
Sigma
}
