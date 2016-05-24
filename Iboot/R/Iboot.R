#Last modified 2013/02/15  by Nicola Lunardon.

Iboot <- function(gradient, B=500, M=500, kB=0.01, alpha=c(0.1, 0.05), control.optim=list(), seed)
{
	###some checks
		if(missing(gradient)) stop("'gradient' is required \n")
		if(B<0 | M<0) stop("'B' and 'M' must be greater than zero \n")
		if(kB<0) warning("To run the algorithm 'kB' must be greater than zero \n")
		if(any(alpha<=0) | any(alpha>=1)) stop("'alpha' must be in the interval (0,1) \n")
	###
	
	gradient <- data.matrix(gradient)
	
	n <- NROW(gradient)
	p <- NCOL(gradient)

	alpha_max <- max(alpha)
	R <- round(alpha_max*(B+1))
		if(R==0) stop("need to increase B\n")

	#control parameters for L-BFGS-B
	con <- list(trace = 0, maxit = 2e4, pgtol=sqrt(.Machine$double.eps), REPORT = 10, lmm = 5, factr = 1e+07)
	con[names(control.optim)] <- control.optim

	if(con$trace<0) stop("trace must be greater than 0 \n")

	if(!missing(seed)) set.seed(seed)

	result <- .C("bstr_iboot_stop",
															oss=double(1), 
															boot=double(B), 
															prep=double(R), 
															gradient=as.double(gradient), 
															n=as.integer(n), 
															p=as.integer(p), 
															B=as.integer(B), 
															M=as.integer(M), 
															R=as.integer(R), 
															kB=as.integer(max(round(kB*B),1)),
															fails.outer=double(1),
															count_good=integer(1),
															ov.fail=as.integer(0), 
															pgtol=as.double(con$pgtol),
															factr=as.double(con$factr),
															lmm=as.integer(con$lmm),
															trace=as.integer(con$trace),
															maxit=as.integer(con$maxit),
															report=as.integer(con$REPORT), 
															PACKAGE="Iboot")
	
	if(result$ov.fail==0)
	{
		##output for double bootstrap
		if(M>0)
		{
			index <- B-round(B*(1-alpha))
			prep <- result$prep[index]/M
			obj.iboot <- list(Call=match.call(), oss=result$oss/n, boot=result$boot/n, map=prep, boot.quant=quantile(result$boot/n, 1-alpha), recalib.quant=quantile(result$boot/n, prep), fails.outer=result$fails.outer, failure=result$ov.fail)
		}
		else ##output for simple bootstrap
		{
			obj.iboot <- list(Call=match.call(), oss=result$oss/n, boot=result$boot/n, map=NULL, boot.quant=quantile(result$boot/n, 1-alpha), recalib.quant=NULL, fails.outer=result$fails.outer, failure=result$ov.fail)
		}
	}	
	else
	{
		if(result$ov.fail==2)
		{
			warning(paste("Algorithm has stopped in the outer level. Only ", result$count_good, " resamplings out of ", B, " have been performed\n", sep=""))
			obj.iboot <- list(Call=match.call(), oss=result$oss/n, boot=result$boot[1:result$count_good], map=NULL, boot.quant=quantile(result$boot[1:result$count_good]/n, 1-alpha), recalib.quant=NULL, fails.outer=result$fails.outer, failure=result$ov.fail)
		}
		else
		{
			warning("Algorithm has stopped at the beginning. No resamplings have been performed\n")
			obj.iboot <- list(Call=match.call(), oss=result$oss/n, boot=NULL, map=NULL, boot.quant=NULL, recalib.quant=NULL, fails.outer=NULL, failure=result$ov.fail)
		}
	}
	class(obj.iboot) <- "Iboot"
	obj.iboot
}


##internal: print for Iboot
print.Iboot <- function(x, ...) 
{
	dgts.def <- getOption("digits")
	options(digits=4)
	
	if(x$failure==0) 
	{
		cat("\nObserved value:", format(x$oss, digits=getOption("digits")),"\n")
		
			if(!is.null(x$boot.quant) & !is.null(x$recalib.quant))
			{
				cat("\nBootstrap quantile(s):\n")
				print(x$boot.quant)
				cat("\nRe-calibrated bootstrap quantile(s):\n")
				print(x$recalib.quant)
				cat("\nAlgorithm ended succesfully.\nActual proportion of convex hull condition failures:", format(x$fails.outer, digits=getOption("digits")),"\n\n")
			}
			else
			{
				cat("\nBootstrap quantile(s):\n")
				print(x$boot.quant)
				cat("\nAlgorithm ended succesfully.\nActual proportion of convex hull condition failures:", format(x$fails.outer, digits=getOption("digits")),"\n\n")
			}

	}
	else
	{
		if(x$failure==1)
		{
			cat("\nObserved value:", format(x$oss, digits=getOption("digits")),"\n")
			cat("\nAlgorithm has stopped at the beginning\n\n")
		}
		else
		{
				cat("\nObserved value:", format(x$oss, digits=getOption("digits")),"\n")
				cat("\nBootstrap quantile(s):\n")
				print(x$boot.quant)
				cat("\nAlgorithm has stopped in the outer level. Only ", length(x$boot), " resamplings out of ", x$Call$B, " have been performed\n", sep="")
		}
	}
	on.exit(options(digits=dgts.def))	
}

##internal: summary for Iboot
summary.Iboot <- function(object, ...) 
{
	dgts.def <- getOption("digits")
	options(digits=4)
	
	if(object$failure==0) 
	{
		cat("\n")
      cat("Call: \n")
		print(object$Call)

		cat("\nBootstrap distribution:\n")
		print(summary(object$boot))
		cat("\nObserved value:", format(object$oss, digits=getOption("digits")),"\n")

		
			if(!is.null(object$boot.quant) & !is.null(object$recalib.quant))
			{
				cat("\nBootstrap quantile(s):\n")
				print(object$boot.quant)
				cat("\nRe-calibrated bootstrap quantile(s):\n")
				print(object$recalib.quant)
				cat("\nAlgorithm ended succesfully.\nActual proportion of convex hull condition failures:", format(object$fails.outer, digits=getOption("digits")),"\n\n")
			}
			else
			{
				cat("\nBootstrap quantile(s):\n")
				print(object$boot.quant)
				cat("\nAlgorithm ended succesfully.\nActual proportion of convex hull condition failures:", format(object$fails.outer, digits=getOption("digits")),"\n\n")
			}

	}
	else
	{
		if(object$failure==1)
		{
			cat("\n")
			cat("Call: \n")
			print(object$Call)

			cat("\nObserved value:", format(object$oss, digits=getOption("digits")),"\n")
			cat("\nAlgorithm has stopped at the beginning\n\n")
		}
		else
		{
				cat("\n")
				cat("Call: \n")
				print(object$Call)

				cat("\nBootstrap distribution:\n")
				print(summary(object$boot))

				cat("\nObserved value:", format(object$oss, digits=getOption("digits")),"\n")
				cat("\nBootstrap quantile(s):\n")
				print(object$boot.quant)
				cat("\nAlgorithm has stopped in the outer level. Only ", length(object$boot), " resamplings out of ", object$Call$B, " have been performed\n", sep="")
		}
	}
	on.exit(options(digits=dgts.def))	
}


