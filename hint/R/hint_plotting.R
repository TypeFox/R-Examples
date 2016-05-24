#########################################################
#
# Plotting functions for R package 'hint'.
#
# Author: Alex T. Kalinka (alex.t.kalinka@gmail.com)
#
#########################################################



plot.hint.test <- function(x, ...)
	# S3 method for "plot" generic function.
	# test is a "hint.test" object.
	{
	plotDistr(x, ...)
	return(invisible())
	}



# Check input for plotting function and return necessary variables.
.check.distr.plotting <- function(distr)
	{
	if(class(distr)=="data.frame" || class(distr)=="matrix"){
		if(ncol(distr) == 1){
			warning(paste("no random variable values given: using 0-",(nrow(distr)-1),sep=""), call. = FALSE)
			rv <- 0:(nrow(distr)-1)
			pv <- distr[,1]
			cutoff <- alt <- NA
		}else if(ncol(distr)==2){
			rv <- distr[,1]
			pv <- distr[,2]
			cutoff <- alt <- NA
		}else{
			stop("distr must be a data frame or matrix with 1 or 2 columns\n", call. = FALSE)
			}
	}else if(class(distr)=="hint.test"){
		dd <- NULL
		params <- distr$parameters
		md <- match("d",names(params))
		if(!is.na(md)){ # distance distribution.
			dd <- .hint.dist.distr(params[2],c(params[3],params[4]),params[6],c(params[7],params[8]),params[5],params[9])
			rv <- dd[,1]
			pv <- dd[,2]
			cutoff <- params[1]
			alt <- distr$alternative
		}else{
			dd <- dhint(params[1],c(params[2],params[3]),params[4])
			rv <- dd[,1]
			pv <- dd[,2]
			cutoff <- params[5]
			alt <- distr$alternative
			}
	}else if(class(distr)=="numeric"){
		warning(paste("no random variable values given: using 0-",(length(distr)-1),sep=""), call. = FALSE)
		rv <- 0:(length(distr)-1)
		pv <- distr
		cutoff <- alt <- NA
	}else{
		stop("distr must be a valid distribution or statistical test\n", call.=FALSE)
		}
	ret <- list()
	ret$distr <- data.frame(random.variable=rv, p=pv)
	ret$cutoff <- cutoff
	ret$alt <- alt
	return(ret)
	}



plotDistr <- function(distr, col = "black", test.col = "red", xlim = NULL, ylim = NULL, xlab = "Intersection size (v)", ylab = "Probability", add = FALSE, ...)
	{
	# distr is a distribution as either a data frame or a matrix containing P-values for a discrete distribution:
	#  v            p
	#  0 6.119727e-08
	#  1 1.876133e-06
	#  2 2.626587e-05
	#  3 2.236148e-04
	#  4 1.300497e-03
	#  5 5.502102e-03
	# or a test object (e.g. class "hint.test").
	# or a vector of numbers (but here we must guess the random variable values).
	# Ultimately, we would like to include continuous distributions and their test objects (t-test, chi-squared, etc)
	#  for pedagogical purposes.
	#
	distr <- .check.distr.plotting(distr)
	# Get extremities for plot region.
	if(is.null(xlim)){
		xmin <- 0
		xmax <- max(distr$distr[,1])
	}else{
		xmin <- xlim[1]
		xmax <- xlim[2]
		}
	if(is.null(ylim)){
		pmin <- 0
		pmax <- max(distr$distr[,2]) + 0.02
	}else{
		pmin <- ylim[1]
		pmax <- ylim[2]
		}
	if(!add){
		dev.hold()
		on.exit(dev.flush())
		plot(0, 0, xlim = range(xmin,xmax), ylim = range(pmin,pmax), type="n", xlab=xlab, ylab=ylab, ...)
		}
	dd <- distr$distr
	cutoff <- distr$cutoff
	alt <- distr$alt
	for(j in 1:nrow(dd)){
		if(!is.na(cutoff)){
			if(alt == "greater"){
				if(dd[j,1] >= cutoff){
					tcol <- test.col
				}else{
					tcol <- col
					}
			}else if(alt == "less"){
				if(dd[j,1] <= cutoff){
					tcol <- test.col
				}else{
					tcol <- col
					}
			}else if(alt == "two.sided"){
				expect <- sum(dd[,1]*dd[,2])
				cd <- cutoff-expect
				if(cd >= 0){
					if(dd[j,1] >= cutoff){
						tcol <- test.col
					}else{
						tcol <- col
						}
				}else{
					if(dd[j,1] <= cutoff){
						tcol <- test.col
					}else{
						tcol <- col
						}
					}
				}
		}else{
			tcol <- col
			}
		segments(dd[j,1], 0, dd[j,1], dd[j,2], col = tcol, lwd = 2)
		}
	return(invisible())
	}



# Add segment plot distribution(s) to an already live plot.
add.distr <- function(..., cols = "blue", test.cols = "red")
	{
	# cols are the colours of each distribution.
	# test.cols are the cutoff cols if needed.
	# ... are one or more distributions or statistical tests.
	alld <- list(...)
	if(length(alld)==0){
		stop("one or more distributions must be supplied in the ... argument\n", call. = FALSE)
		}
	j <- k <- 1
	for(i in 1:length(alld)){
		if(i > length(cols)){
			j <- 1
			}
		if(i > length(test.cols)){
			k <- 1
			}
		plotDistr(alld[[i]], col = cols[j], test.col = test.cols[k], add = TRUE)
		j <- j + 1
		k <- k + 1
		}
	return(invisible())
	}







