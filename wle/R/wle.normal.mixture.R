#############################################################
#                                                           #
#	wle.normal.mixture function                             #
#	Author: Claudio Agostinelli                             #
#	E-mail: claudio@unive.it                                #
#	Date: November, 19, 2003                                #
#	Version: 0.2-1                                          #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

wle.normal.mixture <- function(x, m, boot=5, group, num.sol=1, raf="HD", smooth=0.003, tol=10^(-15), equal=10^(-2), max.iter=1000, all.comp=TRUE, min.size=0.02, min.weights=0.3, boot.start=10, group.start=3, tol.start=10^(-6), equal.start=10^(-3), smooth.start=0.003, max.iter.start=500, max.iter.boot=25, verbose=FALSE) {

raf <- switch(raf,
	HD = 1,
	NED = 2,
	SCHI2 = 3,
	-1)

if (raf==-1) stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

if (missing(m)) {
    m <-2
    if (verbose) cat("wle.normal.mixture: number of component set to default value: 2 \n")
}

if (m<2) {
    m <-2
    if (verbose) cat("wle.normal.mixture: number of component set to default value: 2\n")
}

if (missing(group)) {
    group <- 0
}

x <- as.vector(x)
size <- length(x)
result <- list()

if (size<2) {
    stop("Number of observation must be at least equal to 2")
}

if (group<2) {
    group <- max(2,round(size/4))
    if (verbose) cat("wle.normal.mixture: dimension of the subsample set to default value: ",group,"\n")
}

if (group.start<2 | group.start>=group) {
    group.start <- max(2,round(group/4))
    if (verbose) cat("wle.normal.mixture: dimension of the subsample of the starting values set to default value: ",group.start,"\n")
}


maxboot <- sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))

if (boot<1 | log(boot) > maxboot) {
    stop("Bootstrap replication not in the range")
}

maxboots <- sum(log(1:group))-(sum(log(1:group.start))+sum(log(1:(group-group.start))))

if (boot.start<1 | log(boot.start) > maxboots) {
    stop("Bootstrap replication for the starting values not in the range")
}


if (!(num.sol>=1)) {
    if (verbose) cat("wle.normal.mixture: number of solution to report set to 1 \n")
    num.sol <- 1
}

if (max.iter<1) {
    if (verbose) cat("wle.normal.mixture: max number of iteration set to 1000 \n")
    max.iter <- 1000
}

if (max.iter.boot<1) {
    if (verbose) cat("wle.normal.mixture: max number of iteration in the boot step set to 20 \n")
    max.iter.boot <- 20
}

if (max.iter.start<1) {
    if (verbose) cat("wle.normal.mixture: max number of iteration in the starting process set to 500 \n")
    max.iter.start <- 500
}

if (!is.logical(all.comp)) {
    if (verbose) cat("wle.normal.mixture: all.comp must be logical, set to the default value: all.comp \n")
    all.comp <- TRUE 
}

all.comp <- as.numeric(all.comp)

if (smooth<10^(-5)) {
    if (verbose) cat("wle.normal.mixture: the smooth parameter seems too small \n")
}

if (smooth.start<10^(-5)) {
    if (verbose) cat("wle.normal.mixture: the smooth.start parameter seems too small \n")
}

if (tol<=0) {
    if (verbose) cat("wle.normal.mixture: the accuracy must be positive, using default value: 10^(-6) \n")
tol <- 10^(-6)
}

if (equal<=tol) {
    if (verbose) cat("wle.normal.mixture: the equal parameter can not be less than or equal to tol, using default value: tol+10^(-3)\n")
equal <- tol+10^(-3)
}
if (tol.start<=0) {
    if (verbose) cat("wle.normal.mixture: the accuracy in the starting process must be positive, using default value: 10^(-6) \n")
tol.start <- 10^(-6)
}

if (equal.start<=tol.start) {
    if (verbose) cat("wle.normal.mixture: the equal parameter can not be less than or equal to tol.start, using default value: tol.start+10^(-3)\n")
equal.start <- tol.start+10^(-3)
}

loc.boot <- matrix(0,nrow=boot,ncol=m)
var.boot <- matrix(0,nrow=boot,ncol=m)
prop.boot <- matrix(0,nrow=boot,ncol=m)

iboot <- 0
mboot <- 0
while (iboot < boot & mboot <= max.iter.boot) {

    mboot <- mboot + 1
    mmboot <- 0    

    samp <- sample(x,group,replace=FALSE)
    temp.location <- vector(length=0,mode="numeric")
    temp.scale <- vector(length=0,mode="numeric")
    temp.prop <- vector(length=0,mode="numeric")
    temp.weights <- vector(length=0,mode="numeric")
    tot.size <- 1

    while (tot.size > min.size & mmboot <= max.iter.boot) {

           nsamp <- length(samp)
           group.start <- min(group.start,nsamp)
           maxboots <- sum(log(1:nsamp))-(sum(log(1:group.start))+sum(log(1:(nsamp-group.start))))
           if (log(boot.start)>maxboots) boot.start <- max(1,maxboots)

           temp <- wle.normal(samp,num.sol=(m+1),
                              group=group.start,
                              boot=boot.start,
                              smooth=smooth.start,
                              tol=tol.start,equal=equal.start,
                              raf=raf,
                              max.iter=max.iter.start,
                              verbose=verbose)

#       if (verbose) print(temp)

       if (temp$tot.sol>0) {
           if (temp$tot.sol>1) {
               t.location <- temp$location[!(temp$tot.weights==max(temp$tot.weights))]
               t.scale <- temp$scale[!(temp$tot.weights==max(temp$tot.weights))]
               t.weights <- temp$tot.weights[!(temp$tot.weights==max(temp$tot.weights))]
               m.weights <- matrix(temp$weights[!(temp$tot.weights==max(temp$tot.weights)),],nrow=max(1,(temp$tot.sol-1)),byrow=FALSE)
           } else {
               t.location <- temp$location
               t.scale <- temp$scale
               t.weights <- temp$tot.weights
               m.weights <- matrix(temp$weights,nrow=temp$tot.sol)
           }

           temp.location <- c(temp.location,t.location)
           temp.scale <- c(temp.scale,t.scale)
           temp.weights <- c(temp.weights,t.weights)

           samp <- samp[apply(m.weights,2,max)<min.weights]

           tot.size <- length(samp)/group       
     
           if (length(samp)<3 | length(temp.location)>=m) {
               tot.size <- 0
           }

       } else {
           mmboot <- mmboot + 1
       }
       
    }
# end of while (tot.size > min.size)
    mmm <- length(temp.location)
    if (mmm>=m) {
        iboot <- iboot + 1
        if (verbose) cat("Found ",iboot," starting points \n")
        loc.boot[iboot,] <- (temp.location[order(temp.weights)])[(mmm-m+1):mmm]
        var.boot[iboot,] <- ((temp.scale[order(temp.weights)])[(mmm-m+1):mmm])^2
        temp.prop <- (temp.weights[order(temp.weights)])[(mmm-m+1):mmm]
        prop.boot[iboot,] <- temp.prop/sum(temp.prop)
    }

}

if (iboot>0) {


if (verbose) {
    cat("Starting points: \n")
    cat("Location: \n")
    print.default(format(loc.boot, digits=3),
		  print.gap = 2, quote = FALSE)
    cat("Scale: \n")
    print.default(format(sqrt(var.boot), digits=3),
		  print.gap = 2, quote = FALSE)
    cat("Proportion: \n")
    print.default(format(prop.boot, digits=3),
		  print.gap = 2, quote = FALSE)
}

  z <- .Fortran("wlenmix",
	as.double(x), 
	as.integer(size),
        as.integer(m),
        as.matrix(loc.boot),
        as.matrix(var.boot),
        as.matrix(prop.boot),
	as.integer(boot),
	as.integer(num.sol),
	as.integer(raf),
	as.double(smooth),
	as.double(tol),
	as.double(equal),
	as.integer(max.iter),
        as.integer(all.comp),
        as.double(min.weights),
	mean=mat.or.vec(num.sol,m),
	var=mat.or.vec(num.sol,m),
        pi=mat.or.vec(num.sol,m),
	totweight=double(num.sol),
	weight=mat.or.vec(num.sol,size),
	density=mat.or.vec(num.sol,size),
	model=mat.or.vec(num.sol,size),
	delta=mat.or.vec(num.sol,size),
	same=integer(num.sol),
	nsol=integer(1),
	nconv=integer(1),
	PACKAGE = "wle")

keep.col <- 1:z$nsol

if (z$nsol>0) {

    if (num.sol==1) {
        if (any(is.nan(z$mean))) z$nsol <- 0
    } else {
        del.nan <- apply(is.nan(z$mean),1,any)
        keep <- as.numeric(!del.nan)
        keep.col <- c(keep.col,rep(0,num.sol))[1:num.sol]
        keep.col <- keep*keep.col
        keep.col <- keep.col[keep.col!=0]
        z$nsol <- length(keep.col)
    }
}
          
if (z$nsol>0) {

    if (num.sol>1) {
        result$location <- z$mean[keep.col,]
        result$scale <- sqrt(z$var[keep.col,])
        result$pi <- z$pi[keep.col,]
    } else {
        result$location <- z$mean
        result$scale <- sqrt(z$var)
        result$pi <- z$pi
    }
    result$tot.weights <- z$totweight[keep.col]/size
    result$weights <- z$weight[keep.col,]
    result$f.density <- z$density[keep.col,]
    result$m.density <- z$model[keep.col,]
    result$delta <- z$delta[keep.col,]
    result$freq <- z$same[keep.col]
    result$tot.sol <- z$nsol
    result$not.conv <- z$nconv
} else {
    if (verbose) cat("wle.normal.mixture: No solutions are fuond, checks the parameters\n")

    result$location <- rep(NA,m)
    result$scale <- rep(NA,m)
    result$pi <- rep(NA,m)
    result$tot.weights <- NA
    result$weights <- rep(NA,size)
    result$f.density <- rep(NA,size)
    result$m.density <- rep(NA,size)
    result$delta <- rep(NA,size)
    result$freq <- NA
    result$tot.sol <- 0
    result$not.conv <- boot
}

} else {
    if (verbose) cat("wle.normal.mixture: Not able to find a good starting values\n")

    result$location <- rep(NA,m)
    result$scale <- rep(NA,m)
    result$pi <- rep(NA,m)
    result$tot.weights <- NA
    result$weights <- rep(NA,size)
    result$f.density <- rep(NA,size)
    result$m.density <- rep(NA,size)
    result$delta <- rep(NA,size)
    result$freq <- NA
    result$tot.sol <- 0
    result$not.conv <- boot

}

result$call <- match.call()
class(result) <- "wle.normal.mixture"
return(result)
}

#############################################################
#                                                           #
#	print.wle.normal.mixture function                   #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: October, 3, 2001                              #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

print.wle.normal.mixture <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("Location:\n")
    print.default(format(x$location, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Scale:\n")
    print.default(format(x$scale, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Proportion:\n")
    print.default(format(x$pi, digits=digits),
		  print.gap = 2, quote = FALSE)


    cat("\n")
    cat("\nNumber of solutions ",x$tot.sol,"\n")
    cat("\n")
    invisible(x)
}






     




