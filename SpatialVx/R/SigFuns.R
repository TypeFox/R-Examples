EBS <- function(object, model=1, block.length=NULL, alpha.boot=0.05, field.sig=0.05, bootR=1000, ntrials=1000, verbose=FALSE) {

    a <- attributes(object)
    out <- list()
    attributes(out) <- a

    X <- object[[1]]
    Xhat <- object[[2]]
    if(!is.array(Xhat)) {
	if(!is.numeric(model)) {
	    nf <- a$nforecast
	    dn <- a$data.name
	    if(length(dn) == nf + 2) mod.names <- dn[-(1:2)]
            else mod.names <- dn[-1]
            model <- (1:nf)[dn == model]
            if(is.na(model)) stop("datagrabber: invalid model argument.")
	}
	Xhat <- Xhat[[model]]
    } # end of if 'Xhat' is not an array (i.e., if more than one model in object) stmts.

    X <- t(apply(X, 3, c))
    Xhat <- t(apply(Xhat, 3, c))
    loc <- a$loc

    if(!is.null(a$subset)) {
	id <- a$subset
	X <- X[,id]
	Xhat <- Xhat[,id]
	loc <- loc[id,]
    }

    res <- spatbiasFS(X=X, Y=Xhat, loc=a$loc, block.length=block.length,
		    alpha.boot=alpha.boot, field.sig=field.sig, bootR=bootR,
		    ntrials=ntrials, verbose=verbose)

    out$block.boot.results <- res$block.boot.results
    out$sig.results <- res$sig.results
    sig.values <- c(res$field.significance, res$alpha.boot, bootR, ntrials)
    names(sig.values) <- c("field sig.", "alpha boot", "bootstrap replicates", "number of trials")
    attr(out, "arguments") <- sig.values
    attr(out, "which.model") <- model
    class(out) <- "EBS"
    return(out)
} # end of 'EBS' function.

LocSig <- function(Z, numrep=1000, block.length=NULL, bootfun="mean", alpha=0.05, bca=FALSE, ...) {
   if(bootfun=="mean") bootfun <- function(data) return(colMeans(data, na.rm=TRUE))
   else if(is.character(bootfun)) bootfun <- get(bootfun)
   zdim <- dim(Z)
   n <- zdim[1]
   m <- zdim[2]
   out <- data.frame(Lower = numeric(m), Estimate = numeric(m), Upper = numeric(m))
   if(is.null(block.length)) block.length <- floor(sqrt(n))
   if(block.length==1) booted <- boot(Z, bootfun, R=numrep, ...)
   else booted <- tsboot(Z, bootfun, l=block.length, R=numrep, sim="fixed", ...)
   out$Estimate <- booted$t0
   if((block.length==1) & bca) {
	for(i in 1:m) {
	   tmp <- boot.ci( booted, conf=1-alpha, type="bca", index=i)
	   out$Lower[i] <- tmp$bca[,4]
	   out$Upper[i] <- tmp$bca[,5]
	} # end of for 'i' loop.
   } else {
	if(bca) warning("LocSig: You chose to use the BCa method, but block.length != 1.  Using percentile method with circular block bootstrap instead.")
	for(i in 1:m) {
	   tmp <- boot.ci( booted, conf=1-alpha, type="perc", index=i)
	   out$Lower[i] <- tmp$perc[,4]
	   out$Upper[i] <- tmp$perc[,5]
	} # end of for 'i' loop.
   } # end of if else do "BCa" (IID bootstrap only) or percentile confidence limits.
   class(out) <- "LocSig"
   return(out)
} # end of 'LocSig' function.

MCdof <- function(x, ntrials=5000, field.sig=0.05, zfun="rnorm", zfun.args=NULL, which.test=c("t", "Z", "cor.test"), verbose=FALSE, ...) {
   if(verbose) begin.time <- Sys.time()
   if(length(which.test)>1) which.test <- "t"
   xdim <- dim(x)
   tlen <- xdim[1]
   B.dof.test <- numeric(ntrials)
   if(which.test=="cor.test") cortester <- function(x,y,...) return(cor.test(x=x,y=y,...)$p.value)
   if(verbose) cat("\n", "Looping through ", ntrials, " times to simulate data and take correlations.  Enjoy!\n")
   for(i in 1:ntrials) {
	if(verbose & (i < 100 | i%%100==0)) cat(i, " ")
	z <- do.call(zfun, args=c(list(n=tlen), zfun.args))
	if(which.test=="cor.test") tmp <- apply(x, 2, cortester, y=z, ...)
	else {
	   cor.value <- abs(cor(x, z, use = "pairwise.complete.obs"))
	   if(which.test=="t") tmp <- sig.cor.t(cor.value, len=tlen, ...)
	   else if(which.test=="Z") tmp <- sig.cor.Z(cor.value, len=tlen, ...)
	}
	B.dof.test[i] <- mean(field.sig > tmp, na.rm = TRUE)
   } # end of for 'i' loop.
   if(verbose) print(Sys.time() - begin.time)
   return(list(MCprops=B.dof.test, minsigcov=quantile(B.dof.test, probs=1-field.sig, na.rm=TRUE)))
} # end of 'MCdof' function.

sig.cor.t <- function(r, len = 40, ...)
{
  #
  # A function to determine if a correlation is significantly different from x
  #
  # Input -
  # r - the (unsigned) correlation coefficient
  # len - the length of the vectors used to generate the correlation (default = 40)
  # Alpha - the significance level (default = 0.05)
  #
  # Output -
  # palpha.cor - the p-level of the correlation
  #
  # v1.0
  # KLE 2/3/2009
  # Ported to R -- 03/30/2011
  #
  t <- abs(r) * sqrt((len - 2)/(1 - r^2))
  palpha.cor <- 1 - pt(t, len - 2)
  return(palpha.cor)
} # end of 'sig.cor.t' function.

sig.cor.Z <- function(r, len = 40, H0 = 0)
{
  #
  # A function to find the significance of a correlation from H0 using Fisher's Z transform.
  #
  # Input -
  # r - the correlation
  #  len - the length of the crrelated vectors
  # Ho0- the null hypothesis correlation default = 0
  #
  # Output -
  # palpha.cor - the p-value of the correlation.
  #
  # KLE 02/20/2009
  # Ported to R -- 03/30/2011
  #
  W <- fisherz(abs(r))
  stderr <- 1/sqrt(len - 3)
  zscore <- (W - H0)/stderr
  palpha.cor <- 1 - pnorm(zscore)
  return(palpha.cor)
} # end of 'sig.cor.Z' function.

fisherz <- function(r)
{
  #
  # The sampling distribution of Pearson's r is not normally distributed. Fisher
  # developed a transformation now called "Fisher's z' transformation" that converts
  # Pearson's r's to the normally distributed variable z'.
  #
  # A function to perfoem Fisher's Z transformation on a correlation value, allowing
  # a t-test for significant correlation.
  #
  # Input -
  #  r - the correlation
  #
  # Output -
  #  W - the transformed corelation
  #
  # v1.0
  # KLE 2 Feb 2009
  # Ported to R -- 03/30/2011
  #
  W <- 0.5 * (log((1 + r)/(1 - r)))
  return(W)
}

plot.LocSig <- function(x, loc=NULL, nx=NULL, ny=NULL, ...){
  n <- dim(x)[1]
  if(is.null(loc)) {
     if(is.null(nx) | is.null(ny)) stop("plot.LocSig: must specify either loc or both nx and ny")
     loc <- cbind(rep(1:nx, ny), rep(1:ny, each=nx))
  } 
  mean.i <- as.image(x$Estimate, x=loc)
  thk.i <- as.image((x$Upper - x$Lower), x=loc)
  output <- list(mean.i = mean.i, thk.i = thk.i)
  par(mfrow=c(1,2))
  image.plot(mean.i, main="Mean of Estimate", ...)
  image.plot(thk.i, main="CI range of Estimate", ...)
  invisible(output)
}

plot.EBS <- function(x, ..., set.pw=FALSE, col, horizontal) {

    if(missing(col)) col <- c("gray", tim.colors(64))
    if(missing(horizontal)) horizontal <- TRUE

    if(set.pw) par(mfrow=c(1,2), oma=c(0,0,2,0))
    else par(oma=c(0,0,2,0))

    Zest <- x$block.boot.results$Estimate
    ZciR <- x$block.boot.results$Upper - x$block.boot.results$Lower

    a <- attributes(x)
    loc.byrow <- a$loc.byrow
    xd <- a$xdim

    if(is.null(a$subset)) {
        if(!is.matrix(Zest)) Zest <- matrix(Zest, xd[1], xd[2])
        if(!is.matrix(ZciR)) ZciR <- matrix(ZciR, xd[1], xd[2])
    } 

    if(a$projection && is.null(a$subset)) {
	xloc <- matrix(a$loc[,1], xd[1], xd[2], byrow=loc.byrow)
	yloc <- matrix(a$loc[,2], xd[1], xd[2], byrow=loc.byrow)
    }

    if(!is.null(a$subset)) {
	if(is.logical(a$subset)) Ns <- sum(a$subset, na.rm=TRUE)
	else Ns <- length(a$subset)
	Zest <- as.image(Zest, nx=ceiling(Ns/2), ny=ceiling(Ns/2), x=a$loc[a$subset,], na.rm=TRUE)
	ZciR <- as.image(ZciR, nx=ceiling(Ns/2), ny=ceiling(Ns/2), x=a$loc[a$subset,], na.rm=TRUE)
    } else if(!a$reg.grid) {
	Zest <- as.image(Zest, nx=xd[1], ny=xd[2], x=a$loc, na.rm=TRUE)
	ZciR <- as.image(ZciR, nx=xd[1], ny=xd[2], x=a$loc, na.rm=TRUE)
    }

    if(a$map) {

	ax <- list(x=pretty(round(a$loc[,1], digits=2)), y=pretty(round(a$loc[,2], digits=2)))

	if(is.null(a$subset)) r <- apply(a$loc, 2, range, finite=TRUE)
	else r <- apply(a$loc[a$subset,], 2, range, finite=TRUE)

	map(xlim=r[,1], ylim=r[,2], type="n")
	axis(1, at=ax$x, labels=ax$x)
	axis(2, at=ax$y, labels=ax$y)
	if(a$projection && a$reg.grid && is.null(a$subset)) image.plot(xloc, yloc, Zest, add=TRUE, col=col, main="Mean of Estimate", horizontal=horizontal, ...)
	else if(a$reg.grid && is.null(a$subset)) image.plot(Zest, add=TRUE, col=col, main="Mean of Estimate", horizontal=horizontal, ...)
	else image.plot(Zest, add=TRUE, col=col, main="Mean of Estimate", horizontal=horizontal, ...)
	map(add=TRUE, lwd=1.5)
	map(database="state", add=TRUE)

	map(xlim=r[,1], ylim=r[,2], type="n")
	axis(1, at=ax$x, labels=ax$x)
	axis(2, at=ax$y, labels=ax$y)
        if(a$projection && a$reg.grid && is.null(a$subset)) image.plot(xloc, yloc, ZciR, add=TRUE, col=col, main="CI Range", horizontal=horizontal, ...)
        else if(a$reg.grid && is.null(a$subset)) image.plot(ZciR, add=TRUE, col=col, main="Mean of Estimate", horizontal=horizontal, ...)
        else image.plot(ZciR, add=TRUE, col=col, main="CI Range", horizontal=horizontal, ...)
        map(add=TRUE, lwd=1.5)
        map(database="state", add=TRUE)

    } else {

	if(a$projection && a$reg.grid && is.null(a$subset)) image.plot(xloc, yloc, Zest, col=col, main="Mean of Estimate", horizontal=horizontal, ...)
        else if(a$reg.grid && is.null(a$subset)) image.plot(Zest, col=col, main="Mean of Estimate", horizontal=horizontal, ...)
        else image.plot(Zest, col=col, main="Mean of Estimate", horizontal=horizontal, ...)

	if(a$projection && a$reg.grid && is.null(a$subset)) image.plot(xloc, yloc, ZciR, col=col, main="CI Range", horizontal=horizontal, ...)
        else if(a$reg.grid && is.null(a$subset)) image.plot(ZciR, col=col, main="Mean of Estimate", ...)
        else image.plot(ZciR, col=col, main="CI Range", horizontal=horizontal, ...)

    } # end of if else 'map' stmts.

    if(length(a$data.name) == a$nforecast + 2) {
	msg <- paste(a$data.name[1], ": ", a$data.name[2], " vs ", a$data.name[a$which.model+2], sep="")
    } else msg <- paste(a$data.name[2], " vs ", a$data.name[a$which.model+1], sep="")
    if(a$field.type != "" && a$units != "") msg <- paste(msg, "\n", a$field.type, " (", a$units, ")", sep="")
    else if(a$field.type != "") msg <- paste(msg, a$field.type, sep="\n")
    else if(a$units != "") msg <- paste(msg, "\n(", a$units, ")", sep="")
    mtext(msg, line=0.05, outer=TRUE)

    invisible()
} # end of 'plot.EBS' function.

inside <- function(DF)
{
  #
  # A function to determine if the mean error at a data point is
  # inside the confidence limits.
  #
  # Input
  #    DF - data frame containing mean error values and upper and lower
  #  confidence limits
  #
  # Output
  #    logical T if outside, F if inside.
  #
  result <- !is.na(CI.fun(DF$Upper - DF$Lower, DF$Estimate))
  return(result)
} # end of 'inside' function.

CI.fun <- function(CI.field, est.field, replacement = NA)
{
  test <- ifelse((0.5 * CI.field < abs(est.field)), est.field, replacement)
  return(test)
} # end of 'CI.fun' function.

sig.coverage <- function(DF)
{
  tmp <- inside(DF)
  out <- sum(tmp,na.rm=TRUE)/(length(tmp[!is.na(tmp)])) - sum(is.na(DF$Estimate))
  return(out)
} # end of 'sig.coverage' function.

is.sig <- function(X, blockboot.results.df, n = 3000, fld.sig = 0.05, verbose=FALSE)
{
  #
  # Input -
  #  X - matrix of errors at gridpoints. n gridpoints by m days
  #  blockboot.results.df - dataframe of errors and CI at each
  #       gridpoint.
  #  n - number of Monte Carlo trials.
  #  field.sig - significnace alpha for field significance.
  #
  # Output -
  #  List
  #   name - Name of the data beig tested
  #   results - The minimum amount of areal coverage needed for
  #        significance at
  #             the given fld.sig
  #   actual - The actual coverage of sigificant gridpoint
  #        results.
  #   issig - Logical variable: T for significant results, F for
  #     non-significant results.
  #
  # KLE 01/2005
  #
  # Remove missing rows
  #
  sig.results <- MCdof(X, ntrials = n, field.sig = fld.sig, verbose=verbose)$minsigcov
  actual.coverage <- sig.coverage(blockboot.results.df)
  sig <- (actual.coverage > sig.results)
  output <- list(name = as.character(deparse(substitute(X))),
    required = as.numeric(sig.results), actual = as.numeric(
    actual.coverage), issig = as.logical(sig))
  return(output)
}

spatbiasFS <- function(X, Y, loc=NULL, block.length=NULL, alpha.boot=0.05, field.sig=0.05, bootR=1000, ntrials=1000, verbose=FALSE) {
   out <- list()
   if(!is.null(loc)) {
      data.name <- c(as.character(substitute(X)),as.character(substitute(Y)),as.character(substitute(loc)))
      names(data.name) <- c("verification","forecast","locations")
   } else {
      data.name <- c(as.character(substitute(X)),as.character(substitute(Y)))
      names(data.name) <- c("verification","forecast")
   }
   out$data.name <- data.name
   errfield <- Y - X
   hold <- LocSig(Z=errfield, numrep=bootR, block.length=block.length, alpha=alpha.boot)
   res <- is.sig(errfield, hold, n=ntrials, fld.sig=field.sig, verbose=verbose)
   out$block.boot.results <- hold
   out$sig.results <- res
   out$field.significance <- field.sig
   out$alpha.boot <- alpha.boot
   out$bootR <- bootR
   out$ntrials <- ntrials
   class(out) <- "spatbiasFS"
   return(out)
} # end of 'spatbiasFS' function.

summary.spatbiasFS <- function(object, ...) {
   cat("\n")
   msg <- paste("Results for ", object$data.name[2], " compared against ", object$data.name[1], sep="")
   print(msg)
   cat("\n", "\n")
   cat("Field significance level: ", object$field.significance, "\n")
   cat("Observed coverage of significant difference: ", object$sig.results$actual, "\n")
   cat("Required coverage for field significance: ", object$sig.results$required, "\n")
   invisible()
} # end of 'summary.spatbiasFS' function.

plot.spatbiasFS <- function(x, ...) {

# TO DO: Try to make a better plot (without using as.image).

   msg <- paste("Mean Error: ", x$data.name[2], " vs ", x$data.name[1], sep="")
   # X <- get(x$data.name[1])
   # Y <- get(x$data.name[2])
   if(length(x$data.name)==3) loc <- get(x$data.name[3])
   else stop("plot.spatbiasFS: No entry loc.  Must supply location information.")
   est.i <- as.image(x$block.boot.results$Estimate, x=loc)
   CIrange <- as.image(x$block.boot.results$Upper - x$block.boot.results$Lower, x=loc)
   par(mfrow=c(1,2))
   image.plot(est.i, col=tim.colors(64), axes=FALSE, main=msg)
   # map(add=TRUE)
   # map(add=TRUE,database="state")
   image.plot(CIrange, col=tim.colors(64), axes=FALSE, xlab=paste("Req. Coverage: ", round(x$sig.results$required,digits=2), " vs Obs. coverage: ",
			round(x$sig.results$actual,digits=2),  sep=""),
			main=paste((1-x$alpha.boot)*100, "% CI range", sep=""))
   invisible()
} # end of 'plot.spatbiasFS' function.
