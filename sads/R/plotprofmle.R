setGeneric("plotprofmle", 
    def=function(object, nseg=20, ratio=log(8), which=NULL, ask=NULL, col.line="blue", varname=NULL, ...) standardGeneric("plotprofmle")
    )
setMethod("plotprofmle", "profile.mle2",
function(object, nseg, ratio, which, ask, col.line, varname, ...){
  mleprof <- object@profile
  npar <- length(mleprof)
  if(missing(which))
    which <- 1:npar
  if(missing(ask))
    ask <- (prod(par("mfcol")) < length(which)) && dev.interactive()
  dots <- list(...)
  if(!"ylab" %in% names(dots)) dots$ylab <- "Negative relative log-likelihood"
  if(!"type" %in% names(dots)) dots$type <- "l"
  if(!"col" %in% names(dots)) dots$col <- "red"
  vname <- names(mleprof)
  if( is.null(which) ){
    if(!missing(varname)){
      if(length(varname)!=length(mleprof))stop("Length of 'varname' should match number os mles in mle.prof object")
      vname <- varname
    }
    parseq = 1:npar
  }
  else{
    if(!missing(varname)){
      if(length(varname)!=length(which))stop("Length of 'which' should match length of 'varnames'")
      vname[which] <- varname
    }
    parseq = which
  }
  if(ask){
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  for(i in parseq)
    {
      tmp <- mleprof[i][[1]]
      y <- tmp[,1]^2/2
      x <- (tmp[,2][,i])
      interpol = spline(x, y, n=nseg*length(x) )
      do.call(plot, c(list(x=interpol,xlab=vname[i]),dots))
      if(!is.null(ratio)){
        l <- length(interpol$y)
	  # Finds where the interpolation crosses the "y = ratio" line
        change <- (interpol$y - ratio)[2:l] * (interpol$y - ratio)[1:(l-1)]
        endpoints <- which(change < 0)
		# Adds the borders, if any of them is lower than ratio
		if(interpol$y[1] < ratio) endpoints <- c(1, endpoints)
		if(interpol$y[l] < ratio) endpoints <- c(endpoints, l)
        corr <- (interpol$x[2]-interpol$x[1])/2
		if(length(endpoints) > 0)
        for (j in 1:(length(endpoints)/2)) {
          lower <-interpol$x[endpoints[(2*j)-1]]+corr
          upper <- interpol$x[endpoints[2*j]]+corr
          lines(c(lower,upper ),c(ratio, ratio), col=col.line, lty=2)
		  if(endpoints[(2*j-1)] != 1) # dont draw vertical lines at the borders
	          lines(rep(lower,2), c(-1, ratio), col=col.line, lty=2)
		  if(endpoints[(2*j)] != l) # dont draw vertical lines at the borders
			  lines(rep(upper,2), c(-1, ratio), col=col.line, lty=2)
        }
      }
    }
})
setMethod("plotprofmle", "mle2",
    function(object, ...) {
    cat("NOTICE: Running a profile on the object. You should consider storing the profile\n")
    cat("in a different variable\n")
    plotprofmle(profile(object), ...)
})
