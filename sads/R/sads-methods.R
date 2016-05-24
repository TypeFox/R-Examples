setMethod("plot", "rad",
          function(x, ...){
            dots <- list(...)
            if(!"log" %in% names(dots)) dots$log <- "y"
            if(!"xlab" %in% names(dots)) dots$xlab = "Species Rank"
            if(!"ylab" %in% names(dots)) dots$ylab = "Species Abundance"
            if(!"frame.plot" %in% names(dots)) dots$frame.plot = TRUE
            if(!"axes" %in% names(dots)){ 
              do.call(plot, c(list(x = x[, 1], y = x[, 2], axes=FALSE), dots))
              axis(2)
              sc <- axisTicks(range(x[, 1]),nint=10,log=FALSE)
              sc[sc==0] <- 1
              axis(1,at=sc)
            }
            if("axes" %in% names(dots)){ 
              do.call(plot, c(list(x = x[, 1], y = x[, 2]), dots))
            }
            
          }
            )

setMethod("points", "rad",
          function(x, ...){
            dots <- list(...)
            if(!"type" %in% names(dots)) dots$type = "p"
            if(!"col" %in% names(dots)) dots$col = "blue"
            do.call(points, c(list(x = x[, 1], y = x[, 2]), dots)) 
          }
          )

setMethod("lines", "rad",
          function(x, ...){
            dots <- list(...)
            if(!"type" %in% names(dots)) dots$type = "l"
            if(!"col" %in% names(dots)) dots$col = "blue"
            do.call(lines, c(list(x = x[, 1], y = x[, 2]), dots)) 
          }
)

setMethod("plot","octav",
          function(x, prop=FALSE, x.oct=FALSE, par.axis=list(), ...){
            dots <- list(...)
            x.hist <- rep(as.integer(as.character(x$octave)), as.integer(as.character(x$Freq)))
            h1 <- hist(x=x.hist,
                           breaks = c((min(as.integer(as.character(x$octave)))-1),as.integer(as.character(x$octave))),
                           plot=FALSE)
            if(prop) h1$counts <- h1$counts/sum(h1$counts)
            if(x.oct) xlab <- x[seq(1,length(x[,1]),2),1]
            if(!x.oct) xlab <- x[seq(1,length(x[,1]),2),2]
            if(!"col" %in% names(dots)) dots$col = "gray"
            if(!"main" %in% names(dots)) dots$main = ""
            if(!"ylab" %in% names(dots) & !prop) dots$ylab = "N of species"
            if(!"ylab" %in% names(dots) & prop) dots$ylab = "Proportion of species"
            if(!"xlab" %in% names(dots) & !x.oct) dots$xlab = "Abundance class"
            if(!"xlab" %in% names(dots) & x.oct) dots$xlab = "Abundance class (log2)"
            if(!"axes" %in% names(dots)){
                do.call(plot, c(list(x=h1, axes=FALSE),dots))
                n <- as.numeric(as.character(x[,1]))
                do.call(axis, c(list(side=2), par.axis))
                do.call(axis, c(list(side=1,at=n[seq(1,length(x[,1]),2)],
                     labels=xlab),par.axis))
            }
            else
              do.call(plot, c(list(x=h1,dots)))
          }
          )

setMethod("points","octav",
          function(x, prop=FALSE, ...){
            dots <- list(...)
            if(!"type" %in% names(dots)) dots$type="b"
            if(!"col" %in% names(dots)) dots$col="blue"
            X <- c((min(as.integer(as.character(x$octave)))-1), as.integer(as.character(x$octave)))
            X <- X[-length(X)]+diff(X)/2
            if(prop) Y <- x$Freq/sum(x$Freq)
            if(!prop) Y <- x$Freq
            do.call(points, c(list(x = X, y = Y), dots))
          }
          )

setMethod("lines","octav",
          function(x, prop=FALSE, ...){
            dots <- list(...)
            if(!"type" %in% names(dots)) dots$type="b"
            if(!"col" %in% names(dots)) dots$col="blue"
            X <- c((min(as.integer(as.character(x$octave)))-1), as.integer(as.character(x$octave)))
            X <- X[-length(X)]+diff(X)/2
            if(prop) Y <- x$Freq/sum(x$Freq)
            if(!prop) Y <- x$Freq
            do.call(lines, c(list(x = X, y = Y), dots))
          }
)

setMethod("plot","fitsad",
          function(x, which=1:4, ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...){
            oct.df <- octav(x)
            rad.df <- rad(x)
            oct.pred <- octavpred(x)
            oct.ymax <- max(c(oct.df[, 3], oct.pred[, 3]), na.rm = TRUE)
            rad.pred <- radpred(x)
            rad.ylim <- range(c(rad.df[, 2], rad.pred[, 2]), na.rm = TRUE)
            if (ask) {
              oask <- devAskNewPage(TRUE)
              on.exit(devAskNewPage(oask))
            }
            if(1 %in% which){
              plot(oct.df, ylim = c(0, oct.ymax), ...)
              points(oct.pred, ...)
            }
            if(2 %in% which){
              plot(rad.df, ylim = rad.ylim, ...)
              lines(rad.pred, ...)
            }
            if(3 %in% which){
              qqsad(x, ...)
            }
            if(4 %in% which){
              ppsad(x, ...)
            }
          }
          )

setMethod("plot","fitrad",
          function(x, which=1:4, ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...){
            oct.df <- octav(x)
            rad.df <- rad(x)
            oct.pred <- octavpred(x)
            oct.ymax <- max(c(oct.df[, 3], oct.pred[, 3]), na.rm = TRUE)
            rad.pred <- radpred(x)
            rad.ylim <- range(c(rad.df[, 2], rad.pred[, 2]), na.rm = TRUE)
            if (ask) {
              oask <- devAskNewPage(TRUE)
              on.exit(devAskNewPage(oask))
            }
            if(1 %in% which){
              plot(oct.df, ylim = c(0, oct.ymax), ...)
              points(oct.pred, ...)
            }
            if(2 %in% which){
              plot(rad.df, ylim = rad.ylim, ...)
              lines(rad.pred, ...)
            }
            if(3 %in% which){
              qqrad(x, ...)
            }
            if(4 %in% which){
              pprad(x, ...)
            }
          }
          )

## copy of the methods in bbmle 1.0.17, tweaked to work better with fitsad/fitrad classed
## which do not have an explicit df attribute. Also fixes some inconsistencies
## in the handling of parameters
setMethod("AIC", "mle2",
		  function (object, ..., k = 2) {
			  L <- list(object, ...)
			  if (!all(sapply(L, function(x) inherits(x, "mle2")))) 
				  stop("all objects in list must be class mle2 or inherit from mle2")
			  if (length(L) > 1) {
				  logLiks <- lapply(L, logLik)
				  AICs <- sapply(logLiks,AIC,k=k)
				  df <- sapply(logLiks,attr,"df")
			  data.frame(AIC=AICs,df=df)
			  } else AIC(logLik(object), k=k)
		  })

setMethod("AICc","mle2",
		  function (object, ..., nobs, k = 2){
			  L <- list(object, ...)
			  if (!all(sapply(L, function(x) inherits(x, "mle2")))) 
				  stop("all objects in list must be class mle2 or inherit from mle2")
			  if (missing(nobs)) {
				  nobs <- sapply(L,match.fun("nobs"))
			  }
			  if (length(L) > 1) {
				  if (length(unique(nobs)) > 1) 
					  stop("nobs different: must have identical data for all objects")
				  logLiks <- lapply(L, logLik)
				  df <- sapply(logLiks,attr,"df")
				  val <- -2*unlist(logLiks)+k*df+2*df*(df+1)/(nobs-df-1)
				  data.frame(AICc = val, df = df)
			  } else {
				  df <- attr(logLik(object), "df")
				  c(-2 * logLik(object)+k*df+2*df*(df+1)/(nobs-df-1))
			  }
		  }
          )

setMethod("nobs", "fitsad",
		  function(object) length(object@data$x)
		  )
setMethod("nobs", "fitrad",
		  function(object) length(object@rad.tab$abund)
		  )

### Copy of functions from mle2, including some error-checking, slots specific to fitrad/fitsad and
### truncating the display of the call
showmle2 <- function(object) {
    cat("Maximum likelihood estimation\nType: ")
	if (inherits(object, "fitsad")) {
		cat (distr(object@sad), " species abundance distribution")
		my.x <- object@data$x
	}
	else {
		cat (distr(object@rad), "rank abundance distribution")
		my.x <- object@rad.tab$abund
	}
	cat("\nSpecies:",length(my.x),"individuals:", sum(my.x), "\n")
    cat("\nCall:\n")
	# Summarizes the call to avoid printing pages of data
	d <- object@call.orig$data$x
	if (length(d) > 6) { 
		d <- c(as.list(as.numeric(d[1:5])), "etc")
		object@call.orig$data$x <- d
	}
	print(object@call.orig)
    cat("\nCoefficients:\n")
    print(coef(object))
    if(!is.nan(object@trunc)) {
		cat(paste("\nTruncation point:", object@trunc, "\n"))
	}
    cat("\nLog-likelihood: ")
    cat(round(as.numeric(logLik(object)),2),"\n")
    if (object@optimizer=="optimx" && length(object@method)>1) {
      cat("Best method:",object@details$method.used,"\n")
    }
	if (!is.null(object@details$convergence))
		if(object@details$convergence > 0)
	      cat("\nWarning: optimization did not converge (code ",
          object@details$convergence,": ",object@details$message,")\n",sep="")
  }
setMethod("show", "fitsad", function(object){showmle2(object)})
setMethod("show", "fitrad", function(object){showmle2(object)})

#### summary class dealing with fixed parameters (such as fitls, fitvolkov, etc)
#' @rdname summary.sads-class
#' @param object An object of class fitsad/fitrad is required to generate a summary.sads object.
setMethod("show", "summary.sads", function(object){
          cat("Maximum likelihood estimation\n\nCall:\n")
          print(object@call)
          cat("\nCoefficients:\n")
          printCoefmat(object@coef)
          if (length(object@fixed) > 0) {
            cat("\nFixed parameters:\n")
            print(object@fixed)
          }
          cat("\n-2 log L:", object@m2logL, "\n")
          })
sumle2 <- function(object, ...){
  cmat <- cbind(Estimate = object@coef,
                `Std. Error` = sqrt(diag(object@vcov)))
  zval <- cmat[,"Estimate"]/cmat[,"Std. Error"]
  pval <- 2*pnorm(-abs(zval))
  coefmat <- cbind(cmat,"z value"=zval,"Pr(z)"=pval)
  m2logL <- 2*object@min
  fixed <- numeric()
  if (! all(object@fullcoef %in% object@coef))
    fixed <- object@fullcoef [! object@fullcoef %in% object@coef]
  new("summary.sads", call=object@call.orig, coef=coefmat, fixed=fixed, m2logL= m2logL)
}
#' @rdname summary.sads-class
setMethod("summary", "fitsad", function(object){sumle2(object)})
#' @rdname summary.sads-class
setMethod("summary", "fitrad", function(object){sumle2(object)})

## radpred generic functions and methods ###
setGeneric("radpred",
def = function(object, sad, rad, coef, trunc , distr=NA, S, N) standardGeneric("radpred")
           )

## if object is of class fitsad (no other argument should be provided)
# Extracts information from object and uses method below
setMethod("radpred",signature(object="fitsad", sad="missing", rad="missing",
                              coef="missing", trunc="missing", distr="missing", S="missing", N="missing"),
          function (object){
			  ab = object@data$x
			  radpred(sad=object@sad, coef=as.list(bbmle::coef(object)),
					  trunc=object@trunc, S=length(ab), N=sum(ab))
		  }
		  )

## if object is of class fitrad (no other argument should be provided)
# Extracts information from object and uses method below
setMethod("radpred",signature(object="fitrad", sad="missing", rad="missing",
                              coef="missing", trunc="missing", distr="missing", S="missing", N="missing"),
          function(object){
			  ab = object@rad.tab$abund
			  radpred(rad=object@rad, coef=as.list(bbmle::coef(object)), 
					  trunc=object@trunc, S=length(ab), N=sum(ab))
		  }
		  )

## if object is a numeric vector of abundances and rad argument is given (sad, S, N, distr,  arguments should be missing)
# Extracts information from object and uses method below
setMethod("radpred",signature(object="numeric", sad="missing", rad="character",
                              coef="list", trunc="ANY", distr="missing", S="missing", N="missing"),
          function(object, sad, rad, coef, trunc){
			  if(missing(trunc)) trunc <- NaN
			  radpred(rad=rad, coef=coef, trunc=trunc, S=length(object), N= sum(object))
		  }
		  )

## if object is a numeric vector of abundances and sad argument is given (rad, S, N,  arguments should be missing)
setMethod("radpred",signature(object="numeric", sad="character", rad="missing",
                              coef="list", trunc="ANY", distr="ANY", S="missing", N="missing"),
          function(object, sad, rad, coef, trunc, distr=NA){
        if(!is.na(distr)) warning("The parameter distr has been deprecated and is ignored, see ?distr")
			  if(missing(trunc)) trunc <- NaN
			  radpred(sad=sad, coef=coef, trunc=trunc, S=length(object), N= sum(object))
		  }
		  )

## if object is missing and rad is given. sad should not be given. All other arguments except distr should be given,
## except trunc (optional). This is the base method for all signatures using "rad" or "fitrad" 
setMethod("radpred", signature(object="missing", sad="missing", rad="character",
                              coef="list", trunc="ANY", distr="missing", S="numeric", N="numeric"),
          function(object, sad, rad, coef, trunc, distr, S, N){
            y <- 1:S
            if(missing(trunc)) trunc <- NaN
            if(!is.nan(trunc))
              ab <- do.call(dtrunc, c(list(rad, x = y, coef = coef, trunc = trunc)))*N
            else{
              drad <- get(paste("d", rad, sep=""),  mode = "function")
              ab <- do.call(drad, c(list(x = y), coef))*N
            }
            new("rad", data.frame(rank=1:S, abund=ab))
          }
          )

## if object is missing and sad is given. rad should not be given.
## All other arguments except distr should be given, except trunc (optional)
# This is the base method for all signatures using "sad" or "fitsad" 
setMethod("radpred", signature(object="missing", sad="character", rad="missing",
                               coef="list", trunc="ANY", distr="ANY", S="numeric", N="numeric"),
          function(object, sad, rad, coef, trunc, distr=NA, S, N){
              if(!is.na(distr)) warning("The parameter distr has been deprecated and is ignored, see ?distr")
              distribution <- distr(sad)
              if(missing(trunc)) trunc <- NaN
              if (distribution == "discrete"){
                  ## Approximates the [q] function instead of calling it directly to save some
                  ## computational time (as [q] is inneficiently vectorized)
                  y <- 1:N
                  Y <- ppoints(S)
                  if(!is.nan(trunc))
                    X <- do.call(ptrunc, list(sad, q = y, coef = coef, lower.tail=F, trunc = trunc))
                  else {
                    psad <- get(paste("p", sad, sep=""), mode = "function")
                    qsad <- get(paste("q", sad, sep=""), mode = "function")
                    X <- do.call(psad, c(list(q = y, lower.tail = F), coef))
                  }
                  ab <- approx(x=c(1, X), y=c(0, y), xout=Y, method="constant")$y
                  ## Extreme values of abundance are out of bounds for approx. Explicit form:
                  for (i in 1:length(ab)) {
                      if (!is.na(ab[i])) break;
                      cat("Note: extreme values generated by radpred. Calculations will take a while...\n")
                      if(! is.nan(trunc))
                        ab[i] <- do.call(qtrunc, list(sad, p = Y[i], coef = coef, lower.tail=FALSE, trunc = trunc))
                      else
                        ab[i] <- do.call(qsad, c(list(p = Y[i], lower.tail=FALSE), coef))
                  }
              }
              else if(distribution == "continuous"){
                Y <- ppoints(S)
                if(!is.nan(trunc))
                  ab <- do.call(qtrunc, list(sad, p = Y, coef = coef, lower.tail=F, trunc = trunc))
                else{
                  qsad <- get(paste("q", sad, sep=""), mode = "function")
                  ab <- do.call(qsad, c(list(p = Y, lower.tail = F), coef))
                }
              } else
                stop("Please provide a valid distribution") 
                new("rad", data.frame(rank=1:S, abund=ab))
          }
          )

# Helper function for octavpred
genoct <- function (x) {
  oct <- 0:(ceiling(max(log2(x)))+1)
  if(any(x < 1)){
    octlower <- floor(min(log2((x)))):-1
    oct <- c(octlower, oct)
  }
  oct
}
          
## octavpred generic functions and methods ###
setGeneric("octavpred",
def = function(object, sad, rad, coef, trunc, oct, S, N, preston=FALSE, ...) standardGeneric("octavpred"))

## if object is of class fitsad (no other argument should be provided)
setMethod("octavpred", signature(object="fitsad",sad="missing", rad="missing",
                                 coef="missing", trunc="missing", oct="ANY",
                                 S="missing", N="missing"),
          function(object, sad, rad, coef, trunc, oct, S, N, preston, ...){
            x <- object@data$x
            if(missing(oct)) oct <- genoct(x)
            octavpred(sad = object@sad, coef = as.list(bbmle::coef(object)),
                      trunc = object@trunc, oct = oct, S=length(x), N=sum(x), preston=preston, ...)
          }
          )
## if object is a numeric vector of abundances and sad argument is given (rad, S, N,  arguments should be missing)
setMethod("octavpred", signature(object="numeric",sad="character", rad="missing",
                                 coef="list", oct="ANY", trunc="ANY", S="missing", N="missing"),
          function(object, sad, rad, coef, trunc, oct, S, N, preston, ...){
            if(missing(oct)) oct <- genoct(object)
            if(missing(trunc)) trunc<-NaN
            octavpred(sad=sad, coef=coef, trunc=trunc, oct=oct, S = length(object), N = sum(object),
                      preston=preston, ...)
          }
          )
## Octavpred workhorse for "sads"
setMethod("octavpred", signature(object="missing",sad="character", rad="missing",
                                 coef="list", trunc="ANY", oct="ANY", S="numeric", N="numeric"),
          function(object, sad, rad, coef, trunc, oct, S, N, preston, ...){
            dots <- list(...)
            if(missing(oct)) oct <- genoct(N)
            if(missing(trunc)) trunc <- NaN
            oct <- unique(oct)
            if (preston) {
              return(octav(radpred(sad=sad, coef=coef, trunc=trunc, distr=NA, S=S, N=N)$abund, preston=TRUE))
            } else {
              n <- 2^oct
              if(!is.nan(trunc)){
                Y <- do.call(ptrunc, c(list(sad, q = n, coef = coef, trunc = trunc), dots))
              }
              else{
                psad <- get(paste("p",sad,sep=""),mode="function")
                Y <- do.call(psad, c(list(q = n),coef,dots))
              }
              Y <- c(Y[1], diff(Y))*S
              return(new("octav", data.frame(octave = oct, upper = n, Freq = Y)))
            }
          }
          )

## if object is of class fitrad (no other argument should be provided, except oct (optional))
setMethod("octavpred", signature(object="fitrad",sad="missing", rad="missing",
                                 coef="missing", trunc="missing", oct="ANY",
                                 S="missing", N="missing"),
          function(object, sad, rad, coef, trunc, oct, S, N, preston, ...){
            x <- object@rad.tab$abund
            if(missing(oct)) oct <- genoct(x)
            octavpred(rad = object@rad, coef = as.list(bbmle::coef(object)),
                      trunc = object@trunc, oct = oct, S=length(x), N=sum(x),
                      preston=preston, ...)
          }
          )
## if object is a numeric vector of abundances and rad argument is given (sad, S, N,  arguments should be missing)
setMethod("octavpred", signature(object="numeric",sad="missing", rad="character",
                                 coef="list", trunc="ANY", oct="ANY", S="missing", N="missing"),
          function(object, sad, rad, coef, trunc, oct, S, N, preston, ...){
            if(missing(oct)) oct <- genoct(object)
            if(missing(trunc)) trunc<-NaN
            octavpred(rad=rad, coef=coef, trunc=trunc, oct=oct, S = length(object), N = sum(object),
                      preston=preston, ...)
          }
)
## Octavpred workhorse for "rads"
setMethod("octavpred", signature(object="missing",sad="missing", rad="character",
                                 coef="list", trunc="ANY", oct="ANY", S="numeric", N="numeric"),
          function(object, sad, rad, coef, trunc, oct, S, N, preston, ...){
            dots <- list(...)
            if(missing(oct)) oct <- genoct(N)
            if(missing(trunc)) trunc<-NaN
            oct <- unique(oct)
            n <- 2^oct
            if(!is.nan(trunc)){
              ab <- do.call(dtrunc, c(list(f=rad, q = 1:S, coef=coef,trunc = trunc),dots))*N
            }
            else{
              drad <- get(paste("d",rad,sep=""),mode="function")
              ab <- do.call(drad, c(list(x=1:S),coef,dots))*N
            }
            tryCatch({Y = hist(ab, breaks=c(2^(min(oct)-2),n), plot=FALSE)},
                     error = function(cond) stop("Octaves do not span the entire range, try using a larger oct argument (maybe negative octaves?)")
                     )
            res <- data.frame(octave = oct, upper = n, Freq = Y$count)
            if(preston) res <- prestonfy(res, ceiling(ab))
            new("octav", res)
          }
)

## Generic and methods for qqsad
setGeneric("qqsad",
def = function(x, sad, coef, trunc=NA, distr=NA, plot=TRUE, line=TRUE, ...) standardGeneric("qqsad"))

## method for class numeric
## if x is numeric (abundances), all other arguments should be given.
## Only trunc, plot and line are optional because they have default values
setMethod("qqsad",
          signature(x="numeric", sad="character", coef="list", distr="ANY"),
          function(x, sad, coef, trunc=NA, distr=NA, plot=TRUE, line=TRUE, ...){
        if(!is.na(distr)) warning("The parameter distr has been deprecated and is ignored, see ?distr")
        distribution <- distr(sad)
              x.sorted <- sort(x)
              S <- length(x)
              if(distribution == "discrete"){
                  q <- 1:sum(x)
                  if(!is.na(trunc)){
                      p <- do.call(ptrunc, list(sad, q = q, coef=coef, trunc=trunc))
                  }
                  else{
                      psad <- get(paste("p", sad, sep=""), mode = "function")
                      p <- do.call(psad, c(list(q = q), coef))
                  }
                  f1 <- approxfun(x=c(1, p), y=c(0, q), method="constant")
                  q <- f1(ppoints(S))
              }
        else if(distribution == "continuous"){
            p <- ppoints(S)
            if(!is.na(trunc))
                q <- do.call(qtrunc, list(sad, p = p, trunc = trunc, coef=coef))
            else{
                qsad <- get(paste("q", sad, sep=""), mode = "function")
                q <- do.call(qsad, c(list(p = p), coef))
            }
        }
        else
            stop("Please provide a valid distribution") 
        if(plot){
            dots <- list(...)
            if(!"main" %in% names(dots)) dots$main = "Q-Q plot"
            if(!"xlab" %in% names(dots)) dots$xlab = "Theoretical Quantile"
            if(!"ylab" %in% names(dots)) dots$ylab = "Sample Quantiles"
            do.call(graphics::plot, c(list(x=q, y=x.sorted),dots))
            if(line) abline(0, 1, col = "red", lty = 2)
        }
        return(invisible(data.frame(theoret.q=q, sample.q=x.sorted)))
          }
        )


## If x is of the class fitsad all other arguments should be ommited
## plot and line have default values and are optional
setMethod("qqsad",
          signature(x="fitsad", sad="missing", coef="missing",
                    trunc="missing", distr="missing"),
          function(x, sad, coef, trunc, distr, plot=TRUE, line=TRUE, ...){
              qqsad(x=x@data$x, sad=x@sad, coef=as.list(bbmle::coef(x)), 
                    trunc=x@trunc, plot=plot, line=line, ...)
          }
          )


## Generic and methods for qqrad
setGeneric("qqrad",
def = function(x, rad, coef, trunc=NA, plot=TRUE, line=TRUE, ...) standardGeneric("qqrad"))

## If x is an object of class rad
setMethod("qqrad",
          signature(x="rad", rad="character", coef="list"),
          function(x, rad , coef, trunc=NA, plot=TRUE, line=TRUE, ...){
              pr <- cumsum(x$abund/sum(x$abund))
              if(!is.na(trunc))
                  q <- do.call(qtrunc, list(rad, p = pr, coef = coef, trunc = trunc))
              else{
                  qrad <- get(paste("q", rad, sep=""), mode = "function")
                  q <- do.call(qrad, c(list(p = pr), coef))
              }
              if(plot){
                  dots <- list(...)
                  if(!"main" %in% names(dots)) dots$main = "Q-Q plot"
                  if(!"xlab" %in% names(dots)) dots$xlab = "Theoretical Quantile"
                  if(!"ylab" %in% names(dots)) dots$ylab = "Sample Quantiles"
                  do.call(graphics::plot, c(list(x=q, y=x$rank),dots))
                  if(line) abline(0, 1, col = "red", lty = 2)
              }
              return(invisible(data.frame(theoret.q=q, sample.q=x$rank)))
          }
          )

## If object is of class numeric arguments rad and coef should be provided
setMethod("qqrad",
          signature(x="numeric", rad="character", coef="list"),
          function(x, rad , coef, trunc=NA, plot=TRUE, line=TRUE, ...){
              y <- rad(x)
              qqrad(x=y, rad=rad, coef=coef, trunc=trunc, plot=plot, line=line, ...)
          }
          )

## If object is of class integer arguments rad and coef should be provided
## setMethod("qqrad",
##           signature(x="integer", rad="character", coef="list",
##                     trunc="ANY", plot="ANY", line="ANY"),
##           function(x, rad , coef, trunc=NA, plot=TRUE, line=TRUE, ...){
##               y <- as.numeric(x)
##               qqrad(x=y, rad=rad, coef=coef, trunc=trunc, plot=plot, line=line, ...)
##           }
##           )

## If object is of class fitrad arguments rad or coef should be missing
setMethod("qqrad",
          signature(x="fitrad", rad="missing", coef="missing", trunc="missing"),
          function(x, rad , coef, trunc, plot=TRUE, line=TRUE, ...){
              rad <- x@rad
              coef <- as.list(bbmle::coef(x))
              trunc <- x@trunc
              y <- x@rad.tab
              qqrad(x=y, rad=rad, coef=coef, trunc=trunc, plot=plot, line=line, ...)
          }
          )


## Generic function and methods for ppsad ##
setGeneric("ppsad",
def = function(x, sad, coef, trunc=NA, plot=TRUE, line=TRUE, ...) standardGeneric("ppsad"))

## If x is numeric arguments sad and coef should be provided
setMethod("ppsad",
          signature(x="numeric", sad="character", coef="list"),
          function (x, sad, coef, trunc=NA, plot=TRUE, line=TRUE, ...) {
              x.sorted <- sort(x)
              S <- length(x)
              z <- ppoints(S)
              if(!is.na(trunc)){
				  p <- do.call(ptrunc, list(sad, q = x.sorted, coef = coef, trunc = trunc))
              }
			  else{
				  psad <- get(paste("p", sad, sep=""), mode = "function")
				  p <- do.call(psad, c(list(q = x.sorted), coef))
              }
              if(plot){
                  dots <- list(...)
                  if(!"main" %in% names(dots)) dots$main = "P-P plot"
                  if(!"xlab" %in% names(dots)) dots$xlab = "Theoretical Percentiles"
                  if(!"ylab" %in% names(dots)) dots$ylab = "Sample Percentiles"
                  do.call(graphics::plot, c(list(x=p, y=z, ylim=c(0,1)),dots) )
                  if(line) abline(0, 1, col = "red", lty = 2)
              }
              return(invisible(data.frame(theoret.p=p, sample.p=z)))
          }
          )

## If object is of class integer arguments rad and coef should be provided
## setMethod("ppsad",
##           signature(x="integer", sad="character", coef="list",
##                     trunc="ANY", plot="ANY", line="ANY"),
##           function(x, sad , coef, trunc=NA, plot=TRUE, line=TRUE, ...){
##               y <- as.numeric(x)
##               ppsad(x=y, sad=sad, coef=coef, trunc=trunc, plot=plot, line=line, ...)
##           }
##           )

## If argument x is fitsad class, arguments sad and coef should be missing
setMethod("ppsad",
          signature(x="fitsad", sad="missing", coef="missing", trunc="missing"),
          function (x, sad, coef, trunc=NA, plot=TRUE, line=TRUE, ...) {          
              sad <- x@sad
              coef <- as.list(bbmle::coef(x))
              trunc <- x@trunc
              y <- x@data$x
              ppsad(x=y, sad=sad, coef=coef, trunc=trunc, plot=plot, line=line, ...)
          }
          )

## Generic function and methods for pprad ##
setGeneric("pprad",
def = function(x, rad, coef, trunc=NA, plot=TRUE, line=TRUE, ...) standardGeneric("pprad"))

## If argument is of class rad arguments rad and coef should be provided
setMethod("pprad",
          signature(x="rad", rad="character", coef="list"),
          function (x, rad, coef, trunc=NA, plot=TRUE, line=TRUE, ...) {
              rad.tab <- x
              pr <- cumsum(rad.tab$abund/sum(rad.tab$abund))
              if(!is.na(trunc)){
                  p <- do.call(ptrunc, list(rad, q = rad.tab$rank, coef = coef, trunc = trunc))
              }
              else{
                  prad <- get(paste("p", rad, sep=""), mode = "function")
                  p <- do.call(prad, c(list(q = rad.tab$rank), coef))
              }
              if(plot){
                  dots <- list(...)
                  if(!"main" %in% names(dots)) dots$main = "P-P plot"
                  if(!"xlab" %in% names(dots)) dots$xlab = "Theoretical Percentiles"
                  if(!"ylab" %in% names(dots)) dots$ylab = "Sample Percentiles"
                  do.call(graphics::plot, c(list(x=p, y=pr, ylim=c(0,1)),dots) )
                  if(line) abline(0, 1, col = "red", lty = 2)
              }
              return(invisible(data.frame(theoret.p=p, sample.p=pr)))
          }
)

## If argument is of class numeric arguments rad and coef should be provided
setMethod("pprad",
          signature(x="numeric", rad="character", coef="list"),
          function (x, rad, coef, trunc=NA, plot=TRUE, line=TRUE, ...) {
              y <- rad(x)
              pprad(x=y, rad=rad, coef=coef, trunc=trunc, plot=plot, line=line, ...)
          }
          )

## If argument is of class integer arguments rad and coef should be provided
## setMethod("pprad",
##           signature(x="integer", rad="character", coef="list",
##                     trunc="ANY", plot="ANY", line="ANY"),
##           function (x, rad, coef, trunc=NA, plot=TRUE, line=TRUE, ...) {
##               y <- as.numeric(x)
##               pprad(x=y, rad=rad, coef=coef, trunc=trunc, plot=plot, line=line, ...)
##           }
##           )

## If argument is of class fitrad arguments rad and coef should be missing
setMethod("pprad",
          signature(x="fitrad", rad="missing", coef="missing"),
          function (x, rad, coef, trunc=NA, plot=TRUE, line=TRUE, ...) {
              rad <- x@rad
              coef <- as.list(bbmle::coef(x))
              trunc <- x@trunc
              y <- x@rad.tab
              pprad(x=y, rad=rad, coef=coef, trunc=trunc, plot=plot, line=line, ...)
          }
          )

### Providing standard stats methods
#' Standard stats methods
#' 
#' Provide the standard interface for fitted objects
#' 
#' These methods are provided to allow for standard manipulation of \code{\link{fitsad}}
#' and \code{\link{fitrad}} objects using the generic methods defined in the "stats" package.
#' Please see the original man pages for each method.
#' 
#' \code{coefficients} is an alias to \code{\link[stats]{coef}} (implemented in package "bbmle").
#' 
#' \code{fitted} and \code{fitted.values} provide an alternative interface to \code{\link{radpred}};
#' these are also used to calcutate \code{residuals}.
#' 
#' Notice that radpred is a preferred interface for most calculations, specially if there are several
#' ties.
#' 
#' @param object An object from class fitsad or fitrad
#' @param \dots Other arguments to be forwarded for the lower level function
#' @rdname stats
setMethod("coefficients", signature(object="fitsad"),
          function(object, ...) bbmle::coef(object, ...))
#' @rdname stats
setMethod("coefficients", signature(object="fitrad"),
          function(object, ...) bbmle::coef(object, ...))
#' @rdname stats
setMethod("fitted.values", signature(object="fitsad"),
          function(object, ...) fitted(object, ...))
#' @rdname stats
setMethod("fitted.values", signature(object="fitrad"),
          function(object, ...) fitted(object, ...))
#' @rdname stats
setMethod("fitted", signature(object="fitsad"),
          function(object, ...) {
            rad <- radpred(object)$abund
            rad <- rad[rev(order(object@data$x))]
            names(rad) <- as.character(1:length(rad))
            return(rad)
          }
          )
#' @rdname stats
setMethod("fitted", signature(object="fitrad"),
          function(object, ...) {
            rad <- radpred(object)$abund
            rad <- rad[order(object@rad.tab$abund)]
            names(rad) <- as.character(1:length(rad))
            return(rad)
          }
          )
#' @rdname stats
setMethod("residuals", signature(object="fitsad"),
          function(object, ...) object@data$x - fitted(object, ...))
#' @rdname stats
setMethod("residuals", signature(object="fitrad"),
          function(object, ...) object@rad.tab$abund - fitted(object, ...))
