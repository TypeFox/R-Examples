##' Function to test whether the object is of class "spectral"
##' 
##' Returns T or F depending on whether the object is of class "spectral"
##' 
##' 
##' @param dat An R object
##' @return A single element logical vector: T or F
##' @author Jonathan Harrington
##' @seealso \code{\link{as.spectral}}
##' @keywords attribute
##' @examples
##' 
##' 
##' is.spectral(vowlax.dft.5)
##' is.spectral(fric.dft)
##' is.spectral(fric.dft$data)
##' is.spectral(vowlax.dft.5[1,])
##' is.spectral(fric.dft[1,1])
##' 
##' 
##' 
##' @export is.spectral
"is.spectral" <- function(dat)
{
  if(!is.trackdata(dat))
    return(any(class(dat) %in% "spectral"))
  else
    return(any(class(dat$data) %in% "spectral"))
}











##' Function to convert an object into an object of class 'spectral'.
##' 
##' The function converts a vector, matrix, or EMU-trackdata object into an
##' object of the same class and of class 'spectral'
##' 
##' If fs is a single element numeric vector, then the frequencies of trackdata
##' are defined to extend to fs/2. If fs is missing, then the frequencies are
##' 0:(N-1) where N is the length of trackdata.
##' 
##' @param trackdata A vector, matrix, or EMU-trackdata object.
##' @param fs Either a single element numeric vector, or a numeric vector of
##' the same length as the length of trackdata if trackdata is a vector, or of
##' the same number of rows as trackdata
##' @return The same object but of class 'spectral'.
##' @author Jonathan Harrington
##' @seealso \code{\link{is.spectral}} \code{\link{plot.spectral}}
##' @keywords attribute
##' @examples
##' 
##' vec = 1:10
##' as.spectral(vec, 2000)
##' mat = rbind(1:10, 1:10)
##' as.spectral(mat)
##' # turn a spectral trackdata object into a trackdata object
##' tr = as.trackdata(rbind(fric.dft$data), fric.dft$index, fric.dft$ftime)
##' # turn it into a spectral trackdata object with sampling freq 16 kHz
##' tr = as.spectral(tr, 16000)
##' # list the frequencies
##' trackfreq(tr)
##' # Notice that only the $data is made into a spectral matrix,
##' # not the entire trackdata object
##' # so this is trackdata
##' class(tr)
##' # this is a spectral matrix
##' class(tr$data)
##' 
##' 
##' 
##' 
##' @export as.spectral
"as.spectral" <- function(trackdata, fs)
{
  if(is.trackdata(trackdata)){
    
    if(is.spectral(trackdata$data)) {
      warning("matrix is already of class spectral")
      return(trackdata)
    }
    N <- ncol(trackdata$data)
    if(missing(fs))
      fs <- 0: (ncol(trackdata$data)-1)
    else{
      if(length(fs)==1)
      {
        fs <- fs/2
        fs <- seq(0, fs, length=N)
      }
    }
    attr(trackdata$data, "fs") <- fs
    class(trackdata$data) <- c(class(trackdata$data), "spectral")
  }
  
  else if (is.matrix(trackdata)){
    if(is.spectral(trackdata)) {
      warning("matrix is already of class spectral")
      return(trackdata)
    }
    N <- ncol(trackdata)
    if(missing(fs))
      fs <- 0: (ncol(trackdata)-1)
    else{
      if(length(fs)==1)
      {
        fs <- fs/2
        fs <- seq(0, fs, length=N)
      }
    }
    attr(trackdata, "fs") <- fs
    class(trackdata) <- c(class(trackdata), "spectral")
  }
  else
  {
    
    if(is.spectral(trackdata)){
      warning("matrix is already of class spectral")
      return(trackdata)
    }
    N <- length(trackdata)
    if(missing(fs))
      fs <- 0: (length(trackdata)-1)
    else{
      if(length(fs)==1)
      {
        fs <- fs/2
        fs <- seq(0, fs, length=N)
      }
    }
    attr(trackdata, "fs") <- fs
    class(trackdata) <- c(class(trackdata), "spectral")
  }
  trackdata
}












##' Plot spectra from EMU spectral objects
##' 
##' The function plots spectrum of any EMU spectral object.
##' 
##' This function is implemented when a spectral trackdata object is called
##' with the 'plot' function.
##' 
##' @param x An EMU object of class 'spectral'
##' @param labs An optional vector character labels. Must be the same length as
##' specdata
##' @param ylim A two-element numeric vector for the y-axis range (see 'par')
##' @param xlim A two-element numeric vector for the x-axis range (see 'par')
##' @param col Specify a color - see 'mu.colour')
##' @param lty Specify a linetype - see 'mu.colour'
##' @param lwd Specify line thickness - see 'mu.colour'
##' @param fun An R function name e.g., mean, var, sum, etc. The function is
##' applied separately to each category type specified in labs
##' @param freq A numeric vector the same length as the number of columns in
##' specdata specifying the frequencies at which the spectral data is to be
##' plotted. If not supplied, defaults to trackfreq(specdata)
##' @param type A single element character vector for the linetype
##' @param power Logical. If T, then specdata (or specdata\$data if specdata is
##' a trackdata object, is converted to a *
##' specdata\eqn{\mbox{\textasciicircum}}{^}b, where a and b have the values
##' given in powcoeffs. This operation is applied before b
##' @param powcoeffs A two-element numeric vector. Defaults to c(10, 10)
##' @param dbnorm Logical. If T, apply dB-level normalization per spectrum as
##' defined by dbcoeffs below. Defaults to F.
##' @param dbcoeffs A two element numeric vector (x, y). The spectra are
##' normalised in such a way that the values of each spectrum at a frequency of
##' y are set to a dB level of x. For example, to normalise the spectrum to 10
##' dB at 2000 Hz, set dbnorm to T and dbcoeffs to c(2000, 10)
##' @param legend Parameters for defining the legend. See 'mu.legend' for
##' further details
##' @param axes A logical vector indicating whether the axes should be plotted
##' @param \dots Further graphical parameters may be supplied.
##' @note To plot spectral data from a spectral trackdata object, then call the
##' function explicitly with 'plot/spectral' rather than with just 'plot'
##' @author Jonathan Harrington
##' @seealso \code{\link{plot}} \code{\link{plot.trackdata}}
##' \code{\link{as.spectral}}
##' @keywords dplot
##' @examples
##' 
##' plot(vowlax.dft.5[1,])
##' 
##' # with label types
##' plot(vowlax.dft.5[1:20,], vowlax.l[1:20])
##' 
##' # As above but averaged after converting to power ratios.
##' plot(vowlax.dft.5[1:20,], vowlax.l[1:20], fun=mean, power=TRUE)
##' 
##' # All the spectra of one segment in a trackdata object
##' plot(fric.dft[1,])
##' 
##' 
##' 
##' 
##' @export
"plot.spectral" <- function (x, labs, ylim, xlim,  col, lty, 
                             lwd, fun, freq, type = "l", 
                             power = FALSE, powcoeffs = c(10, 10), 
                             dbnorm = FALSE, dbcoeffs = c(0, 0), 
                             legend = TRUE, axes=TRUE,  ...) 
{
  specdata = x
  if (is.trackdata(specdata)) 
    specdata <- specdata$data
  if (!is.spectral(specdata)) 
    stop("specdata must be of class spectral")
  if (dbnorm) 
    specdata <- dbnorm(specdata, dbcoeffs[1], dbcoeffs[2])
  if (missing(freq)) 
    f <- trackfreq(specdata)
  else f <- freq
  if (is.matrix(specdata)) 
    N <- nrow(specdata)
  else {
    N <- 1
    specdata <- rbind(specdata)
  }
  if (missing(labs)) 
    labs <- rep(".", N)
  if (!missing(fun)) {
    if (power) 
      specdata <- dbtopower(specdata, powcoeffs[1], powcoeffs[2])
    mat <- list(NULL)
    for (j in unique(labs)) {
      temp <- labs == j
      v <- apply(specdata[temp, ], 2, fun)
      mat$fn <- rbind(mat$fn, v)
      mat$l <- c(mat$l, j)
    }
    dimnames(mat$fn) <- list(mat$l, dimnames(specdata)[[2]])
    specdata <- mat$fn
    if (power) 
      specdata <- dbtopower(specdata, powcoeffs[1], powcoeffs[2], 
                            inv = TRUE)
    if (length(unique(labs)) > 1) 
      labs <- dimnames(specdata)[[1]]
    else {
      labs <- unique(labs)
      specdata <- rbind(specdata)
    }
  }
  if (missing(ylim)) 
    ylim <- range(specdata)
  if (missing(xlim)) 
    xlim <- range(f)
  if (missing(col)) 
    col <- TRUE
  if (missing(lty)) 
    lty <- FALSE
  if (missing(lwd)) 
    lwd <- NULL
  cols <- mu.colour(labs, col, lty, lwd)
  for (j in 1:nrow(specdata)) {
    graphics::plot(f, specdata[j, ], type = type, col = cols$colour[j], 
         lty = cols$linetype[j], lwd = cols$lwd[j], xlim = xlim, 
         ylim = ylim, xlab = "", ylab = "", axes = FALSE)
    graphics::par(new = TRUE)
  }
  if (is.logical(legend)) {
    if (legend & length(unique(labs)) > 1) {
      legend <- "topright"
      legend(legend, NULL, cols$legend$lab, col = cols$legend$col, 
             lty = as.numeric(cols$legend$lty), lwd = as.numeric(cols$legend$lwd))
    }
  }
  else legend(legend, NULL, cols$legend$lab, col = cols$legend$col, 
              lty = as.numeric(cols$legend$lty), lwd = as.numeric(cols$legend$lwd))
  if(axes)
  {
    graphics::axis(side = 1)
    graphics::axis(side = 2)
  }
  graphics::title(...)
  graphics::box(...)
}


##' @export
"bark.spectral" <- function (f, ...) 
{
  specobject = f
  if (!is.trackdata(specobject)) {
    if (!is.matrix(specobject)) 
      specobject <- as.spectral(rbind(specobject), attr(specobject, 
                                                        "fs"))
  }
  f <- trackfreq(specobject)
  b <- bark(f)
  temp <- b < 0
  if (any(temp)) 
    specobject <- specobject[, !temp]
  f <- trackfreq(specobject)
  b <- bark(f)
  N <- length(b)
  ba <- seq(min(b), max(b), length = N)
  if (is.trackdata(specobject)) 
    spec <- specobject$data
  else if (is.matrix(specobject)) 
    spec <- specobject
  else spec <- as.spectral(rbind(specobject), attr(specobject,"fs"))
  res <- NULL
  for (j in 1:nrow(spec)) {
    v = approx(b, c(spec[j, ]), ba)
    res <- rbind(res, v$y)
  }
  if (is.trackdata(specobject)) {
    specobject$data <- res
    if (!is.null(tracktimes(spec))) 
      rownames(specobject$data) <- tracktimes(spec)
    specobject <- as.spectral(specobject, ba)
  }
  else {
    specobject <- res
    specobject <- as.spectral(specobject, ba)
  }
  specobject
}

##' @export
"mel.spectral" <- function (a) 
{
  specobject = a
  if (!is.trackdata(specobject)) {
    if (!is.matrix(specobject)) 
      specobject <- as.spectral(rbind(specobject), attr(specobject, "fs"))
  }
  f <- trackfreq(specobject)
  b <- mel(f)
  N <- length(b)
  ba <- seq(min(b), max(b), length = N)
  if (is.trackdata(specobject)) 
    spec <- specobject$data
  else if (is.matrix(specobject)) 
    spec <- specobject
  else spec <- as.spectral(rbind(specobject), attr(specobject, 
                                                   "fs"))
  res <- NULL
  for (j in 1:nrow(spec)) {
    v = approx(b, c(spec[j, ]), ba)
    res <- rbind(res, v$y)
  }
  if (is.trackdata(specobject)) {
    specobject$data <- res
    if (!is.null(tracktimes(spec))) 
      rownames(specobject$data) <- tracktimes(spec)
    specobject <- as.spectral(specobject, ba)
  }
  else {
    specobject <- res
    specobject <- as.spectral(specobject, ba)
  }
  specobject
}
