##' closest
##' 
##' see function
##' 
##' 
##' @keywords internal
##' @export closest
"closest" <- function(vec, val)
{
  ## return the index of the value in vec which is closest to val
  newval <- min(abs((vec - val)))
  z <- abs(vec - val)
  temp <- z == newval
  nums <- c(1:length(vec))
  nums[temp]
}











##' num label
##' 
##' see function
##' 
##' 
##' @keywords internal
##' @export label_num
"label_num" <- function(labs)
{
  ## labs: a vector of labels
  ## convert labels to integers
  qq <- unique(labs)
  nums <- c(1:length(qq))
  for(j in 1:length(qq)) {
    temp <- labs == qq[j]
    labs[temp] <- nums[j]
  }
  as.numeric(labs)
}










##' convert label
##' 
##' see function
##' 
##' 
##' @keywords internal
##' @export label_convert
"label_convert" <- function(segs.or.labels, match, replace)
{
  ## label_convert --
  ##  map one set of labels to another in a segment
  ##  list or label vector
  ##
  if (is.seglist( segs.or.labels ) ) {
    labs <- label(segs.or.labels)
  } else {
    labs <- segs.or.labels
  }
  if( length(match) != length(replace) ) {
    ## this is only ok if length(replace) == 1 so that we
    ## replace anything in match with replace
    if( length(replace) != 1 ) {
      stop("Lengths of match and replace vectors differ in label_convert")
    } 
  }
  
  if (length(replace) == 1) {
    temp <- muclass(labs, match)
    labs[temp] <- replace
  } else {
    for( i in 1:length(match) ) {
      labs[labs==match[i]] <- replace[i]
    }
  }
  if (is.seglist(segs.or.labels)) {
    return( modify.seglist(segs.or.labels, labels=labs))
  } else {
    return( labs )
  }
}










##' Read matrix data from a file
##' 
##' Reads matrix data from a file
##' 
##' This function has been partially superceeded by the introduction of data
##' frames and the read.table function.  It is still useful however for reading
##' data into Splus matrix objects.
##' 
##' @param file A filename.
##' @param num.cols The number of columns of data in the file.
##' @param what A template for the data elements in the file, it should be a
##' number for numeric data (the default) or a string for string data. Note
##' that an Splus matrix can only hold one type of data (string or numeric),
##' for mixed types use data tables and the \code{read.table} function.
##' @param sk The number of leading lines of the file to skip.
##' @return A matrix corresponding to the data in \code{file}.
##' @seealso read.table
##' @keywords misc
##' @export matscan
"matscan"<- function(file, num.cols=utils::count.fields(file)[1], what = 0, sk = 0)
{
  ## first make a template, a list of num.cols what's
  template <- as.list(rep(what, num.cols))
  data <- scan(file, template, skip=sk,quote = "")
  num.rows <- length(data[[1]])
  mat <- matrix(0, num.rows, num.cols)
  for(i in 1:num.cols) {
    mat[, i] <- data[[i]]
  }
  mat
}











##' Find common elements in vectors
##' 
##' Finds common elements in vectors
##' 
##' 
##' @param labels A vector of labels.
##' @param class A label or vector of labels.
##' @return A logical vector which is T for each element in \code{labels} which
##' matches \code{class} or an element of \code{class}.
##' @seealso match
##' @keywords misc
##' @examples
##' 
##' muclass(c("a", "b", "c"), c("a", "c"))
##' #[1] T F T
##' 
##' @export muclass
"muclass"<- function(labels, class)
{
  !(is.na(match(labels, class)))
}









##' Randomise or Reverse items in a segment list
##' 
##' Randomises or Reverses items in a segment list
##' 
##' 
##' @param segs An Emu segment list.
##' @param bwd If T, reverse the order of the segment list.
##' @param rand If T, randomise the order of the segment lists (default).
##' @return A segment list containing the original elements in random or
##' reversed order. This is useful if the segment list is to be used as the
##' source for a set of stimuli in a perception experiment.
##' @seealso \code{\link{query}}
##' @keywords misc
##' @examples
##' 
##' data(vowlax)
##' ## assumes a database called demo is available on your system and that
##' ## the Emu system is installed. 
##' 
##' # all Phonetic vowels in the database
##' segs <- vowlax
##' 
##' # randomise the segment list
##' rsegs <- randomise.segs( segs )
##' 
##' 
##' @export randomise.segs
"randomise.segs" <-  function( segs, rand = TRUE, bwd=FALSE )
{
  if( bwd ){
    ## reverse the segment list
    segs[nrow(segs):1,]
  } else {
    ## randomise the segment list
    segs[sample(1:nrow(segs)),]
  }
}









##' Converts degrees to radians
##' 
##' Converts degrees to radians
##' 
##' There are 360 degrees or 2 * PI radians in one full rotation.
##' 
##' @param degrees Angular measurement for conversion.
##' @return Angular measurement in radians.
##' @keywords misc
##' @export radians
"radians"<- function(degrees) (degrees * 2 * pi)/360









##' Sort matrix by label
##' 
##' Sorts matrix by label
##' 
##' 
##' @param mat A mu+ sement matrix.
##' @param labs A label vector which has the same number of columns as
##' \code{mat}.
##' @return Returns a sorted matrix by label, created from \code{mat}.
##' @seealso label, phon
##' @keywords misc
##' @export sortmatrix
"sortmatrix" <- function(mat, labs = dimnames(mat)[[2]])
{
  ## labs can also be a vector of labels, which has
  ## the same number of columns as mat; e.g.
  ## if labs is "p" "t" "k", then mat
  ## will be sorted with these three labels in the first
  ## three columns
  b1 <- labs
  b2 <- dimnames(mat)[[1]]
  c1 <- match(b2, b1)
  d1 <- cbind(c1, mat)
  newmat <- d1[sort.list(d1[, 1]),  ]
  newmat <- newmat[, -1]
  b1 <- dimnames(newmat)[[1]]
  b2 <- dimnames(newmat)[[2]]
  c1 <- match(b2, b1)
  d1 <- rbind(c1, newmat)
  newmat2 <- d1[, (sort.list(d1[1,  ],  ))]
  newmat2[-1,  ]
}











##' Function to convert between Hertz and Radians
##' 
##' convert between Hertz and Radians
##' 
##' 
##' @param vec A numerical vector of frequencies in Hz or radians
##' @param samfreq A single element numerical vector of the sampling frequency.
##' Defaults to 20000 Hz
##' @param hz Logical. If T, convert from Hz to radians otherwise from radians
##' to hz
##' @author Jonahtan Harrington
##' @seealso \code{\link{help}}
##' @keywords math
##' @examples
##' 
##' # 4000 Hz in radians at a sampling frequency of 8000 Hz
##' rad(4000, 8000)
##' # pi/2 and pi/4  radians in Hz at a sampling frequency of 10000 Hz
##' rad(c(pi/2, pi/4), 10000, FALSE)
##' 
##' 
##' @export rad
"rad" <- function(vec, samfreq = 20000, hz = TRUE)
{
  # hz: if T, vec is a vector in Hertz, otherwise it's radians
  # convert from radians to Hz, or Hz to radians
  if(hz) vals <- (vec * 2 * pi)/samfreq else vals <- (vec * samfreq)/(2 * 
                                                                        pi)
  vals
}









##' Function to find the column number corresponding to frequencies of a
##' spectral object
##' 
##' Find the column number corresponding to frequencies of a spectral object.
##' 
##' This function is used in conjunction with object oriented programming of
##' EMU spectral objects. It should not in general be called from inside a
##' function. Its principal use is to determine the column number(s)
##' corresponding to frequencies for spectral trackdata objects or spectral
##' matrices or the element number for spectral vectors.
##' 
##' @param trackdata A spectral object
##' @param j A vector of frequencies
##' @author Jonathan Harrington
##' @keywords math
##' @examples
##' 
##' freqtoint(fric.dft,1000:2000)
##' # all frequencies except 1000-2000
##' freqtoint(vowlax.dft.5, -(1000:2000))
##' # all frequencies except 1000 Hz
##' freqtoint(e.dft, -1000)
##' # the d.c. offset - i.e. column 1
##' freqtoint(vowlax.dft.5, 0)
##' # all freqs except the d.c. offset - i.e. not column 1
##' freqtoint(vowlax.dft.5, -1)
##' 
##' 
##' 
##' 
##' @export freqtoint
"freqtoint" <- function(trackdata, j){
  # note to remove the dc offset, set j to -1
  sg <- sign(j)
  zerowhich <- sg==0
  if(all(sg[!zerowhich] > 0))
    sg[zerowhich] = 1
  else sg[zerowhich] = -1
  j <- abs(j)
  fs <- trackfreq(trackdata)
  N <- length(fs)
  res <- 1+ (((j - fs[1]) * (N-1))/(fs[N] - fs[1]))
  res[res < 1] = 1
  res[res > N] = N
  res <- sg * res
  unique(round(res))
}










##' Function to dB-normalise spectral objects
##' 
##' The function can be used to rescale a spectrum to a dB value at a
##' particular frequency - for example, to rescale the spectrum so that 3000 Hz
##' has 0 dB and all other values are shifted in relation to this.
##' 
##' 
##' @param specdata An object of class 'spectral'
##' @param f A single element vector specifying the frequency. Defaults to 0
##' @param db A single element vector specifying the dB value to which the
##' spectrum is to be rescaled. Defaults to zero
##' @return An object of the same class with rescaled dB values. The default is
##' to rescale the dB-values of the spectrum to 0 dB at 0 Hz.
##' @author Jonathan Harrington
##' @seealso \code{\link{dbtopower}} \code{\link{plot.spectral}}
##' @keywords manip
##' @examples
##' 
##' # normalise to - 40 dB at 1500 Hz
##' res = dbnorm(e.dft, 1500, 0)
##' # compare the two
##' ylim = range(c(res, e.dft))
##' plot(e.dft, ylim=ylim, type="l")
##' par(new=TRUE)
##' plot(res, ylim=ylim, type="l", col=2)
##' 
##' 
##' @export dbnorm
"dbnorm" <- function(specdata, f=0, db=0)
{
  
  
  if(is.trackdata(specdata))
    dat <- specdata$data
  else
    dat <- specdata
  
  
  # normalise to dbnorm
  minfun <- function(specvals, f, db)
  {
    specvals <- specvals - specvals[f]+db
  }
  if(is.matrix(dat))
    dat <- fapply(dat, minfun, f, db)
  else
    dat = dat-dat[f]+db
  
  if(is.trackdata(specdata))
  {
    specdata$data <- dat
    return(specdata)
  }
  else return(dat)
  
}











##' Function for inter-converting between decibels and a linear scale
##' 
##' The function converts from decibels to a linear scale
##' 
##' The function returns base\eqn{\mbox{\textasciicircum}}{^}(specdata/const)
##' if inv=F, otherwise, const * log(dat, base=base). If the object to which
##' this function is applied is of class 'trackdata' then this function is
##' applied to \$data.
##' 
##' @param specdata A numeric object or an object of class trackdata
##' @param const A single element numeric vector. Defaults to 10
##' @param base A single element numeric vector. Defaults to 10
##' @param inv Logical. If T, then the conversion is from a logarithmic to an
##' anti-logarithmic form, otherwise the other way round
##' @return An object of the same class.
##' @author Jonathan Harrington
##' @seealso \code{\link{dbtopower}} \code{\link{plot.spectral}}
##' @keywords math
##' @examples
##' 
##' 
##' # convert 10 dB to a power ratio
##' vec = dbtopower(10)
##' # convert dB-data to a power ratio and back to decibels
##' res = dbtopower(vowlax.dft.5)
##' res = dbtopower(res, inv=TRUE)
##' 
##' @export dbtopower
"dbtopower" <- function(specdata, const = 10, base=10, inv=FALSE)
{
  # function for converting from db to power and back
  if(is.trackdata(specdata))
    dat <- specdata$data
  else dat <- specdata
  if(!inv)
    result <- base^(dat/const)
  else
    result <- const * log(dat, base=base)
  if(is.trackdata(specdata))
  {
    specdata$data <- result
    return(specdata)
  }
  else return(result)
}










##' Function to shift the elements of a vector.
##' 
##' The function makes use of the function 'fitler' to delay or advance a
##' signal by k points.
##' 
##' The function makes use of the function 'filter' for linear filtering to
##' carry out the shifting.
##' 
##' @param x A numeric vector
##' @param delta A single element numeric vector. Defines the number of points
##' by which the signal should be shifted.
##' @param circular Logical. If T, the signal is wrapped around itself so that
##' if delta = 1, x[n] becomes x[1]. Otherwise, if delta is positive, the same
##' number of zeros are prepended to the signal
##' @return The signal shifted by a certain number of points.  ...
##' @author Jonathan Harrington
##' @seealso filter
##' @keywords manip
##' @examples
##' 
##' vec = 1:10
##' shift(vec, 2)
##' shift(vec, -2)
##' shift(vec, 2, circular=FALSE)
##' 
##' 
##' 
##' @export shift
"shift" <- function(x, delta = 1,  circular = TRUE)
{
  ## converts x[n] into x[n-1] by multiplying the Fourier
  ## transform on x[n] by z^-1 i.e. by e^-iw
  N <- length(x)
  if(delta < 0)
    delta <- N + delta
  h <- c(rep(0, delta), 1)
  if(!circular) {
    N <- length(x) + length(h) - 1
    x <- c(x, rep(0, N - length(x)))
  }
  filter(x, h, sides = 1, circular = TRUE)[1:N]
  
}













##' Split a string into words.
##' 
##' Splits a string into words.
##' 
##' 
##' @param str A string.
##' @param char A character to split on
##' @return A vector of strings.  The original \code{str} is split at ever
##' occurrence of \code{char} to generate a vector of strings.
##' @keywords misc
##' @examples
##' 
##' splitstring("/home/recog/steve/foo", "/")
##' #[1] "home" "recog" "steve" "foo"
##' 
##' @export splitstring
splitstring <- function(str,char) {
  if(str == "")
    mat <- c(str)
  else {
    mat <- NULL
    ind <- 1
    cont <- TRUE
    while(TRUE) {
      ministr <- NULL
      length <- 0
      while(TRUE) {
        ch <- substring(str, ind, ind)
        if(ch == char) {
          ind <- ind + 1
          break
        }
        if(ch == "") {
          break
        }
        ministr <- c(ministr, ch)
        ind <- ind + 1
        length <- length + 1
      }
      ## now concatenate string
      if(length > 0)
        mat <- c(mat, paste(ministr, collapse = ""))
      if(ch == "")
        break
    }
  }
  mat
}
