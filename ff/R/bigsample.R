# Sampling from big pools for ff
# (c) 2007 Daniel Adler & Jens Oehlschlägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2007-08-24
# Last changed: 2007-10-25

# source("d:/mwp/eanalysis/ff/R/bigsample.R")

#! \name{bigsample}
#! \alias{bigsample}
#! \alias{bigsample.default}
#! \alias{bigsample.ff}
#! \title{ Sampling from large pools }
#! \description{
#!    \command{bigsample} samples quicker from large pools than \command{\link{sample}} does.
#! }
#! \usage{
#! bigsample(x, ...)
#! \method{bigsample}{default}(x, size, replace = FALSE, prob = NULL, negative = FALSE, ...)
#! \method{bigsample}{ff}(x, size, replace = FALSE, prob = NULL, ...)
#! }
#! \arguments{
#!   \item{x}{ the pool to sample from }
#!   \item{size}{ the number of elements to sample }
#!   \item{replace}{ TRUE to use sampling with replacement }
#!   \item{prob}{ optional vector of sampling probabilities (recyled to pool length) }
#!   \item{negative}{ \code{negative} }
#!   \item{\dots}{ \code{\dots} }
#! }
#! \details{
#!    For small pools \command{\link{sample}} is called.
#! }
#! \note{
#!   Note that \command{bigsample} and \command{sample} do not necessarily return the same sequence of elements when \command{set.seed} is set before.
#! }
#! \value{
#!   a vector of elements sampled from the pool (argument 'x')
#! }
#! \author{ Daniel Adler, Jens Oehlschlägel, Walter Zucchini}
#! \seealso{ \code{\link{sample}}, \code{\link{ff}} }
#! \examples{
#! message("Specify pool size")
#! bigsample(1e8, 10)
#! message("Sample ff elements (same as x[bigsample(length(ff(1:100 / 10)), 10)])")
#! bigsample(ff(1:100 / 10), 10)
#!  \dontrun{
#!    message("Speed factor")
#!      (system.time(for(i in 1:10)sample(1e8, 10))[3]/10) 
#!    / (system.time(for(i in 1:1000)bigsample(1e8, 10))[3]/1000)
#!  }
#! }
#! \keyword{ distribution }
#! \keyword{ data }

bigsample.default <-
function (x, size, replace = FALSE, prob = NULL, negative=FALSE
, ... # dummy to keep R CMD check quiet
){
  # delegate to sample if x is not a total length
  if (length(x)!=1)
    return(sample(x, size, replace, prob))

  total <- as.integer(x)
  if (missing(size))
    size <- total
  else
    size <- as.integer(size)

  if (!replace && negative && size>total/2){
    neg <- TRUE
    size <- total - size
  }else{
    neg <- FALSE
  }

  # delegate to sample if prob needed or size is >1% but too large because then sample is faster
  if (!is.null(prob) || (size>total/100 && size < 1e5 )){
    if (neg){
      ret <- -sample(total, size, replace, prob)
      setattr(ret ,"maxindex", total)  # attr(ret, "maxindex") <- total
      return(ret)
    }else{
      return(sample(total, size, replace, prob))
    }
  }

  # only here do something different from sample
  stopifnot(replace || size <= total)
  if (replace){
    if (neg){
      ret <- -as.integer(ceiling(runif(size, max=total)))
      setattr(ret ,"maxindex", total)  # attr(ret, "maxindex") <- total
      return(ret)
    }else{
      return(as.integer(ceiling(runif(size, max=total))))
    }
  }else{
    n <- size
    x <- integer()
    while (TRUE) {
      x <- unique(c(x, as.integer(ceiling(runif(n, max=total)))))
      n <- size - length(x)
      if (length(x) == size){
        if (neg){
          x <- -x
          setattr(ret ,"maxindex", total)  # attr(ret, "maxindex") <- total
          return(x)
        }else{
          return(x)
        }
      }
    }
  }
}

# removed sorting here,
# sorting belongs to general optimization for contiguous reads
# sorting is encapsuled in class hi
bigsample.ff <-
function (x, size, replace = FALSE, prob = NULL
, ... # dummy to keep R CMD check quiet
){
    x[bigsample.default(length(x), size, replace, prob)]
}
