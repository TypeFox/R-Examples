

### longitudinal.R  (2005-05-06)
###
###    Data structure for longitudinal data
###
### Copyright 2005 Korbinian Strimmer
###
###
###
### This file is part of the `GeneTS' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA





######### PUBLIC ######## 

## basic functions

is.longitudinal = function(x)
{
  if(inherits(x, "longitudinal")) 
    return(TRUE)
  else
    return(FALSE)
}


as.longitudinal = function(x, repeats=1, time)
{ 
   if (!is.matrix(x) )
    stop("only matrices can be coerced to longitudinal")
 
  
  # create repeats vector
  if (length(repeats) == 1) 
  {
    if (dim(x)[1] %% repeats != 0)
      stop("number of repeats incompatible with number of rows of data matrix")
    repeats = rep(repeats, dim(x)[1] %/% repeats) 
  } 
  if (sum(repeats) != dim(x)[1])
    stop("sum of repeats must equal the number of rows in data matrix")

  # create time vector  
  if (missing(time)) time = 1:length(repeats)
  if (any(duplicated(time))) 
    stop("duplicated entries in time vector") 
  if (length(time) != length(repeats))
    stop("length of time vector must be equal to the length of repeats vector")
  if (any( diff(time) <= 0 ))  
    stop("entries time vector must be monotonically increasing")
  
 
  # construct longitudinal object
  attr(x, "class") = "longitudinal"
  attr(x, "time") = time
  attr(x, "repeats") = repeats
  
  rn = NULL
  for (i in 1:length(repeats))
  {
    rn = c(rn, paste( time[i], seq(1, repeats[i]), sep="-") )
  }
  rownames(x) = rn
   
  return(x)
}




## basic output stuff

summary.longitudinal = function(object, ...)
{
  x = object
  
  if (has.repeated.measurements(x))
     tmp = "yes"
  else
      tmp = "none"
  
  cat("Longitudinal data:\n ")
  cat(paste(dim(x)[2], "variables measured at",
   length(get.time.repeats(x)$time), "different time points\n"))
  
  cat(paste(" Total number of measurements per variable:", dim(x)[1]), "\n")
  
  cat(paste(" Repeated measurements:", tmp, "\n"))
  
  cat("\n")
  cat(" To obtain the measurement design call 'get.time.repeats()'.\n")

}

print.longitudinal = function(x, ...)
{
  summary(x)
  cat("\n")
   
  attr(x, "class") = attr(x, "time") = attr(x, "repeats") = NULL

  NextMethod("print", x, quote = FALSE, right = TRUE)
}

plot.longitudinal = function(x, series=1, type=c("median", "mean"), autolayout=TRUE, ...)
{
   type=match.arg(type)
   if (type=="median") func=median
   if (type=="mean") func=mean

   
   cx = condense.longitudinal(x, series, func)
   lc = length(series)
   tr = get.time.repeats(x)
   
   nrow = ceiling(sqrt(lc))
   ncol = ceiling(lc/nrow)
   
   if(autolayout) par(mfrow=c(nrow,ncol))
   
   xx = rep(tr$time, tr$repeats)
   for (i in 1:lc)
   {
     name = colnames(x)[series[i]]
     if(is.null(name))
       name = paste("Series", series[i])
     plot(xx, x[,series[i]], main=name, ylab="value", xlab="time", col=gray(0.7), ...)
     lines(tr$time, cx[,i], col=4)
   }
   if(autolayout) par(mfrow=c(1,1))
}



