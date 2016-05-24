#' Functions to generate violin scatter plots
#'
#' Arranges data points using quasirandom noise (van der Corput sequence) to create a plot resembling a cross between a violin plot (showing the density distribution) and a scatter plot (showing the individual points). The development version of this package is on \url{http://github.com/sherrillmix/vipor}
#'
#' The main functions are: 
#'      \describe{
#'        \item{\code{\link{offsetX}}:}{calculate offsets in X position for plotting (groups of) one dimensional data}
#'        \item{\code{\link{vpPlot}}:}{a simple wrapper around plot and offsetX to generate plots of grouped data}
#'      }
#'
#' @docType package
#' @name vipor
#' @author Scott Sherrill-Mix, \email{shescott@@upenn.edu}
#' @seealso \url{http://github.com/sherrillmix/vipor}
#' @examples
#' dat<-list(rnorm(100),rnorm(50,1,2))
#' ids<-rep(1:length(dat),sapply(dat,length))
#' offset<-offsetX(unlist(dat),ids)
#' plot(unlist(dat),ids+offset)
NULL

#' Plot data using offsets by quasirandom noise to generate a violin point plot
#' 
#' Arranges data points using quasirandom noise (van der Corput sequence), pseudorandom noise or alternatively positioning extreme values within a band to the left and right to form beeswarm/one-dimensional scatter/strip chart style plots. That is a plot resembling a cross between a violin plot (showing the density distribution) and a scatter plot (showing the individual points) and so here we'll call it a violin point plot.
#' 
#' @param x a grouping factor for y (optional)
#' @param y vector of data points
#' @param xaxt if 'n' then no x axis is plotted
#' @param offsetXArgs a list with arguments for offsetX
#' @param ... additional arguments to plot
#' @return invisibly return the adjusted x positions of the points
#' @export
#' @seealso \code{\link{offsetX}}
#' @examples
#' dat<-list(
#'   'Mean=0'=rnorm(200),
#'   'Mean=1'=rnorm(50,1),
#'   'Bimodal'=c(rnorm(40,-2),rnorm(60,2)),
#'   'Gamma'=rgamma(50,1)
#' )
#' labs<-factor(rep(names(dat),sapply(dat,length)),levels=names(dat))
#' vpPlot(labs,unlist(dat))
vpPlot<-function(x=rep('Data',length(y)),y,xaxt='y',offsetXArgs=NULL,...){
  x<-as.factor(x)
  ids<-as.numeric(x)
  labels<-levels(x)
  labelIds<-1:length(labels)
  names(labelIds)<-labels
  xPos<-ids+do.call(offsetX,c(list(y),list(x),offsetXArgs))
  graphics::plot(xPos,y,...,xaxt='n',xlab='')
  if(xaxt!='n')graphics::axis(1,labelIds,labels,...)
  return(invisible(xPos))
}

#' Offset data using quasirandom noise to avoid overplotting
#' 
#' Arranges data points using quasirandom noise (van der Corput sequence), pseudorandom noise or alternatively positioning extreme values within a band to the left and right to form beeswarm/one-dimensional scatter/strip chart style plots. That is a plot resembling a cross between a violin plot (showing the density distribution) and a scatter plot (showing the individual points). This function returns a vector of the offsets to be used in plotting.
#' 
#' @param y vector of data points
#' @param x a grouping factor for y (optional)
#' @param width the maximum spacing away from center for each group of points. Since points are spaced to left and right, the maximum width of the cluster will be approximately width*2 (0 = no offset, default = 0.4)
#' @param varwidth adjust the width of each group based on the number of points in the group
#' @param ... additional arguments to offsetSingleGroup
#' @return a vector with of x-offsets of the same length as y
#' @export
#' @examples
#' ## Generate fake data
#' dat <- list(rnorm(50), rnorm(500), c(rnorm(100), rnorm(100,5)), rcauchy(100))
#' names(dat) <- c("Normal", "Dense Normal", "Bimodal", "Extremes")
#' 
#' ## Plot each distribution with a variety of parameters
#' par(mfrow=c(4,1), mar=c(2,4, 0.5, 0.5))
#' sapply(names(dat),function(label) {
#'   y<-dat[[label]]
#'   
#'   offsets <- list(
#'     'Default'=offsetX(y),
#'     'Smoother'=offsetX(y, adjust=2),
#'     'Tighter'=offsetX(y, adjust=0.1),
#'     'Thinner'=offsetX(y, width=0.1)
#'   )
#'   ids <- rep(1:length(offsets), sapply(offsets,length))
#'   
#'   plot(unlist(offsets) + ids, rep(y, length(offsets)), 
#'        ylab=label, xlab='', xaxt='n', pch=21, las=1)
#'   axis(1, 1:4, c("Default", "Adjust=2", "Adjust=0.1", "Width=10%"))
#' })
#' 
offsetX <- function(y, x=rep(1, length(y)), width=0.4, varwidth=FALSE,...) {
  if (length(x)!=length(y)) stop(simpleError('x and y not the same length in offsetX'))

  maxLength<-max(table(x))

  # Apply the van der Corput noise to each x group to create offsets
  new_x <- aveWithArgs(y,x, FUN=offsetSingleGroup,maxLength=if(varwidth){maxLength}else{NULL},...)
  new_x <- new_x*width

  return(new_x)
}




# Offset data to avoid overplotting for a single subgroup of data
# 
# Arranges data points using quasirandom noise (van der Corput sequence), pseudorandom noise or alternatively positioning extreme values within a band to the left and right to form beeswarm/one-dimensional scatter/strip chart style plots. Returns a vector of the offsets to be used in plotting. This function is mostly used as a subroutine of \code{\link{offsetX}}
# @param y y values for a single group for which offsets should be calculated
#' @param maxLength multiply the offset by sqrt(length(y)/maxLength) if not NULL. The sqrt is to match boxplot (allows comparison of order of magnitude different ns, scale with standard error)
#' @param method method used to distribute the points
#' @param nbins the number of points used to calculate density (defaults to 1000 for quasirandom and pseudorandom and 100 for others)
#' @param adjust adjust the bandwidth used to calculate the kernel density (smaller values mean tighter fit, larger values looser fit, default is 1)
#' @export
#' @rdname offsetX 
# @seealso \code{\link{offsetX}}, \code{\link[stats]{density}}
# @return a vector with of x-offsets between -1 and 1 of the same length as y
offsetSingleGroup<-function(y,maxLength=NULL,method=c('quasirandom','pseudorandom','smiley','frowney'),nbins=NULL,adjust=1) {
  method<-match.arg(method)
  if(is.null(nbins))nbins<-ifelse(method %in% c("pseudorandom","quasirandom"),2^10,ceiling(length(y)/5))
  #catch 0 length inputs
  if (length(y) == 0) return(NULL)
  # If there's only one value in this group, leave it alone
  if (length(y) == 1) return(0)

  #sqrt to match boxplot (allows comparison of order of magnitude different ns, scale with standard error)
  if(is.null(maxLength)||maxLength<=0)subgroup_width <- 1
  else subgroup_width <- sqrt(length(y)/maxLength)

  dens <- stats::density(y, n = nbins, adjust = adjust)
  dens$y <- dens$y / max(dens$y)
  offset <- switch(method,
    'quasirandom'=vanDerCorput(length(y))[rank(y, ties.method="first")],
    'pseudorandom'=stats::runif(length(y)),
    'smiley'=stats::ave(y,as.character(cut(y,dens$x)),FUN=topBottomDistribute),
    'frowney'=stats::ave(y,as.character(cut(y,dens$x)),FUN=function(x)topBottomDistribute(x,frowney=TRUE)),
    stop(simpleError('Unrecognized method in offsetSingleGroup'))
  )

  pointDensities<-stats::approx(dens$x,dens$y,y)$y

  #*2 to get -1 to 1
  out<-(offset-.5)*2*pointDensities*subgroup_width

  return(out)
}

#' Produce offsets such that points are sorted with most extreme values to right and left
#' 
#' Produce offsets to generate smile-like or frown-like distributions of points. That is sorting the points so that the most extreme values alternate between the left and right e.g. (max,3rd max,...,4th max, 2nd max). The function returns either a proportion between 0 and 1 (useful for plotting) or an order
#' 
#' @param x the elements to be sorted
#' @param frowney if TRUE then sort minimums to the outside, otherwise sort maximums to the outside
#' @param prop if FALSE then return an ordering of the data with extremes on the outside. If TRUE then return a sequence between 0 and 1 sorted by the ordering
#' @return a vector of the same length as x with values ranging between 0 and 1 if prop is TRUE or an ordering of 1 to length(x)
#' @export
#' @examples
#' topBottomDistribute(1:10)
#' topBottomDistribute(1:10,TRUE)
topBottomDistribute<-function(x,frowney=FALSE,prop=TRUE){
  if(length(x)==1)return(.5)
  if(frowney)x<- -x
  newOrder<-rank(x,ties.method='first')
  newOrder[newOrder%%2==1]<- -newOrder[newOrder%%2==1]
  newOrder<-rank(newOrder)
  if(prop){
    props<-seq(0,1,length.out=length(newOrder))
    newOrder<-props[newOrder]
  }
  return(newOrder)
}

#' Generate van der Corput sequences
#' 
#' Generates the first (or an arbitrary offset) n elements of the van der Corput low-discrepancy sequence for a given base
#' 
#' @param n the first n elements of the van der Corput sequence
#' @param base the base to use for calculating the van der Corput sequence
#' @param start start at this position in the sequence
#' @return a vector of length n with values ranging between 0 and 1
#' @references \url{https://en.wikipedia.org/wiki/Van_der_Corput_sequence}
#' @export
#' @examples
#' vanDerCorput(100)
vanDerCorput <- function(n, base=2,start=1){
  #generate n first digits of the van der Corput sequence
  if(n==0)return(c())
  if(n<0)stop(simpleError('n < 0 in vanDerCorput'))
  if(base<=1)stop(simpleError('base <=1 in vanDerCorput'))
  if(start<1)stop(simpleError('start < 1 in vanDerCorput'))
  out<-sapply(1:n+start-1,function(ii)digits2number(rev(number2digits(ii,base)),base,TRUE))
  return(out)
}

#' Convert an integer to an arbitrary base
#' 
#' Takes an integer and converts it into an arbitrary base e.g. binary or octal. Note that the first digit in the output is the least significant.
#' 
#' @param n the integer to be converted
#' @param base the base for the numeral system (e.g. 2 for binary or 8 for octal)
#' @return a vector of length ceiling(log(n+1,base)) respresenting each digit for that numeral system
#' @references \url{https://en.wikipedia.org/wiki/Radix}
#' @export
#' @examples
#' number2digits(100)
#' number2digits(100,8)
number2digits <- function(n, base=2){
  if(n==0)return(c())
  if(n<0)stop(simpleError('negative number in number2digits'))
  if(base<=1)stop(simpleError('base <=1 in number2digits'))
  nDigits<-ceiling(log(n+1,base))
  powers<-base^(0:nDigits)
  out<-diff(n %% powers)/powers[-length(powers)]
  return(out)
}

#' Convert a vector of integers representing digits in an arbitrary base to an integer
#' 
#' Takes a vector of integers representing digits in an arbitrary base e.g. binary or octal and converts it into an integer (or the integer divided by base^length(digits) for the number of digits if fractional is TRUE). Note that the first digit in the input is the least significant.
#' 
#' @param digits a vector of integers representing digits in an arbitrary base
#' @param base the base for the numeral system (e.g. 2 for binary or 8 for octal)
#' @param fractional divide the 
#' @return an integer
#' @references \url{https://en.wikipedia.org/wiki/Radix}
#' @export
#' @examples
#' digits2number(c(4,4,1),8)
#' digits2number(number2digits(100))
digits2number<-function(digits,base=2,fractional=FALSE){
  if(length(digits)==0)return(0)
  if(base<=0)stop(simpleError('base <= 0 in digits2number'))
  if(any(digits<0))stop(simpleError('digit < 0 in digits2number'))
  powers<-0:(length(digits)-1)
  out<-sum(digits*base^powers)
  if(fractional)out<-out/base^(length(digits))
  return(out)
}

#' the ave() function but with arguments passed to FUN
#' 
#' A function is applied to subsets of \code{x} where each subset consist of those observations with the same groupings in \code{y}
#'
#' @param x a vector to apply FUN to
#' @param y a vector or list of vectors of grouping variables all of the same length as \code{x}
#' @param FUN function to apply for each factor level combination.
#' @param ... additional arguments to \code{FUN}
#' @return A numeric vector of the same length as \code{x} where an each element contains the output from \code{FUN} after \code{FUN} was applied on the corresponding subgroup for that element (repeated if necessary within a subgroup).
#' @seealso \code{\link{ave}}
#' @export
#' @examples
#' aveWithArgs(1:10,rep(1:5,2))
#' aveWithArgs(c(1:9,NA),rep(1:5,2),max,na.rm=TRUE)
aveWithArgs<-function(x, y, FUN = mean,...){
  if (missing(y))
    x[] <- FUN(x,...)
  else {
    g <- interaction(y)
    split(x, g) <- lapply(split(x, g), FUN,...)
  }
  x
}


#' Data on HIV integration sites from several studies
#'
#' A dataset containing data from a meta-analysis looking for differences between active and inactive HIV integrations. Each row represents a provirus integrated somewhere in a human chromosome with whether viral expression was detectd, the distance to the nearest gene and the number of reads from H4K12ac ChIP-Seq mapped to within 50,000 bases of the integration.
#'
#' @format A data frame with 12436 rows and 4 variables:
#' \describe{
#'   \item{study}{the cell population infected by HIV}
#'   \item{latent}{whether the provirus was active (expressed) or inactive (latent)}
#'   \item{nearestGene}{distance to nearest gene (transcription unit) (0 if in a gene)}
#'   \item{H4K12ac}{number of reads aligned within +- 50,000 bases in a H4K12ac ChIP-Seq}
#' }
#' @references \url{http://www.retrovirology.com/content/10/1/90}
#' @source \url{http://www.retrovirology.com/content/10/1/90/additional}, system.file("data-raw", "makeIntegrations.R", package = "vipor")
"integrations"
