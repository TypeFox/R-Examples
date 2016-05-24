## ************************************************************************
##
## Compute GM-Distance between discrete distributions, using method from:
##
## Systematic Clustering of Transcription Start Site Landscapes
## Author(s): Zhao et al, 2011
## Source: PLoS ONE 6(8): e23409. doi:10.1371/journal.pone.0023409 
## 
## URL: http://dx.plos.org/10.1371/journal.pone.0023409
##      http://www.plosone.org/article/info:doi/10.1371/journal.pone.0023409
##
##
## (c) Xiaobei Zhao 2011 <xiaobei@binf.ku.dk>
## 
## ************************************************************************

##' Generalized Minimum Distance between a pair of distributions
##'
##' Generalized Minimum Distance between a pair of distributions
##' @title Generalized Minimum Distance between a pair of distributions
##' @aliases gmdp print.gmdp summary.gmdp
##' @usage
##'
##' gmdp(v1, v2, labels=c("v1","v2"), pseudocount=0, sliding=TRUE,
##' resolution=1)
##'
##' \method{print}{gmdp}(x, mode=c("brief","detailed","full"),
##' digits=3, ...)
##'
##' \method{summary}{gmdp}(object, ...)
##'
##' 
##' @param v1 a numeric vector, giving positional counts as a discrete distribution.
##' @param v2 a numeric vector, giving positional counts as a discrete distribution.
##' @param labels a string vector of length 2, giving the names of v1 and v2 respectively.
##' @param pseudocount a numeric value to be allocated for each position to reduce bias;
##' by default \code{pseudocount = 0}.
##' @param sliding logical, indicating whether sliding is allowed or not for an optimal solution; 
##' by default \code{sliding = TRUE}.
##' @param resolution relative resolution, numeric (>=1), changing the size of the bin by
##' multiplying the value. A larger value (lower resolution) is more computational efficiet
##' but missing details.
##' @param x an object of class \code{gmdp}.
##' @param object an object of class \code{gmdp}.
##' @param mode a string of the following: \code{c("brief","detailed","full")},
##' indicating whether to print in \emph{full} mode (\emph{default}). 
##' @param digits integer, indicating the number of decimal places to be printed.
##' @param ... arguments to be passed to method.
##' @return \code{gmdp} returns an object of class \code{gmdp}, a numeric with an attribute of
##' \emph{meta} in a list with components: 
##'
##' labels: a string vector, giving the names of distributions
##' 
##' v1.ori: a numeric vector, the first input distribution
##' 
##' v2.ori: a numeric vector, the second input distribution
##' 
##' v1: a numeric vector, the normalized version of the first input distribution
##' 
##' v2: a numeric vector, the normalized version of the second input distribution
##' 
##' distance: numeric, the \emph{GM-Distance} (\emph{GMD})
##' 
##' sliding: logical, indicating whether sliding is performed
##' 
##' pseudocount: a numeric value that is allocated at each position in addition to original values
##' 
##' gap.pair: a numeric matrix, giving one gap pair per row:
##' i.e. relative shifts between distributions of one optimal hit
##' 
##' n.hit: numeric, the number of (equally good) optimal hits
##' @references See \code{citation("GMD")}
##' @seealso \code{\link{print.gmdp}}, \code{\link{summary.gmdp}}, \code{\link{plot.gmdp}}
##' \code{\link{gmdm}} 
##' @keywords classes
##' @examples
##' require(GMD)
##' gmdp(c(4,1,1,0,0,0,3,1),c(2,1,1,0,0,0,3,3),sliding=FALSE)
##' x <- gmdp(c(4,1,1,0,0,0,3,1), c(1,1,2,1,1,0,0,0,3,3,5,5),
##' pseudocount=1, sliding=TRUE)
##' print(x)
##' print(x, "full")
gmdp <-
  function(v1, v2, labels=c("v1","v2"), pseudocount=0, sliding=TRUE, resolution=1)
{
  if (!is.numeric(v1) | !is.numeric(v2)){
    stop("`v1' and `v2' should be `numeric' vectors.")
  }
  if (resolution<1){
    stop("`resolution' should be a `numeric' value no less than 1.")
  }

  ## 
  meta <- list()
  meta$labels <- labels
  meta$v1.ori <- v1
  meta$v2.ori <- v2
  meta$resolution <- resolution
  
  ## resolve
  if(resolution>1){
    v1 <- .resolveHist(v1,resolution)
    v2 <- .resolveHist(v2,resolution)
  }

  ##
  l1 <- length(v1)
  l2 <- length(v2)

  ##
  if (sliding) {
    res <- .gmd(v1, v2, pseudocount)
    sliding <- TRUE
    gap <- which(res$position==1) - min(l1,l2) - 1
    n.hit <- sum(res$position==1)
  } else {
    if (l1 != l2){
      stop("The lengths of two input vectors should be equal.
Otherwise, please allow `sliding` (sliding=TRUE).\n"
           )
    } else{
      d <- .gmd0(v1, v2, pseudocount)
      p <- 0
      res <- list(distance=d, position=c(1))
      sliding <- FALSE
      gap <- 0
      n.hit <- 1 
    }
  }

  if (l1<=l2){
    gap.pair <- cbind(gap,0)
  } else {
    gap.pair <- cbind(0,gap)
  }
  colnames(gap.pair) <- c()
  gap.pair <- t(apply(gap.pair,1,function(x) if(any(x<0)) x-min(x) else x))

  meta$v1 <- v1/sum(v1)
  meta$v2 <- v2/sum(v2)
  meta$sliding <- sliding
  meta$pseudocount <- pseudocount
  meta$gap.pair <- gap.pair
  meta$n.hit <- n.hit

  ## return
  ret <- res$d
  class(ret) <- c("gmdp",class(ret))
  attr(ret, "meta") <- meta
  return(ret)
}



print.gmdp <- function(x, mode=c("brief","detailed","full"), digits=3, ...)
{
  if (!.is.gmdp(x)){
    stop("`x' should be an object of class `gmdp'.")
  }
  mode <- match.arg(mode)

  if (mode=="brief"){
    attr(x,"meta") <- c()
    class(x) <- class(x)[class(x)!="gmdp"]
    print(x, ...)
  } else if (mode=="detailed" | mode=="full"){
    x.meta <- attr(x,"meta")
    s.gap <-
      paste(sprintf("\t%s\t%s\n",x.meta$labels[1],x.meta$labels[2]),
            paste(sapply(1:nrow(x.meta$gap.pair),
                         function(i) sprintf("Hit%s\t%s\t%s",i, x.meta$gap.pair[i,1],x.meta$gap.pair[i,2])),
                  sep="",
                  collapse="\n"),
            "\n",
            sep=""
            )
    cat("\n")
    if (mode=="full"){
      cat(sprintf("Distribution of %s:\n%s\n",x.meta$labels[1],paste(round(x.meta$v1.ori,digits),sep="",collapse=" ")))
      cat(sprintf("(After normalization)\n%s\n\n",paste(round(x.meta$v1,digits),sep="",collapse=" ")))
      cat(sprintf("Distribution of %s:\n%s\n",x.meta$labels[2],paste(round(x.meta$v2.ori,digits),sep="",collapse=" ")))
      cat(sprintf("(After normalization)\n%s\n\n",paste(round(x.meta$v2,digits),sep="",collapse=" ")))
    }
    cat(sprintf(sprintf("GM-Distance: %%.%sf\n\n",digits),x))
    cat(sprintf("Sliding: %s\n\n",x.meta$sliding))
    cat(sprintf("Number of hits: %s\n\n",x.meta$n.hit))
    cat(sprintf("Gap:\n%s\n\n",s.gap))
    cat(sprintf("Resolution: %s\n\n",x.meta$resolution))
    cat("\n")
  }
}



summary.gmdp <- function(object, ...){
  if (!.is.gmdp(object)){
    stop("`object' should be an object of class `gmdp'.")
  }
  class(object) <- "list"
  summary(object, ...)
}



##' Plot Function for Class \code{gmdp}
##'
##' Plot Function for Class \code{gmdp}
##' @title Plot function for class gmdp
##' 
##' @export plot.gmdp
##' @param x an object of class \code{gmdp}.
##' @param labels a string vector of the same length of \code{x$labels},
##' giving the names of the numeric vectors in \code{x}.
##' @param colors the colors of the discrete distributions.
##' See \code{help("plot.mhist", package="GMD")}.
##' @param main an overall title for the plot. 
##' @param ylab a title for the y axis.
##' See \code{help("plot.mhist", package="GMD")}.
##' @param xlab a title for the x axis.
##' See \code{help("plot.mhist", package="GMD")}.
##' @param xlim numeric vectors of length 2, giving the x coordinates ranges.
##' @param type type of plot, as in \code{help("plot", package="graphics")}.
##' @param if.text.gmd logical, indicating whether \emph{GM-Distance} is reported in the subtitle.
##' @param if.text.gap logical, indicating whether \emph{gap} is reported in the subtitle.
##' @param ... arguments to be passed to methods.
##' See \code{help("plot.mhist", package="GMD")}.
##' @references See \code{help(GMD)}
##' @seealso \code{\link{gmdp}}
##' @keywords methods hplot
##' @examples
##' require("GMD") # load library
##' data(cage)     # load data
##'
##' ## measure pairwise distance
##' x <- gmdp(cage[["Pfkfb3 (T02R00AEC2D8)"]],cage[["Csf1 (T03R0672174D)"]])
##' print(x)                     # print a brief version by default
##' print(x, mode="full")  # print a full version by default
##' 
##' ## show alignment
##' plot(x,labels=c("Pfkfb3","Csf1"),beside=FALSE)
##' 
##' ## show another alignment
##' plot(gmdp(cage[["Hig1 (T09R0743763C)"]],cage[["Cd72 (T04R028B8BC9)"]]),
##'      labels=c("Hig1 (T09R0743763C)","Cd72 (T04R028B8BC9)"),
##'      beside=FALSE)
plot.gmdp <-
  function(x,
           ##
           labels=NULL,
           colors=NULL,
           ##
           main,
           ylab="Fraction",
           xlab="Position",
           ##
           xlim=NULL,
           type=NULL,
           ##
           if.text.gmd=TRUE,
           if.text.gap=TRUE,
           ...
           )
{
  ## --
  ## Top: v1
  ## Middle: alignment
  ## Bottom: v2
  ## --
  if (!.is.gmdp(x)){
    stop("`x' should be an object of class `gmdp'.")
  }

  ## meta
  x.meta <- attr(x,"meta")

  ## labels
  if (.invalid(labels)){
    labels <- x.meta$labels
  }

  ## check multiple hits ##
  if (nrow(x.meta$gap.pair)>1){
    warning("There are multiple optimal hits/alignment; the first one will be plotted.")
  }
  gap.pair <- x.meta$gap.pair[1,] ## pick the first if multiple cases
  s.gap <- sprintf("gap=%s",sprintf("c(%s,%s)",gap.pair[1],gap.pair[2]))

  ## title and subtitle ##
  if (.invalid(main)){
    main <-
      sprintf("Optimal alignment between distributions (%s sliding)",
              ifelse(x.meta$sliding,"with","without")
              )
  }
  
  s.gmd <- sprintf("GMD=%.3f",x)
  
  if (if.text.gmd & if.text.gap){
    s.sub <- sprintf("%s, %s",s.gmd,s.gap)
  } else if (if.text.gmd){
    s.sub <- s.gmd
  } else if (if.text.gap){
    s.sub <- s.gap
  } else {
    s.sub <- ""
  }

  ## plot.mhist ##
  tmp.v1 <- x.meta$v1
  tmp.v2 <- x.meta$v2

  ##tmp.v1 <- c(rep(NA,gap.pair[1]),x.meta$v1,rep(NA,gap.pair[2]))
  ##tmp.v2 <- c(rep(NA,gap.pair[2]),x.meta$v2,rep(NA,gap.pair[1]))
  
  if (gap.pair[1]){
    tmp.v1 <- c(rep(NA,gap.pair[1]),x.meta$v1)
    tmp.v2 <- c(x.meta$v2,rep(NA,gap.pair[1]))
  } else if (gap.pair[2]){
    tmp.v1 <- c(x.meta$v1,rep(NA,gap.pair[2]))
    tmp.v2 <- c(rep(NA,gap.pair[2]),x.meta$v2)
  }
  tmp.x <- list(tmp.v1,tmp.v2)
  names(tmp.x) <- labels

  if (.invalid(xlim)) {
    xlim <- c(1,max(unlist(lapply(tmp.x,length))))
  }
  
  ## plot
  ##print(position.lim)
  tmp.plot <-
    plot.mhist(as.mhist(tmp.x),
               labels=labels,
               colors=colors,
               main=main,
               sub=s.sub,
               ylab=ylab,
               xlab=xlab,
               xlim=xlim,
               type=type,
               ...
               )
  
 
  invisible(tmp.plot)
  
}
