###NAMESPACE ADDITIONS###
# Depends: R (>= 2.10), grDevices, graphics, stats, utils, reader
# Imports: tools, proftools, plyr
# Suggests: KernSmooth, BiocInstaller, Matrix
# importFrom(proftools, readProfileData, flatProfile)
# importFrom(tools, toHTML)
# import(grDevices, graphics, stats, utils)
###END NAMESPACE###

# importFrom(dplyr, filter)


## add check.bio() to internals list
# no longer want to: importFrom(BiocInstaller, biocVersion)


# 7 new not in index - NEW! - 


#' Determine whether a function can be applied to an S4 class/object
#' 
#' Wrapper for 'showMethods', allows easy testing whether a function
#' (can be specified as a string, or the actual function itself (FUN)) can be
#' applied to a specific object or class of objects (CLASS)
#' @param FUN the function to test, can be specified as a string, or the actual function itself
#' @param CLASS  a specific object or a class of objects specified by a string, e.g, "GRanges"
#' @param false.if.error logical, the default value is FALSE, in which case an error is returned
#' when FUN is not an S4 generic function. If this parameter is set to TRUE, 'FALSE' will
#' be returned with a warning instead of an error.
#' @param ... additional arguments to showMethods(), e.g, 'where' to specify the environment
#' @export
#' @return returns logical (TRUE/FALSE), or if the function is not S4 will return an error,
#' although this could potentially be because the function's package has not been loaded.
#' @examples
#' require(Matrix); require(methods)
#' has.method("t","dgeMatrix") # t() is the transpose method for a dgeMatrix object
#' has.method(t,"dgeMatrix") # also works without quotes for the method
#' m.example <- as(matrix(rnorm(100),ncol=5),"dgeMatrix")
#' has.method(t, m.example) # works with an instance of an object type too
#' has.method("band", m.example) # band is a function for a 'denseMatrix' but not 'dgeMatrix'
#' ## not run # has.method("notAFunction","GRanges") # should return error
#' has.method("notAFunction","GRanges",TRUE) # should return FALSE and a warning
has.method <- function(FUN,CLASS, false.if.error=FALSE, ...) {
  if(!is.character(CLASS)) { CLASS <- class(CLASS) }
  if(!is.character(FUN) & !is.function(FUN)) { 
    if(false.if.error) {
      warning("FUN should be an R function, as a string or function, returning 'FALSE'")
      return(FALSE)
    } else {
      stop("FUN must be an R function, as a string or function") 
    }
  }
  test <- showMethods(FUN,classes=CLASS,printTo=FALSE,...)
  if(length(grep("not an S4 generic function",test))>0) {
    if(false.if.error) {
      warning("'",FUN,"' was not an S4 generic function or required package not loaded, returning 'FALSE'")
      return(FALSE)
    } else {
      stop("'",FUN,"' was not an S4 generic function or required package not loaded")
    }
  }
  return(!(length(grep("No methods",test))>0))
}


#' Function to add commas for large numbers
#' 
#' Often for nice presentation of genomic locations it is helpful
#' to insert commas every 3 digits when numbers are large. This function
#' makes it simple and allows specification of digits if a decimal number
#' is in use.
#' @param x a vector of numbers, either as character, integer or numeric form
#' @param digits integer, if decimal numbers are in use, how many digits to display, 
#' same as input to base::round()
#' @return returns a character vector with commas inserted every 3 digits
#' @export
#' @examples
#' comify("23432")
#' comify(x=c(1,25,306,999,1000,43434,732454,65372345326))
#' comify(23432.123456)
#' comify(23432.123456,digits=0)
comify <- function(x,digits=2) {
  if(length(Dim(x))>1) { stop("x must be a vector") }
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           if(length(x)>1) { return(sapply(x, comify, digits=digits)) }
  x <- round(as.numeric(x),digits=digits)
  x <- (paste(x)); dec <- ""
  if(any(grep(".",x))) {
    x.plus.dp <- strsplit(x,".",fixed=TRUE)[[1]]
    if(length(x.plus.dp)>2) { stop("x contained invalid decimal point(s)") }
    xx <- x.plus.dp[1]
    if(length(x.plus.dp)==2) { dec <- paste(".",x.plus.dp[2],sep="") }
  } else { xx <- x }
  splt <- strsplit(xx,"")[[1]]
  nm <- rev(splt)
  cnt <- 0; new <- NULL
  LL <- length(nm)
  for (cc in 1:LL) {
    new <- c(nm[cc],new)
    cnt <- cnt+1
    if(cnt>2 & cc!=LL) { new <- c(",",new); cnt <- 0 }
  }
  return(paste(paste(new,collapse=""),dec,sep=""))
}



#' Convert p-values to Z-scores
#' 
#' Simple conversion of two-tailed p-values to Z-scores. Written
#' in a way that allows maximum precision for small p-values.
#' @param p p-values (between 0 and 1), numeric, scalar, vector or matrix, 
#' or other types coercible using as.numeric()
#' @return Z scores with the same dimension as the input
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{Z.to.p}}
#' @examples
#' p.to.Z(0.0001)
#' p.to.Z("5E-8")
#' p.to.Z(c(".05",".01",".005"))
#' p.to.Z(matrix(runif(16),nrow=4))
p.to.Z <- function(p) { 
  if(!is.numeric(p)) { p <- as.numeric(p) }
  if(!is.numeric(p)) { stop("p was not coercible to numeric type") }
  ll <- length(which(p<0 | p>1))
  if(ll>0) { warning(ll, " invalid p-values set to NA"); p[p<0 | p>1] <- NA }
  O <- qnorm((p/2),F)
  O[!is.finite(O)] <- NA
  return(-O) 
}

#' Convert Z-scores to p-values
#' 
#' Simple conversion of Z-scores to two-tailed p-values. Written
#' in a way that allows maximum precision for small p-values.
#' @param Z Z score, numeric, scalar, vector or matrix, or other types coercible
#'  using as.numeric()
#' @param warn logical, whether to give a warning for very low p-values when
#' precision limits are exceeded.
#' @return p-valuues with the same dimension as the input
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{p.to.Z}}
#' @examples
#' Z.to.p("1.96")
#' Z.to.p(p.to.Z(0.0001))
#' Z.to.p(37, TRUE)
#' Z.to.p(39, TRUE) # maximum precision exceeded, warnings on
#' Z.to.p(39) # maximum precision exceeded, warnings off
Z.to.p <- function(Z, warn=FALSE) {
  if(!is.numeric(Z)) { Z <- as.numeric(Z) }
  if(!is.numeric(Z)) { stop("Z was not coercible to numeric type") }
  if(any(abs(Z)>=38) & warn) { warning("maximum precision exceeded, p < 10^-300") }
  O <- 2*pnorm(-abs(Z))
  O[!is.finite(O)] <- NA
  return(O) 
}


#' Posterior probability of association function
#'
#' Estimate the probability of your hypothesis being true,
#' given the observed p-value and a prior probability of
#' the hypothesis being true.
#' @param p p-value you want to test [p<0.367], or 'bayes factor'
#' @param prior prior odds for the hypothesis (Ha) being tested
#' @param BF logical, set to TRUE if you have entered a bayes factor
#' as 'p' rather than a p-value
#' @param quiet logical, whether to display verbose information for
#' calculation
#' @return prints calculations, then returns the posterior 
#' probability of association given the observed p-value 
#' under the specified prior
#' @export
#' @references
#' Equations 1, 2 from
#' http://www.readcube.com/articles/10.1038/nrg2615
#' Equations 2, 3 from
#' http://www.tandfonline.com/doi/pdf/10.1198/000313001300339950
#' @examples
#' ps <- rep(c(.05,.01),3)
#' prs <- rep(c(.05,.50,.90),each=2)
#' mapply(ps,prs,FUN=ppa)  # replicate Nuzzo 2014 table
#' # try with bayes factors
#' ppa(BF=3,prior=.9)
#' ppa(BF=10,prior=.5)
ppa <- function(p=.05, prior=.5, BF=NULL, quiet=TRUE) {
  if(any(p<=0 | p>=(1/exp(1)))) { stop("invalid p value") }
  if(any(prior<=0 | prior>=(1))) { stop("invalid prior") }
  if(is.null(BF)) { 
    # calculate bayes factors from p, if BF not entered
    if(!quiet) { cat("\np value:",p,"with prior:",prior,"\n") }
    BF <- (-exp(1)*(p)*log(p) )^(-1)
    # NB: ^invert BF so in terms of % support for Ha 
  } else { 
    if(!quiet) { cat("\nprior:",prior,"with ") }
    if(any(BF<0)) { stop("invalid bayes factor (BF)") }
  }
  if(!quiet) { cat("bayes factor:",BF,"\n") }
  P0 <- (prior/(1-prior)) * (BF) 
  if(!quiet) { cat("posterior odds = bayes factor * H1/H0 prior:",P0,"\n") }
  ppa <- (P0/(1+P0)) 
  if(!quiet) { cat("posterior probability of association:",ppa,"\n") }
  return(ppa)
}



#' Return vector indexes of statistical univariate outliers
#'
#' Performs simplistic outlier detection and returns indexes for outliers.
#' Acts like the which() function, return indices of elements of a vector
#' satisfying the condition, which by default are outliers exceeding 2 SD
#' above or below the mean. However, the threshold can be specified, only
#' high or low values can be considered outliers, and percentile and interquartile
#' range thresholds can also be used.
#'
#' @param x numeric, or coercible, the vector to test for outliers
#' @param thr numeric, threshold for cutoff, e.g, when method="sd", standard deviations,
#' when 'iq', interquartile ranges (thr=1.5 is most typical here), or when 'pc', you might
#' select the extreme 1\%, 5\%, etc.
#' @param method character, one of "sd","iq" or "pc", selecting whether to test for outliers
#' by standard deviation, interquartile range, or percentile.
#' @param high logical, whether to test for outliers greater than the mean
#' @param low logical, whether to test for outliers less than the mean
#' @return indexes of the vector x that are outliers according to either
#' a SD cutoff, interquartile range, or percentile threshold, above (high) and/or
#' below (low) the mean/median.
#' @export
#' @examples
#' test.vec <- rnorm(200)
#' summary(test.vec)
#' ii <- which.outlier(test.vec) # 2 SD outliers
#' prv(ii); vals <- test.vec[ii]; prv(vals)
#' ii <- which.outlier(test.vec,1.5,"iq") # e.g, 'stars' on a box-plot
#' prv(ii)
#' ii <- which.outlier(test.vec,5,"pc",low=FALSE) # only outliers >mean
#' prv(ii)
which.outlier <- function(x, thr=2, method=c("sd","iq","pc"), high=TRUE, low=TRUE) {
  if(!is.numeric(x)) { x <- as.numeric(x) }
  if(!is.numeric(thr)) { stop("thr must be numeric") }
  X <- x
  #  x <- X
  X <- narm(X)
  if(!is.numeric(x)) { stop("x must be numeric, or coercible to numeric") }
  X <- X[is.finite(X)]
  if(length(X)>1) {
    method <- substr(tolower(method),1,2)[1]
    if(!method %in% c("sd","iq","pc")) { stop("invalid method, must be sd [std dev], iq [interquartile range], or pc [percentile]") }
    if(method=="sd") {
      stat <- sd(X,na.rm=T)
      hi.thr <- mean(X,na.rm=T) + stat*thr
      lo.thr <- mean(X,na.rm=T) - stat*thr
    } else {
      if(method=="iq") {
        sl <- summary(X)
        stat <- (sl[5]-sl[2])
        hi.thr <- median(X,na.rm=T) + stat*thr
        lo.thr <- median(X,na.rm=T) - stat*thr
      } else {
        stat <- pctile(X,pc=force.percentage(thr))
        hi.thr <- stat[2] ; lo.thr <- stat[1]
      }
    }
    if(high) {
      outz <- X[X>hi.thr]
    } else { outz <- NULL }
    if(low) {
      outz <- unique(c(outz,X[X<lo.thr]))
    }
    outz <- which(x %in% outz) # make sure indexes include the NA, Inf values
    return(outz)
  } else {
    warning("outlier detection requires more than 1 datapoint")
    return(numeric(0))
  }
}



#' Obtain an index of all instances of values with duplicates (ordered)
#' 
#' The standard 'duplicated' function, called with which(duplicated(x)) will 
#' only return the indexes of the extra values, not the first instances. For instance
#' in the sequence: A,B,A,C,D,B,E; it would return: 3,6. This function will also
#' return the first instances, so in this example would give: 1,3,2,6 [note it
#' will also be ordered]. This index can be helpful for diagnosis if duplicates 
#' are unexpected, for instance in a data.frame, and you wish to compare the differences
#' between the rows with the duplicate values occuring. Also, duplicate values are sorted
#' to be together in the listing, which can help for manual troubleshooting of undesired
#' duplicates.
#' @param x a vector that you wish to extract duplicates from
#' @return vector of indices of which values in 'x' are duplicates (including
#' the first observed value in pairs, or sets of >2), ordered by set, then
#' by appearance in x.
#' @export
#' @examples
#' set <- c(1,1,2,2,3,4,5,6,2,2,2,2,12,1,3,3,1)
#' dup.pairs(set) # shows the indexes (ordered) of duplicated values
#' set[dup.pairs(set)] # shows the values that were duplicated (only 1's, 2's and 3's)
dup.pairs <- function(x) {
    dx <- duplicated(x)
    other.dups <- which(dx)
    not.and.first <- which(!dx)
    ind.dups <- not.and.first[which(x[not.and.first] %in% x[other.dups])]
    xo <- x[other.dups]
    vc <- vector()
    for (cc in 1:length(ind.dups)) {
      vc <- c(vc,ind.dups[cc],other.dups[which(xo %in% x[ind.dups[cc]])])
    }
    return(vc)
}


#' Create variables from a list
#' 
#' Places named objects in a list into the working environment as individual variables.
#' Can be particularly helpful when you want to call a function that produces a list of
#' multiple return variables; this gives a way to access them all at once in the environment
#' from which the function was called.
#' @param list list, with named objects, each element will become a named variable in
#' the current environment
#' @return New variables will be added to the current environment. Use with care as any 
#' already existing with the same name will be overwritten.
#' @export
#' @seealso base::list2env
#' @examples
#' list.to.env(list(myChar="a string", myNum=1234, myList=list("list within a list",c(1,2,3))))
#' print(myChar)
#' print(myNum)
#' print(myList)
#' two.arg.return <- function(X) { return(list(Y=X+1,Z=X*10)) } 
#' result <- two.arg.return(11) # function returns list with 2 variables
#' list.to.env(result)
#' print(Y); print(Z)
list.to.env <- function(list) {
  if(!is.list(list)) { stop("this function's sole parameter must be a list object")}
  if(is.null(names(list))) { stop("list elements must be named") }
  if(length(list)>1000) { warning("list contains over 1000 elements, this operation will crowd the workspace") }
  for(cc in 1:length(list)) {
    assign(x=names(list)[cc],value=list[[cc]],pos=parent.frame())
  }
  return(NULL)
}

  
#' Simple representation and retrieval of Date/Time
#'
#' Retrieve a simple representation of date_time or just date, 
#' for generating day/time specific file names, etc.
#' @param sep character, separator to use for the date/time, eg, 
#' underscore or <space> " ".
#' @param long logical, whether to display a longer version of the
#' date and time, or just a simple version
#' @param time logical, whether to include the time, or just the date
#' @export
#' @return A string containing the date: MMMDD and optionally time HRam/pm.
#' Or if long=TRUE, a longer representation: DAY MM DD HH.MM.SS YYYY.
#' @examples
#' simple.date()
#' simple.date(" ",long=TRUE)
#' simple.date(time=FALSE)
simple.date <- function(sep="_",long=FALSE,time=TRUE) {
  myt <- format(Sys.time(), "%a %b %d %X %Y")
  if(long) {return(gsub(":",".",gsub(" ",sep,myt))) }
  dt <- strsplit(myt,":",fixed=TRUE)[[1]][1]
  splt <- strsplit(dt," ")[[1]]
  tm <- as.numeric(tail(splt,1))
  pr.tm <- head(splt,length(splt)-1)
  pr.tm[2] <- toupper(pr.tm[2])
  ampm <- {if(as.numeric(tm)>11) {"PM"} else {"AM"}}
  tm <- {if(tm>12) { tm-12 } else { if(tm<1) { tm+12 } else { tm } }}
  if(nchar(paste(pr.tm[3]))==1) { pr.tm[3] <- paste("0",pr.tm[3],sep="" ) }
  if(!time) { out <- paste(pr.tm[-1],collapse="") } else {
    out <- paste(paste(pr.tm[-1],collapse=""),sep,
                 tm, ampm ,sep="") }
  return(out)
}




#' Easily display fraction and percentages
#' 
#' For a subset 'n' and total 'N', nicely prints text n/N and/or percentage%.
#' Often we want to display proportions and this simple function reduces the
#' required amount of code for fraction and percentage reporting. If 
#' insufficient digits are provided small percentage may truncate to zero.
#' @param n numeric, the count for the subset of N (the numerator)
#' @param N numeric, the total size of the full set (the denominator)
#' @param digits, integer, the number of digits to display in the percentage
#' @param pc, logical, whether to display the percentage of N that n comprises
#' @param oo, logical, whether to display n/N as a fraction
#' @param use.sci, logical, whether to allow scientific notation for small/large
#' percentages.
#' @return A string showing the fraction n/N and percentage (or just one of these)
#' @export
#' @examples
#' out.of(345,12144)
#' out.of(345,12144,pc=FALSE)
#' out.of(3,10^6,digits=6,oo=FALSE)
#' out.of(3,10^6,digits=6,oo=FALSE,use.sci=TRUE)
out.of <- function(n,N=100,digits=2,pc=TRUE,oo=TRUE,use.sci=FALSE) {
  pct <- 100*(n/N)
  outof <- paste(n,"/",N,sep="")
  if(use.sci) {
    percent <- paste(round(pct,digits),"%",sep="")
  } else {
    percent <- paste(format(round(pct,digits),scientific=FALSE),"%",sep="")
  }  
  if(pc & oo) {
    outof <- paste(outof," (",percent,")",sep="")
  } else {
    if(pc) { outof <- percent }
  }
  return(outof)
}




#' Extend an interval by percentage
#' 
#' For various reasons, such as applying windows, setting custom range limits for plots, it may 
#' be desirable to extend an interval by a certain percentage.
#' @param X a numeric range, should be length 2. If a longer numeric, will be coerced with range()
#' @param pc percentage by which to extend X, can be entered in either percentage style: 0<pc<1; 
#' or 1<pc<100
#' @param swap logical, if TRUE, flip the extension directions if X[2]<X[1], ie, not in numerical
#' order
#' @param pos logical, if TRUE, make an extension in the positive direction
#' @param neg logical, if TRUE, make an extension in the negative direction
#' @export 
#' @examples
#' extend.pc(c(2,10),0.25) # extend X symmetrically
#' extend.pc(c(2:10),0.25) # extend the range of X
#' # the following 3 examples extend X by 1% only in the 'positive' direction
#' extend.pc(c(25000,55000),.01,neg=FALSE) # standard positive extension
#' extend.pc(c(55000,25000),.01,neg=FALSE) # ranges in reverse order, not swapped
#' extend.pc(c(55000,25000),.01,neg=FALSE,swap=TRUE) # ranges in reverse order, swapped
extend.pc <- function(X,pc=.5,pos=TRUE,neg=TRUE,swap=FALSE) {
  if(!is.numeric(X)) { stop("X must be numeric") }
  if(length(X)==0) { stop("X was empty") }
  if(length(X)==1) { X <- c(X,X); warning("X was length=1, extended by repeating X twice") }
  if(length(X)>2) { X <- range(X) } #; warning("X was length>2, coerced using X <-range(X)") }
  pc <- force.percentage(pc)
  yn <- yp <- abs(X[2]-X[1])*pc
  if(!pos) { yp <- 0 }; if(!neg) { yn <- 0 }
  if(swap) {
    # flip the extension directions if X[2]<X[1], ie, not in numerical order
    if(X[1]>X[2]) { temp <- yn; yn <- yp; yp <- temp }
  }
  return(c(X[1]-yn,X[2]+yp))
}


#' Draw a scatterplot with a fit line
#'
#' Drawing a fit line usually requires some manual steps requiring several lines of code,
#' such as ensuring the data is sorted by x, and for some functions doesn't contain missing values.
#' This function takes care of these steps and automatically adds a loess fitline, or non-linear 
#' fitline. The type of scatter defaults to 'plot', but other scatter plot functions can be 
#' specified, such as graphics::smoothScatter(), for example. If 'file' is specifed, will 
#' automatically plot to a pdf of that name.
#' @param x data for the horizontal axis (independent variable)
#' @param y data for the vertical axis (dependent variable)
#' @param file file name for pdf export, leave as NULL if simply plotting to the GUI. File 
#' extension will be added automatically if missing
#' @param loess logical, if TRUE, fit using loess(), else use a polynomial fit
#' @param span numeric scalar, argument passed to the 'span' parameter of loess(), see ?loess for details
#' @param scatter function, by default is graphics::plot(), but any scatter-plot function of the 
#' form F(x,y,...) can be used, for example graphics::smoothScatter().
#' @param ylim numeric range for y axis, argument passed to plot(), see ?plot.
#' @param return.vectors logical, if TRUE, do not plot anything, just return the x and y coordinates
#' of the fit line as a list of vectors, x and y.
#' @param fit.col colour of the fit line
#' @param fit.lwd width of the fit line
#' @param fit.lty type of the fit line
#' @param fit.leg whether to include an automatic legend for the fit line (will alter the y-limits
#' to fit)
#' @param fit.r2 logical, whether to display r squared of the fit in the fit legend
#' @param ... further arguments to the plot function specified by 'scatter', e.g, 'main', 'xlab', etc
#' @export
#' @return if file is a character argument, plots data x,y to a file, else will generate a plot to
#' the current plotting environment/GUI. The display of the x,y points defaults to 'plot', but 
#' alternate scatter plot functions can be specified, such as graphics::smoothScatter() which used 
#' density smoothing, for example. Also, another option is to set return.vectors=TRUE, and then
#' the coordinates of the fit line will be returned, and no plot will be produced.
#' @examples
#' library(NCmisc)
#' require(KernSmooth)
#' DD <- sim.cor(1000,4) # create a simulated, correlated dataset
#' loess.scatter(DD[,3],DD[,4],loess=FALSE,bty="n",pch=".",cex=2)
#' loess.scatter(DD[,3],DD[,4],scatter=smoothScatter)
#' xy <- loess.scatter(DD[,3],DD[,4],return.vectors=TRUE)
#' prv(xy) # preview the vectors produced
loess.scatter <- function(x,y,file=NULL,loess=TRUE,span=0.75,scatter=plot,...,ylim=NULL,return.vectors=FALSE,
                          fit.col="red",fit.lwd=2,fit.lty="solid",fit.leg=TRUE,fit.r2=TRUE) {
  if(length(Dim(x))!=1 | length(Dim(y))!=1) { stop("x and y must be vectors") }
  if(length(x)<1 | length(y)<1) { warning("x/y must have positive length"); return(NULL) }
  if(!is.numeric(x) | !is.numeric(y)) { stop("x and y must be numeric") }
  if(length(x)!=length(y)) { stop("x and y must be vectors of the same length") }
  y1 <- y[order(x)]
  x1 <- x[order(x)]
  missing.either <- is.na(x1) | is.na(y1)
  if(length(which(missing.either))>0) { y1 <- y1[!missing.either]; x1 <- x1[!missing.either] }
  if(length(y1)<5) { 
    do.fit=F; warning("not enough points remain to generate plot with fit-line") 
  } else { do.fit <- T }
  if(do.fit) {
    if(!loess) {
      # if(all(x1>0)) {
      #   fit <- "non-linear"
      #   lo <- lm(y1~x1+sqrt(x1)+log(x1))
      # } else {
      fit <- "polynomial"
      lo <- lm(y1~x1 + (x1^2) + (x1^3) + (x1^4))
      # }
    } else {
      fit <- "loess"
      lo <- loess(y1~x1,span=span)
    }
    y2 <- predict(lo)
    if(fit.r2) {
      r2 <- round(cor(y1,y2,use="pairwise.complete"),3)
      leg.txt <- paste(fit,"fit line, r2 =",r2)
    } else {
      leg.txt <- paste(fit,"fit line")
    }
  }
  if(!return.vectors) {
    if(is.character(file)) { fnm <- cat.path("",fn=file[1],ext="pdf"); pdf(fnm) }
    if(fit.leg & do.fit) {
      y.range <- range(y1)
      ## this section of code allows a custom 'ylim' setting to override the internal ylim
      if(is.numeric(ylim) & length(ylim)==2) {
        if(ylim[1]>y.range[1] | ylim[2]<y.range[2]) { warning("ylim will truncate the y vector in the plot") }
        y.range <- ylim
      } 
      y.lims <- extend.pc(y.range,pc=0.25,neg=F)
      y.lims <- extend.pc(y.lims,pc=0.1,pos=F)
    } else { y.lims <- NULL }
    
    y <- y1; x <- x1 # ensures default x,y labels are x,y
    scatter(x,y,...,ylim=y.lims)
    if(do.fit) { lines(x1,y2,col=fit.col,lwd=fit.lwd,lty=fit.lty) }
    if(fit.leg & do.fit) { legend("topright",legend=leg.txt, lty=fit.lty, col=fit.col, lwd=fit.lwd, bty="n") }
    if(is.character(file)) { cat("wrote file",file,"\n"); dev.off() }
  } else {
    return(list(x=x1,y=y2))
  }
}

  



#' Return up to 22 distinct colours.
#' 
#' Useful if you want to colour 22 autosomes, etc, because most R
#' colour palettes only provide 12 or fewer colours, or else provide,
#' a gradient which is not distinguishable for discrete categories.
#' Manually curated so the most similar colours aren't side by side.
#'
#' @param n number of unique colours to return
#' @return returns vector of n colours
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' get.distinct.cols(10)
#' plot(1:22,pch=19,col=get.distinct.cols(22))
get.distinct.cols <- function(n=22)
{
  distinct.cols <- c(38,6,23,11,12,26,30,31,94,134,100,47,139,53,58,68,116,128,172,142,367,656,77)
  # reorder so similar colours aren't adjacent.
  distinct.cols <- distinct.cols[c(16,10,5,1,21,8,13,18,7,11,3,20,22,14,2,6,19,4,17,12,9,15)]
  colz <- colors()[distinct.cols[1:min(n,length(distinct.cols))]]
  if(n>length(distinct.cols)) { warning(paste(n,"requested, but only",length(distinct.cols),"colors were available")) }
  return(colz)
}




#' Unlist a list, starting only from a set depth.
#' 
#' Allows unlisting preserving the top levels of a list.
#' Can specify the number of list depth levels to skip 
#' before running unlist()
#'
#' @param obj the list to unlist
#' @param depth skip to what layer of the list before unlisting; eg.
#'  the base unlist() function would correspond to depth=0
#' @return returns vectors of strings of char, lengths X
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' complex.list <- list(1,1:3,list(2,2:4,list(3,3:4,list(10))),list(4,5:7,list(3)))
#' Unlist(complex.list,0) # equivalent to unlist()
#' Unlist(complex.list,1) # unlist within the top level lists
#' Unlist(complex.list,2) # unlist within the second level lists
#' Unlist(complex.list,10) # once depth >= list-depth, no difference
#' unlist(complex.list,recursive=FALSE) # not the same as any of the above
Unlist <- function(obj,depth=1) {
  if(depth==0) { return(unlist(obj)) }
  if(length(obj)>0) {
    for(cc in 1:length(obj)) {
      if(is.list(obj)) {
        if(is.list(obj[[cc]])) {
          if(depth<=1) {
            names(obj[[cc]]) <- NULL 
            val <- unlist(obj[[cc]])
            if(is.null(val)) {
              obj[cc] <- list(NULL)
            } else {
              obj[[cc]] <- val            
            }
          } else {
            val <- Unlist(obj[[cc]],depth=depth-1)
            if(is.null(val)) {
              obj[cc] <- list(NULL)
            } else {
              obj[[cc]] <- val            
            }
          }
        }
      }
    }
    return(obj)
  } else {
    return(obj) 
  }
}




#' Convert a numeric vector to Z-scores.
#' 
#' Transform a vector to z scores by subtracting its mean
#'  and dividing by its standard deviation
#'
#' @param X numeric vector to standardize
#' @return vector of the same length in standardised form
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' x1 <- rnorm(10,100,15); x2 <- sample(10)
#' print(x1) ;  standardize(x1)
#' print(x2) ;  standardize(x2)
standardize <- function(X)
{
  if(!is.numeric(X)) { stop("x must be numeric") }
  if(length(X)>1) {
    u <- mean(X,na.rm=T)
    s <- sd(X,na.rm=T)
    return((X-u)/s)
  } else {
    warning("X should have length>1")
    return(X)
  }
}



#' Find data thresholds corresponding to percentiles
#' 
#' Finds the top and bottom bounds corresponding to percentile
#' 'pc' of the data 'dat'.
#'
#' @param dat numeric vector of data
#' @param pc the percentile to seek, c(pc, 1-pc)
#' @return returns the upper and lower threshold
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' pctile(rnorm(100000),.025)
#' pctile(sample(100),.9)
pctile <- function(dat,pc=0.01)
{
  # get top/bottom bounds for percentile 'pc' of data 'dat'
  if(!is.numeric(dat)) { stop("dat must be numeric") }
  rr <- rank(dat,na.last=NA)
  tpbt <- round(c(pc,1-pc)*max(rr,na.rm=T))
  ord <- sort(narm(dat))
  if(tpbt[1]==0) { tpbt[1] <- 1 }
  pcts <- ord[tpbt]
  return(pcts)
}




#' Check whether a given system command is installed (e.g, bash)
#' 
#' Tests whether a command is installed and callable by system().
#' Will return a warning if run on windows when linux.more=TRUE
#'
#' @param cmd character vector of commands to test
#' @param linux.mode logical, alternate way of command testing that only works on linux and
#'  mac OS X, to turn this on, set to TRUE.
#' @return returns true or false for each command in 'cmd'
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' check.linux.install("R") # should be standard
#' check.linux.install(c("perl","sed","fake-cmd"))
check.linux.install <- function(cmd=c("plink","perl","sed"),linux.mode=FALSE) {
  # define function to test any OS
  sys.test <- function(x) {
    X <- Sys.which(x); out <- T
    if(is.na(X)) { out <- F } else { if(is.null(X)) { out <- F } else { if(X=="") { out <- F } } } 
    if(out & !is.character(X)) { out <- F }
    return(out)
  }
  if(linux.mode) {
    if(tolower(.Platform$OS.type)=="windows") {
      warning("function only works on OS X and Linux") ; return(F) 
    } else {
      anz <- character(length(cmd))
      for (dd in 1:length(cmd)) {
        anz[dd] <- system(paste("if hash",cmd[dd],"2>/dev/null; then echo 'yes'; else echo 'no'; fi"),intern=T)
      }
      out <- (anz=="yes")
    }
  } else {
    out <- sapply(cmd,sys.test);
  }
  if(any(!out)) { warning(paste("command",if(length(cmd[!out])>1){ "s" }," '",paste(cmd[!out],collapse="', '"),"' not installed",sep="")) }
  names(out) <- cmd
  return(out)
}


#internal
# why was 'utils' necessary???
#head2 <- function(X,...) { if(length(dim(X))==2) { prv.large(X,...,warn=F) } else { print(utils::head(X,...)) } }
head2 <- function(X,...) { if(length(dim(X))==2) { prv.large(X,...,warn=F) } else { print(head(X,...)) } }



#' A good way to preview large lists.
#' 
#' An alternative to head(list) which allows limiting of large list 
#'  components in the console display
#'
#' @param x a list to preview
#' @param n The number of values to display for the deepest nodes
#'  of the list
#' @param skip number of first level elements to display before skipping
#'  the remainder
#' @param skip2 number of subsequent level elements to display before 
#'  skipping the remainder
#' @param ind indent character for first level elements
#' @param ind2 indent character for subsequent level elements
#' @return prints truncated preview of a large list
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' sub1 <- list(list(1:100),list(2:101),list(101:200),list(201:300),list(301:400))
#' big.list <- list(sub1,sub1,sub1,sub1,sub1,sub1)
#' headl(sub1)
#' headl(big.list,skip=2)
headl <- function (x, n = 6, skip = 20, skip2 = 10, ind = "", ind2 = "  ") 
{
  if (!is.list(x)) {
    warning("not a list")
    return(NULL)
  }
  ll <- length(x)
  if (ll > 0) {
    ind.new <- paste(ind, ind2, sep = "")
    if (ll > skip) {
      ll <- skip
      skipped <- T
    }
    else {
      skipped <- F
    }
    for (cc in 1:ll) {
      if (!is.null(names(x))) {
        cat(ind, "$", names(x)[cc], ":\n", sep = "")
      }
      else {
        (cat("[[", cc, "]]\n", sep = ""))
      }
      if (is(x[[cc]], "list")) {
        headl(x[[cc]], n, ind = ind.new, ind2 = ind2, 
              skip = skip2, skip2 = skip2)
      }
      else {
        cat(ind, sep = "")
        head2(x[[cc]], n)
      }
    }
    if (skipped) {
      cat(ind, "... skipped ", length(x) - skip, " ...\n", 
          sep = "")
    }
  }
  else {
    cat(ind, "<empty>", "\n", sep = "")
  }
}




#' Return an object with missing values removed.
#' 
#'
#' Convenience function, removes NAs from most standard objects.
#' Uses function na.exclude for matrices and dataframes. 
#' Main difference to na.exlude is that it simply performs the 
#' transformation, without adding attributes
#' For unknown types, leaves unchanged with a warning.
#'
#' @param X The object to remove NAs, any vector, matrix or data.frame
#' @return Vector minus NA's, or the matrix/data.frame minus NA rows.
#' If it's a character vector then values of "NA" will also be excluded
#' in addition to values = NA, so be careful if "NA" is a valid value
#' of your character vector. Note that "NA" values occur when 'paste(...,NA,...)' is
#' applied to a vector of any type, whereas 'as.character(...,NA,...)'
#' avoids this.
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' narm(c(1,2,4,NA,5))
#' DF <- data.frame(x = c(1, 2, 3), y = c(0, 10, NA))
#' DF; narm(DF)
#' # if a list, will only completely remove NA from the lowest levels
#' # empty places will be left at top levels
#' print(narm(list(1,2,3,NA,list(1,2,3,NA))))
narm <- function(X) { 
  if(is.data.frame(X) | is.matrix(X)) {
    X <- na.exclude(X)
    attr(X,"na.action") <- NULL
    return(X)
  } else {
    if(is.vector(X)) {
      if(is.list(X)) {
        if(any(sapply(X,length)>1)) {
          X <- lapply(X,narm) 
          return(X)
        } else {
          return(X[!is.na(X)])
        }
      } else {
        if(is.character(X)) {
          ## paste(NA) = "NA", whereas as.character(NA) = NA , this fixes the "NA"'s ##
          out <- X[!is.na(X)]
          out <- out[out!="NA"]
          return(out)
        } else {
          return(X[!is.na(X)])            
        }
      }
    } else {
      warning("unsupported type, X unchanged"); return(X)
    }
  }
}


#' Return a string with each first letter of each word in upper case.
#' 
#' @param txt a character string
#' @param strict whether to force non-leading letters to lowercase
#' @return Vector minus NA's, or the matrix/data.frame minus NA rows
#' @export 
#' @author via R Core
#' @examples
#' toheader(c("using AIC for model selection"))
#' toheader(c("using AIC", "for MODEL selection"), strict=TRUE)
toheader <- function(txt, strict = FALSE) {
  if(!is.character(txt)) { stop("text must be character()") }
  cap <- function(txt) { 
    txt.m <- substring(txt,2); if(strict) { txt.m <- tolower(txt.m) } 
    paste(toupper(substring(txt,1,1)),txt.m,sep = "", collapse = " " ) 
  }
  sapply(strsplit(txt, split = " "), cap, USE.NAMES = !is.null(names(txt)))
}



#' Print heading text with a border.
#'
#' Makes highly visible headings, can separately horizontal, 
#' vertical and corner characters
#'
#' @param txt The text to display in the centre
#' @param h the ascii character to use on the horizontal sections of
#'  the border, and used for v,corner too if not specified separately
#' @param v the character to use on vertical sections of the border
#' @param corner the character to use on corner sections of the border
#' @param align alignment of the writing, when there are multiple lines,
#'  e.g, "right", "left", "centre"/"center"
#' @return returns nothing, simply prints the heading to the console
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' Header("Section 1")
#' Header("Section 1",h="-",v="|",corner="*")
#' Header(c("SPACE","The final frontier"))
#' Header(c("MY SCRIPT","Part 1"),align="left",h=".")
Header <- function(txt,h="=",v=h,corner=h,align="center") {
  ## Make a heading with a box for text (can be multiple lines) optional horiz/vert/corner symbols
  if(is.numeric(txt)) { txt <- paste(txt) }
  if(!is.character(txt)) { stop("txt must be character()") }
  nC <- nchar(txt); align <- tolower(align); if(align!="right" & align!="left") { align <- "center" }
  v <- substr(v,1,1); h <- substr(h,1,1); corner <- substr(corner,1,1)
  extend <- function(X,L,align) { 
    nn <- (L-nchar(X))
    switch(align,right=paste(spc(nn),X,sep=""),
           left=paste(X,spc(nn),sep=""),
           center=paste(spc(floor(nn/2)),X,spc(ceiling(nn/2)),sep="")) }
  mC <- max(nC)
  txt <- extend(txt,mC,align)
  aline <- c(corner,rep(h,mC+2),corner)
  cat("\n",aline,"\n",sep="")
  cat(paste(v," ",txt," ",v,"\n",sep=""),sep="")
  cat(aline,"\n\n",sep="")
}



#' Print a character a specified number of times.
#' 
#' Returns 'char' X_i number of times for each element i of X.
#' Useful for padding for alignment purposes.
#'
#' @param X numeric vector of number of repeats
#' @param char The character to repeat (longer will be shortened)
#' @return returns vectors of strings of char, lengths X
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{rmv.spc}}
#' @examples
#' cat(paste(spc(9),"123\n"))
#' cat(paste(spc(8),"1234\n"))
#' spc(c(1:5),".")
spc <- function(X,char=" ") { 
  if(!is.numeric(X)) { stop("X must be numeric") }
  ch <- substr(paste(char)[1],1,1)
  lX <- length(X); out <- rep("",lX)
  for(j in 1:lX) {
    if(X[j]>0) { out[j] <- paste(rep(ch,times=X[j]),collapse="") }
  } 
  return(out) 
}


#' Remove leading and trailing spaces (or other character).
#'
#' @param str character vector, may containing leading or trailing chars
#' @param before logical, whether to remove leading spaces
#' @param after logical, whether to remove trailing spaces
#' @param char an alternative character to be removed instead of spaces
#' @return returns vectors without the leading/trailing characters
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{spc}}
#' @examples
#' rmv.spc("  mid sentence  ")
#' rmv.spc("0012300",after=FALSE,char="0")
#' rmv.spc(" change nothing ",after=FALSE,before=FALSE)
rmv.spc <- function(str,before=TRUE,after=TRUE,char=" ") {
  # remove spaces at start and end of string
  if(!is.character(str)) { warning("not a character() type") ; return(str) }
  ch <- substr(paste(char)[1],1,1)
  kk <- (length(str))
  if(kk<1) { return(str) }
  for (cc in 1:kk) {
    if(before){
      while(substr(str[cc],1,1)==ch) {
        if(nchar(str[cc])>1) {
          str[cc] <- substr(str[cc],2,nchar(str[cc])) 
        } else {
          str[cc] <- gsub(ch,"",str[cc])
        }
      }
    }
    if(after) {
      while(substr(str[cc],nchar(str[cc]),nchar(str[cc]))==ch) {
        if(nchar(str[cc])>1) {
          str[cc] <- substr(str[cc],1,nchar(str[cc])-1)
        } else {
          str[cc] <- gsub(ch,"",str[cc])
        }
      }
    }
  }
  return(str)
}


#' Estimate the memory required for an object.
#'
#' Can enter an existing object or just the dimensions or total length of a proposed object.
#' The estimate is based on the object being of numeric type. Integers use half the space
#' of numeric, raw() use 1/8th of the space. Factors and characters can vary, although
#' factors will always use less than numeric, and character variables may easily use up
#' to twice as much depending on the length [nchar()] of each element.
#'
#' @param dat either a vector/matrix/dataframe object, or else up to 10 dimensions of such an
#' object, or a potential object, i.e; c(nrow,ncol). If entering an object directly,
#' you can leave out the 'integer' and 'raw' arguments as these will be detected from
#' the object type. Any set of dimensions >10 will be assumed to be a vector, so
#' if you have such an object, better to submit the total product [base::prod()].
#' @param integer if the object or potential object is integer or logical type,
#' set this argument to TRUE, if this is TRUE, the parameter 'RAW' will
#' be ignored; integer and logical types use 1/2 of the memory of numeric types
#' @param raw if the object or potential object is of 'raw' type,
#' set this argument to TRUE, note that if 'integer' is TRUE, this parameter 'RAW' will
#' be ignored; raw types use 1/8 of the memory of numeric types
#' @param unit the storage units to use for the result, ie, "gb", "mb","kb", "b" for
#' gigabytes, megabytes, kilobytes, or bytes respectively.
#' @param add.unit logical, whether to append the unit being used to the result,
#' making the result character type instead of numeric.
#' @return returns the minimum memory requirement to store and object of the specified
#' size, as a numeric scalar, in gigabytes (default) or else using the units specified by 'unit',
#' and if add.unit = TRUE, then the result will be character type instead of numeric, with
#' the units appended.
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' myMatrix <- matrix(rnorm(100),nrow=10)
#' myVec <- sample(1:1000)
#' estimate.memory(myMatrix,unit="bytes") # enter a matrix object
#' estimate.memory(myVec,unit="kb" ,add.unit=TRUE) # enter a vector object
#' estimate.memory(c(10,10,10,10,10),unit="kb") # 5 dimensional array
#' estimate.memory(c(10^6,10^4), add.unit=TRUE) # large matrix
#' estimate.memory(5.4*10^8, add.unit=TRUE)  # entering argument as # total cells, rather than dims
#' estimate.memory(5.4*10^8, integer=TRUE, add.unit=TRUE)
#' estimate.memory(5.4*10^8, raw=TRUE, add.unit=TRUE)
#' estimate.memory(5.4*10^8, TRUE, TRUE, add.unit=TRUE) #  'integer' overrides 'raw'
estimate.memory <- function(dat, integer=FALSE, raw=FALSE, unit=c("gb","mb","kb","b"), add.unit=FALSE)
{
  # based on a numeric object, estimate the minimum memory requirement
  cells.per.gb <- 2^27  # size of double() resulting in ~1GB of memory use by R 2.15
  divisor <- make.divisor(unit,"unit")
  multiplier <- 10^9/divisor
  if(!is.null(dim(dat))) { 
    memory.estimate <- as.numeric(object.size(dat))
    memory.estimate <- memory.estimate/divisor
    if(add.unit) { memory.estimate <- paste(memory.estimate,unit[1]) }
    return(memory.estimate)
  } else { dimz <- dat }
  dimz <- narm(dimz)
  if(length(dimz)==1) { dimz[2] <- 1 }
  if(length(dimz)>1 & length(dimz)<11 & is.numeric(dimz)) {
    total.size <- as.double(1)
    for(cc in 1:length(dimz)) { total.size <- as.double(total.size*as.double(dimz[cc])) }
    memory.estimate <- as.double(as.double(total.size)/cells.per.gb)
    memory.estimate <- memory.estimate*multiplier
    if(integer) { memory.estimate <- memory.estimate/2 } else { if(raw) { memory.estimate <- memory.estimate/8 } }
    if(add.unit) { memory.estimate <- paste(round(memory.estimate,6),unit[1]) }
    return(memory.estimate)
  } else {
    # guessing this is a vector
    if(!is.list(dimz) & is.vector(dimz)) { 
      LL <- length(dimz) 
      return(estimate.memory(LL, integer=integer, raw=raw, unit=unit, add.unit=add.unit))
    } else {
      warning("tried to estimate memory for object which is neither a vector, pair of dimension sizes or a dataframe/matrix") 
    }
  }
}


#internal
make.divisor <- function(unit=c("kb","mb","gb","b"), par.name="unit") {
  valid.units <- c("k","m","g","b")
  unit <- tolower(unit[1]);
  unit <- substr(unit,1,1)
  if(!unit %in% valid.units) { warning("invalid entry to ",par.name," defaulting to kilobytes") ; unit <- "k" }
  divisor <- switch(unit,k=1000,m=10^6, g=10^9, b=1)
  return(divisor)
}


#' Summary of RAM footprint for all R objects in the current session.
#' Not my function, but taken from an R-Help response by Elizabeth Purdom,
#' at Berkeley. Simply applies the function 'object.size' to the objects
#' in ls(). Also very similar to an example in the 'Help' for the 
#' utils::object.size() function.
#' @param unit default is to display "kb", but you can also choose
#' "b"=bytes, "mb"= megabyte, or "gb" = gigabytes. Only the first
#' letter is used, and is not case sensitive, so enter units how you
#' like.
#' @return a list of object names with memory usage in bytes
#' @export
#' @examples
#' memory.summary() # shows memory used by all objects in the current session in kb
#' memory.summary("mb") # change units to megabytes
memory.summary <- function(unit=c("kb","mb","gb","b")) {
  divisor <- make.divisor(unit,"unit")
  out <- sapply(ls(envir=parent.frame(n = 1)),function(x){object.size(get(x))/divisor})
  if(is.atomic(out)) {
    out <- sort(out)
  } 
  return(out)
}

#' Wait for a period of time.
#' 
#' Waits a number of hours minutes or seconds (doing nothing).
#' Note that this 'waiting' will use 100% of 1 cpu.
#'
#' @param dur waiting time
#' @param unit time units h/m/s, seconds are the default
#' @param silent print text showing that waiting is in progress
#' @return no return value
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' wait(.25,silent=FALSE) # wait 0.25 seconds
#' wait(0.005, "m")
#' wait(0.0001, "Hours", silent=FALSE)
wait <- function(dur,unit="s",silent=TRUE) {
  ## do nothing for a period of time
  if(!is.numeric(dur)) { stop("dur must be a number") }
  if(!is.logical(silent)) { silent <- F }
  unit <- tolower(substr(paste(unit),1,1))
  jj <- proc.time()[3]; mm <- 1
  if(unit=="s") { mm <- 1 }
  if(unit=="m") { mm <- 60 }
  if(unit=="h") { mm <- 3600 }
  if(!silent) { timz <- c("hour","minute","second");
                cat("waiting ",dur," ",timz[grep(unit,timz)],"s...",sep="") }
  while((proc.time()[3]-jj)< (mm*dur)) { NULL  }
  if(!silent) { cat("done\n") }
}


#' Times an expression, with breakdown of time spent in each function
#' 
#' A wrapper for the proftools package Rprof() function.
#' It is to Rprof() as system.time() is to proc.time() (base)
#' Useful for identifying which functions are taking the
#' most time. This procedure will return an error unless
#' expr takes more than ~0.1 seconds to evaluate. I 
#' could not see any simple way to avoid this limitation. Occassionally
#' other errors are produced for no apparent reason which are due
#' to issues within the proftools package that are out of my
#' control.
#' 
#' @param expr an expression, must take at least 1 second (roughly)
#' @param suppressResult logical, if true, will return timing 
#'   information rather than the result of expr
#' @param total.time to sort by total.time, else by self.time
#' @return returns matrix where rows are function names, and columns
#'  are self.time and total.time. total.time is total time spent 
#'  in that function, including function calls made by that function.
#'  self.time doesn't count other functions within a function
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' # this function writes and removes a temporary file
#' # run only if ok to do this in your temporary directory
#' #not run# timeit(wait(0.1,"s") ,total.time=TRUE)
#' #not run# timeit(wait(0.1,"s") ,total.time=FALSE)
timeit <- function(expr,suppressResult=F,total.time=TRUE) {
  # function to measure in detail which function calls take the most time
  # during the evaluation of an expression. NB: will error with use of a trivial/instant expression
  tf <- cat.path(tempdir(),"Rproftemp.out")
  Rprof(tf)
  #do the stuff
  result <- { expr }
  Rprof()
  rd <- readProfileData(tf)
  tab <- flatProfile(rd,F)
  if(total.time) { col <- 5 } else { col <- 3 }
  summary <- head(tab[rev(order(tab[,col])),],30)[,c(3,5)]
  unlink(tf)
  if(suppressResult) {
    return(summary)
  } else {
    print(summary)
    return(result)
  }
}


#' Creates a progess bar within a loop
#' 
#' Only requires a single line within a loop to run, in contrast
#' with the built-in tracker which requires a line to initialise,
#' and a line to close. Also has option to backup objects during long loops.
#' Ideal for a loop with a counter such as a for loop.
#' Tracks progress as either percentage of time remaining or
#' by intermittently displaying the estimated number of minutes to go
#'  
#' @param cc integer, current value of the loop counter
#' @param max integer, final value of the loop counter
#' @param st.time 'start time' when using 'time to go' mode, taken 
#'  from a call to proc.time()
#' @param sav.obj optionally an object to backup during the course of 
#'  a very long loop, to restore in the event of a crash.
#' @param sav.fn the file name to save 'save.obj'
#' @param sav.freq how often to update 'sav.obj' to file, in terms of 
#'  percentage of run-time
#' @param unit time units h/m/s if using 'time to go' mode
#' @return returns nothing, simply prints progress to the console
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' # simple example with a for-loop
#' max <- 100; for (cc in 1:max) { loop.tracker(cc,max); wait(0.004,"s") }
#' #example using the 'time to go' with a while loop
#' cc <- 0; max <- 10; start <- proc.time()
#' while(cc < max) { cc <- cc + 1; wait(0.05,"s"); loop.tracker(cc,max,start,unit="s") }
#' # example with saving an object, and restoring after a crash
#' X <- matrix(rnorm(5000),nrow=50); max <- nrow(X); sums <- numeric(max)
#' for (cc in 1:max) { 
#'   sums[cc] <- sum(X[cc,])
#'   wait(.05) # just so this trivial loop doesn't finish so quickly
#'   loop.tracker(cc,max, sav.obj=sums, sav.fn="temp.rda", sav.freq=5);
#'   if(cc==29) { warning("faked a crash at iteration 29!"); rm(sums); break }
#' }
#' cat("\nloaded latest backup from iteration 28:",paste(load("temp.rda")),"\n")
#' print(sav.obj); unlink("temp.rda")
loop.tracker <- function(cc, max, st.time=NULL,  
                         sav.obj=NULL, sav.fn=NA,
                         sav.freq=10, unit=c("m","s","h")[1])
{
  # insert this function into any loop with a counter and it will provide
  # a progress tracker either as a percentage or estimate of time remaining
  ## cc is loop counter, max is the last value the loop should take
  ## st.time is the result of a call to proc.time().
  cc <- round(as.numeric(cc)); max <- round(as.numeric(max))
  freq <- round(max/min(max,50)); 
  if(cc<1) { return() }
  if(cc>max) { cc <- max; warning("count passed to loop.tracker() exceed 'max'") }
  if(cc==1)
  {
    if(is.null(st.time)) {
      scale <- "0         25%         50%         75%         100%\n"
      cat(scale)
    } else {
      cat("Processing: time left (",unit,"):\n",sep="")
    }
  }
  if (cc %% freq == 0) {
    if(is.null(st.time))
    {
      intv <- diff(round(seq(from=1,to=51,length.out=(max/freq))))[cc %/% freq]
      if(!is.na(intv)) { if(intv>0) { cat(rep(".",intv),sep="") } } else {
        if(max==1) { cat(rep(".",50),sep="") }
      }
    } else {
      time.now <- proc.time()[3]-st.time[3]; time.per <- time.now/cc
      tm.u <- switch(unit,m=60,s=1,h=3600)
      to.go <- round(((max-cc)*time.per/tm.u))
      cat(to.go,"..",sep="") 
    }
    if((cc+freq)>max) { cat("\n") }
    ## save as we go - in case of crash
    if(abs(sav.freq)<1) { sav.freq <- abs(sav.freq)*100 } # allow decimal or integer percentage
    sav.freq <- round(max(1,min(50,(sav.freq/2))))
    if ((cc %% (sav.freq*freq)) == 0)
    {
      if(!is.null(sav.obj) & !is.na(sav.fn) & ((max-cc)>1) ) {
        save(sav.obj,file=sav.fn)
      } 
    }
  }
}




#' Make an ascii histogram in the console.
#' 
#' Uses a call to base::hist(...) and uses the densities to make a
#' a text histogram in the console
#' Particularly useful when working in the terminal without graphics.
#'
#' @param X numeric vector of data
#' @param range optional sub-range of X to test; c(low,high)
#' @param ... additional arguments passed to base::hist()
#' @return outputs an ascii histogram to the console
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' textogram(runif(100000))
#' textogram(rnorm(10000),range=c(-3,3))
textogram <- function(X,range=NA,...)
{
  # print a text based histogram of data, extra args may be passed to hist()
  if(!is.numeric(X)) { warning("X must be numeric") ; return(NULL) }
  if(all(!is.na(range))) { X <- X[X>=min(range) & X<=max(range)] }
  hdat <- hist(X,plot=F,...)
  dens <- round(100*hdat$density/sum(hdat$density))
  if(max(dens)>90) {
    cat(" using halved %'s as maximum freq is too big for terminal\n")
    dens <- dens/2
  }
  label <- pad.left(hdat$breaks,char=" ")
  for (cc in 1:length(dens))
  {
    cat(label[cc]," ",paste(rep(".",times=dens[cc]),collapse=""),"\n")
  }
}




#' Print a vector with appropriate padding so each has equal char length.
#' 
#' @param X vector of data to pad to equal length
#' @param char character to pad with, space is default, but zero might
#'  be a desirable choice for padding numbers
#' @param numdigits if using numeric data, the number of digits to keep
#' @return returns the vector in character format with equal nchar()
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' pad.left(1:10)
#' phone.numbers <- c("07429719234","7876345123","7123543765")
#' pad.left(phone.numbers,"0")
#' pad.left(rnorm(10),numdigits=3)
pad.left <- function(X, char=" ", numdigits=NA)
{
  pad <- function(X,L,char=" ") { paste(spc(L-nchar(X),char=char),X,sep="") }
  if (!is.na(numdigits)) { X <- round(X,numdigits)}
  max.len <- max(nchar(X))
  padn <- pad(X,L=max.len,char=char)
  return(padn)
}


#' Do everything possible to load an R package.
#' 
#' Like 'require()' except it will attempt to install a package if
#' necessary, and will also deal automatically with bioconductor
#' packages too. Useful if you wish to share code with people who
#' may not have the same libraries as you, you can include a call to
#' this function which will simply load the library if present, or
#' else install, then load, if they do not have it.
#'
#' @param pcknms list of packages to load/install, shouldn't mix 
#'  bioconductor/CRAN in one call
#' @param bioC whether the listed packages are from bioconductor
#' @param reload indicates to reload the package even if loaded
#' @param avail when bioC=FALSE, see whether pcknms are in the list 
#'  of available CRAN packages
#' @param ask whether to get the user's permission to install a
#'  required package, or just go ahead and do it
#' @param quietly passed to library/require, display installation
#'  text or not
#' @return nothing, simply loads the packages specified if possible
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' # not run : run if you are ok to install/already have these packages
#' # must.use.package(c("MASS","nlme","lme4"),ask=FALSE)
#' # must.use.package("limma",bioC=TRUE)
#' # search() # show packages have loaded, then detach them again:
#' # sapply(paste("package",c("limma","MASS","nlme","lme4"),sep=":"),detach,character.only=TRUE)
must.use.package <- function(pcknms,bioC=FALSE,ask=FALSE,reload=FALSE,avail=FALSE,quietly=FALSE)  
{
  ## like 'base::library()' but can easily work for bioconductor packages, and will
  # automatically attempt to install any required packages not found
  # reload is for when there might a conflict, so will detach and reload packages 
  # to force their version of a function with a duplicate name
  if(!bioC) { 
    repos <- getOption("repos")
    if(any(repos=="@CRAN@" | repos=="")) { repos <- getRepositories(1) }
    if(is.null(repos) | is.na(repos)) { repos <- getRepositories(1) }
    if(avail) {
      goty <- getOption("pkgType"); 
      av.pk <- available.packages(type=goty,
          contrib.url(repos=repos, type=goty)) 
    }
  }
  for (cc in 1:length(pcknms))
  {
    nxt.pck <- pcknms[cc]
    if(!bioC & avail) { if(!nxt.pck %in% av.pk) { 
      warning(nxt.pck,
              " is not in the list of CRAN packages for the current version of R. ",
              "Either it has not compiled successfully for this version, or the name (",
              nxt.pck,") is wrong") } }
    if(reload) {
      fp <- paste("package:",nxt.pck,sep="")
      if(fp %in% search())  { 
        detach(name=fp,character.only=T)
      }
    }
    checklib <- function(package,character.only=FALSE,warn.conflicts=TRUE,quietly=FALSE) { 
      do.call("require",args=list(package=package,character.only=character.only,
                             warn.conflicts=warn.conflicts,quietly=quietly)) 
    }
    got1 <- suppressWarnings(checklib(nxt.pck,character.only=TRUE,warn.conflicts=FALSE))
    if(!got1) {
      if(ask) {
        # ask whether to install a package
        ans <- select.list(c("yes","no"),"yes",FALSE,paste("ok to install",nxt.pck," (required)?"))
      } else { 
        ans <- "yes" 
      }
      if(ans=="yes") {
        if(bioC) {
          source("http://bioconductor.org/biocLite.R",local=TRUE) # biocLite() should now be replaced
          if(!exists("biocLite",mode="function")) {
            biocLite <- function(x,...) { 
              cat("please load biocLite function from http://bioconductor.org/biocLite.R and install",nxt.pck,"manually") 
            } 
          }
          repos <- getOption("repos")
          if(any(repos=="@CRAN@" | repos=="")) { repos <- getRepositories(1) }
          if(is.null(repos) | is.na(repos)) { repos <- getRepositories(1) }
          biocLite(nxt.pck,siteRepos=repos)
          suppressWarnings(worked <- checklib(nxt.pck,character.only=TRUE,warn.conflicts=FALSE,quietly=quietly))
        } else {
          install.packages(nxt.pck,repos=repos,dependencies=TRUE); 
          suppressWarnings(worked <- checklib(nxt.pck,character.only=TRUE,warn.conflicts=FALSE,quietly=quietly)) 
        }
        if(!worked) { warning("automated installation of required package: ",nxt.pck," failed") }
      } else {
        warning("please manually install package ",nxt.pck," to continue")
      }
    } 
  }
}


#' Search all CRAN packages for those containing keyword(s).
#' 
#' Can be useful for trying to find new packages for a particular
#' purpose. No need for these packages to be installed or loaded.
#' Further searching can be done using utils::RSiteSearch()
#'
#' @param txt text to search for, a character vector, not case-sensitive
#' @param repos repository(s) (CRAN mirror) to use, "" defaults to getOption("repos")
#' @param all.repos logical, if TRUE, then use all available repositories from getRepositories()
#' @return list of hits for each keyword (txt)
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' repos <- "http://cran.ma.imperial.ac.uk/" # OR: repos <- getOption("repos")
#' # setRepositories(ind=1:2) # for the session will by default search bioconductor packages too
#' search.cran("useful",repos)
#' search.cran(c("hmm","markov","hidden"),repos=repos)
#' require(BiocInstaller)
#' search.cran(c("snpStats","genoset","limma"),all.repos=TRUE)
search.cran <- function(txt,repos="",all.repos=FALSE) {
  goty <- getOption("pkgType"); 
  if(all.repos) {
    repos <- getRepositories() # use all available
  }
  if(all(repos=="")) { 
    repos <- getOption("repos") 
  }
  if(any(repos=="@CRAN@")) { repos <- "http://cran.ma.imperial.ac.uk/" }
  av.pk <- available.packages(type=goty,
           contrib.url(repos=repos, type=goty))
  if(is.matrix(av.pk)) { 
    if("Package" %in% colnames(av.pk)) {
      av.pk <- av.pk[,"Package"]; dim(av.pk) <- NULL
    } else { av.pk <- av.pk[[1]] }
  } else { 
    warning("lookup did not return table with header 'Package'")
  }
  if(is.character(av.pk) & is.character(txt)) {
    if(!is.null(names(av.pk))) { names(av.pk) <- NULL }
    if(length(txt)>0) {
      out <- vector("list",length(txt)); names(out) <- txt
      for(cc in 1:length(txt)) {
        out[[cc]] <- av.pk[grep(txt[cc],av.pk,ignore.case=T)]
      }
    }
  } else {
    warning("txt must be character, and must be online to search for available.packages()")
  }
  return(out)
}


#' Detect all available R repositories.
#' 
#' In addition to the default CRAN repository, there are other repositories such
#' as R-Forge, Omegahat, and bioConductor (which is split in to software, annotation,
#' experiments and extras). This function allows you to retrieve which are available.
#' This function complements (and takes code from) utils::setRepositories(), which
#' will just set, not return which are available, but see there for more information
#' about how this works. Detecting the available repositories can be useful to precede
#' a call to setRepositories, and allows you to utilise these repositories without
#' calling setRepositories (which is hard to reverse). This function can be used to
#' expand the search space of the function search.cran() to include bioconductor packages.
#' @param ind index, same as for 'setRepositories', if NULL this function returns all available
#' repositories, or if an index, returns a subset.
#' @param table logical, if TRUE, return a table of information, else just return the URLs, 
#' which are the required input for the 'repos' argument for relevant functions, 
#' e.g, available.packages() or search.cran()
#' @return list of repositories with URLS, note that it is the URL that works best for
#' use for passing a value for 'repos' to various functions.
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' require(BiocInstaller)
#' repos <- "http://cran.ma.imperial.ac.uk/" # OR: repos <- getOption("repos")
#' getRepositories(table=TRUE) # shows all available
#' getRepositories(2:5,FALSE) # returns index for all bioconductor repositories (on my system at least)
#' search.cran("genoset",repos=getRepositories(1)) # does not find this bioconductor package on CRAN
#' search.cran("genoset",repos=getRepositories()) # should now, because all repositories are used
getRepositories <- function(ind = NULL,table=FALSE) {
  p <- file.path(Sys.getenv("HOME"), ".R", "repositories")
  if (!file.exists(p)) 
    p <- file.path(R.home("etc"), "repositories")
  a <- tools_read_repositories(p)  ## had to hack this function together as internal ::: to tools package
  pkgType <- getOption("pkgType")
  if (pkgType == "both") 
    pkgType <- .Platform$pkgType
  if (length(grep("^mac\\.binary", pkgType))) 
    pkgType <- "mac.binary"
  thisType <- a[[pkgType]]
  a <- a[thisType, 1L:3L]
  repos <- getOption("repos")
  if ("CRAN" %in% row.names(a) && !is.na(CRAN <- repos["CRAN"])) 
    a["CRAN", "URL"] <- CRAN
  a[(a[["URL"]] %in% repos), "default"] <- TRUE
  new <- !(repos %in% a[["URL"]])
  if (any(new)) {
    aa <- names(repos[new])
    if (is.null(aa)) 
      aa <- rep("", length(repos[new]))
    aa[aa == ""] <- repos[new][aa == ""]
    newa <- data.frame(menu_name = aa, URL = repos[new], 
                       default = TRUE)
    row.names(newa) <- aa
    a <- rbind(a, newa)
  }
  if(is.numeric(ind)) { 
    ind[ind<1] <- 1; ind[ind>nrow(a)] <- nrow(a); ind <- unique(ind)
    a <- (a[ind,]) 
  } 
  if(table){
    return(a)
  } else {
    return(a[,2])
  }
}


# internal
check.bio <- function() {
	if("BiocInstaller" %in% installed.packages()) {
		do.call("require",args=list("BiocInstaller"))
		return(do.call("biocVersion"))
	} else {
		warning("bioconductor does not appear to be installed - this function works better if it is")
		stop("deliberately throw error for 'tryCatch' to catch")
	}
}



# internal function stolen from 'tools'
tools_read_repositories <- function (file) 
{
  # try to replicate the constant 'tools:::.BioC_version_associated_with_R_version'
  get.bioc.version <- function() {
    biocVers <- tryCatch({
      check.bio() # recent BiocInstaller
    }, error=function(...) {         # no / older BiocInstaller
      numeric_version(Sys.getenv("R_BIOC_VERSION", "2.13"))
    })
    return(biocVers)
  }
  tools_expand_BioC_repository_URLs <- function (x) 
  {
    x <- sub("%bm", as.character(getOption("BioC_mirror", "http://www.bioconductor.org")), 
             x, fixed = TRUE)
    sub("%v", as.character(get.bioc.version()), 
        x, fixed = TRUE)
  }
  db <- utils::read.delim(file, header = TRUE, comment.char = "#", 
                          colClasses = c(rep.int("character", 3L), rep.int("logical", 
                                                                           4L)))
  db[, "URL"] <- tools_expand_BioC_repository_URLs(db[, "URL"])
  return(db)
}



#' Find the mode of a vector.
#'
#' The mode is the most common value in a series.
#' This function can return multiple values if there are equally
#' most frequent values, and can also work with non-numeric types.
#'
#' @param x The data to take the mode from. Dimensions and NA's are 
#'  removed if possible, strings, factors, numeric all permitted
#' @param multi Logical, whether to return multiple modes if values
#'  have equal frequency
#' @param warn Logical, whether to give warnings when multiple values
#'  are found (if multi=FALSE)
#' @return The most frequent value, or sorted set of most frequent
#'  values if multi==TRUE and there are more than one. Numeric if x 
#'  is numeric, else as strings
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' Mode(c(1,2,3,3,4,4)) # 2 values are most common, as multi=FALSE, 
#' # selects the last value (after sort)
#' Mode(c(1,2,3,3,4,4),multi=TRUE) # same test with multi=T, 
#' # returns both most frequent
#' Mode(matrix(1:16,ncol=4),warn=TRUE) # takes mode of the entire
#' # matrix treating as a vector, but all values occur once
#' Mode(c("Tom","Dick","Harry"),multi=FALSE,warn=TRUE) # selects last
#' # sorted value, but warns there are multiple modes
#' Mode(c("Tom","Dick","Harry"),multi=TRUE,warn=TRUE) # multi==TRUE so
#' # warning is negated
Mode <- function(x,multi=FALSE,warn=FALSE) {  
  ## mode function that should work for vectors or arrays of numeric, character, factor
  dim(x) <- NULL; x <- narm(x); tt <- 1
  if(length(x)>0) {
    ii <- table(x); ii <- sort(ii); ll <- length(ii)
    if(length(unique(as.numeric(ii)))==1 & length(as.numeric(ii))>1) {
      if(multi) {
        tt <- length(ii)
      } else { if(warn) { warning("all values of x had equal frequency, returning greatest") } }
    } else {
      if(ll>1) {
        if(ii[ll]==ii[ll-1]) { 
          if(multi) {
            tt <- length(which(ii==ii[ll]))
          } else { if(warn) { warning("multiple values of x had the maximum frequency, returning greatest") } }
        }
      }
    }
    result <- tail(names(ii),tt)
    nresult <- suppressWarnings(as.numeric(result))
    if(all(!is.na(nresult))) { result <- nresult }
    return(result) 
  } else {
    warning("no valid values passed to Mode")
    return(NA)
  }
}




#' Create an index file for an R function file
#'
#' Create a html index for an R function file by looking for functions,
#'  add descriptions using comments directly next to the function()
#'  command. Note that if too much code other than well-formatted
#'  functions is in the file then the result is likely not to be
#'  a nicely formatted index.
#' 
#' @param fn an R file containing functions in standard R script
#' @param below whether to search for comment text below or above
#'  the function() calls
#' @param fn.out optional name for the output file, else will be 
#'  based on the name of the input file
#' @param skip.indent whether to skip functions that are indented,
#'  the assumption being they are functions within functions
#' @return creates an html file with name and description of each function
#' @seealso \code{\link{list.functions.in.file}}
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' # not run:  rfile <- file.choose() # choose an R script file with functions
#' # not run:  out <- Rfile.index(rfile,fn.out="temp.htm")
#' # unlink("temp.htm") # run once you've inspected this file in a browser
Rfile.index <- function(fn,below=TRUE,fn.out="out.htm", skip.indent=TRUE) 
{
  # makes html index of each function in a large functions file
  #require(reader)
  grp <- function(what,ins) { grep(what,ins,fixed=T) }
  if(toupper(get.ext(fn))!="R") { warning("expecting 'fn' to be a *.R file, if not, expect failure") }
  #if(is.null(fn.out)) { fn.out <- cat.path(fn=fn,suf="index",ext="htm") }
  neg <- 1; if(!below) { neg <- -1 }
  if(file.exists(fn))  {
    fl <- readLines(fn)
    fn.lines <- unique(c(grp("<- function",fl),grp("<-function",fl)))
    if(skip.indent) { 
      indz <- which(substr(fl[fn.lines],1,1) %in% c(" ","\t"))
      if(length(indz)>0) { fn.lines <- fn.lines[-indz] }
    }
    fl <- rmv.spc(fl)
    indz2 <- grep("#",fl[fn.lines]); 
    indz3 <- grep("function",sapply(strsplit(fl[fn.lines][indz2],split="#"),"[",1))
    if(length(indz3)>0) {  indz2 <- indz2[-indz3] }
    if(length(indz2)>0) {  fn.lines <- fn.lines[-indz2] }
    nfn <- length(fn.lines)
    fn.list <- vector("list",nfn)
    if(nfn<1) { warning("no functions found in R file"); return(NULL) }
    for (kk in 1:nfn) {
      first.ln <- fl[fn.lines[kk]]
      n <- 1; while(substr(first.ln,n,n)!="<" & substr(first.ln,n,n)!=" ") { n <- n+1 }
      fn.nm <- substr(first.ln,1,n-1);# cat(fn.nm,"\n")
      names(fn.list)[kk] <- paste("<p></p><b>",fn.nm,"</b>",sep=""); descr <- c()
      lnn <- fn.lines[kk]; 
      if(below) { while(length(grp("{",fl[lnn]))==0) { lnn <- lnn+neg } }
      #print(fl[lnn])
      lnn <- lnn+neg ; 
      while(length(grp("#",fl[lnn]))>0) { 
        descr <- c(descr,gsub("#","",gsub("#'","",fl[lnn],fixed=T),fixed=T))
        lnn <- lnn+neg 
      }
      if(!below) { descr <- rev(descr) }
      # remove lines after @ characters (for roxygen)
      roxy <- grep("@",descr); if(length(roxy)>0) { 
        descr <- descr[-c(min(roxy):length(descr))] }
      fn.list[[kk]] <- rmv.spc(paste(descr))
    }
  } else {
    warning("could not find function file to index")
    return(NULL)
  }
  fil.cont <- sapply(fn.list,paste,collapse="\n")
  #return(fil.cont)
  write.table(fil.cont,file=fn.out,quote=F,col.names=F)
  cat("wrote html index to: ",fn.out,"\n")
  return(fn.list)
}


#' Show all functions used in an R script file, by package
#'
#' Parses all functions called by an R script and then lists
#' them by package. Wrapper for 'getParseData'. Inspired by
#' 'hrbrmstr', on StackExchange 3/1/2015. May be of great
#' use for those developing a package to help see what 
#' namespace 'importsFrom' calls will be required.
#' @param filename path to an R file containing R code.
#' @param alphabetic logical, whether to list functions alphabetically.
#' If FALSE, will list in order of appearance.
#' @return Returns a list. Parses all functions called by an R script 
#' and then lists them by package. Those from the script itself are listed
#' under '.GlobalEnv' and any functions that may originate
#' from multiple packages have all possibilities listed. Those listed under
#' 'character(0)' are those for which a package could not be found- may be
#' functions within functions, or from packages that aren't loaded.
#' @seealso \code{\link{Rfile.index}}
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' # not run:  rfile <- file.choose() # choose an R script file with functions
#' # not run:  list.functions.in.file(rfile)
list.functions.in.file <- function(filename,alphabetic=TRUE) {
  # from hrbrmstr, StackExchange 3/1/2015
  if(!file.exists(filename)) { stop("couldn't find file ",filename) }
  if(!get.ext(filename)=="R") { warning("expecting *.R file, will try to proceed") }
 # requireNameSpace("dplyr")
  tmp <- getParseData(parse(filename, keep.source=TRUE))
  crit <- quote(token == "SYMBOL_FUNCTION_CALL")
  tmp <- dplyr::filter(tmp, .dots = crit)
  #tmp <- dplyr::filter(tmp,token=="SYMBOL_FUNCTION_CALL")
  tmp <- unique(if(alphabetic) { sort(tmp$text) } else { tmp$text })
  src <- paste(as.vector(sapply(tmp, find)))
  outlist <- tapply(tmp,factor(src),c)
  return(outlist)
}



#' INTERNAL: Get the file extension from a file-name.
get.ext <- function(fn) {
  # get file extension from a filename character string
  if(length(fn)<1) { warning("fn had length of zero"); return(fn) }
  if(all(is.na(fn)) | !is.character(fn)) { stop("fn should not be NA and should be of type character()") }
  strip.file.frags <- function(X) {
    file.segs <- strsplit(X,".",fixed=T)[[1]]
    lss <- length(file.segs)
    if (lss>1) { out <- paste(file.segs[lss]) } else { out <- "" }
    return(out)
  }
  return(sapply(fn,strip.file.frags))
}




#' Tidy display function for matrix objects
#'
#' This function prints the first and last columns and rows of a matrix, and
#' more, if desired. Allows previewing of a matrix without 
#' overloading the console. Most useful when data has row and column names.
#'
#' @param largeMat a matrix
#' @param rows number of rows to display
#' @param cols number of columns to display
#' @param digits number of digits to display for numeric data
#' @param rL row label to describe the row names/numbers, e.g, row number, ID, etc
#' @param rlab label to describe the data rows
#' @param clab label to describe the data columns
#' @param rownums logical, whether to display rownumbers or ignore them
#' @param ret logical, whether to return the result as a formatted object, or just print to console
#' @param warn logical, whether to warn if the object type is not supported
#' @export
#' @examples
#' mat <- matrix(rnorm(1000),nrow=50)
#' rownames(mat) <- paste("ID",1:50,sep="")
#' colnames(mat) <- paste("Var",1:20,sep="")
#' prv.large(mat)
#' prv.large(mat,rows=9,cols=4,digits=1,rlab="samples",clab="variables",rownums=FALSE)
prv.large <- function(largeMat,rows=3,cols=2,digits=4,rL="Row#",
                      rlab="rownames",clab="colnames",rownums=T,ret=FALSE,warn=TRUE) 
{
  # nicely print a large matrix without overloading the output space
  # can return result as lines of text instead of printing to screen (for printing to file)
  # allows customization of row and column labels
  # only worth using with data that has row/col names
  # DEFINE INTERNAL FUNCTIONS #
  pad <- function(X,L) { X<-paste(X); if(is.character(X)) { paste(spc(L-nchar(X)),X,sep="") } else { stop(X) } }
  RND <- function(X,...) { if (is.numeric(X)) { round(X,...) } else { X }}
  #
  if(packages.loaded("bigmemory",cran.check=F)) { TF <- !do.call("is.big.matrix",args=list(largeMat)) } else { TF <- TRUE }
  if(!is.data.frame(largeMat) & !is.matrix(largeMat) & TF ) { 
    if(warn) { warning("unsupported object type, using 'head'") }
    print(head(largeMat))
    return()
  }
  if(length(dim(largeMat))!=2) { stop("expected largeMat to have 2 dimensions") }
  nC <- ncol(largeMat); nR <- nrow(largeMat); 
  if(nC<2 | nR<3) { 
    if(warn) {
      warning("prv.large only works for matrices with dims >= c(3,2), passed to print(head())")
    }
    print(head(largeMat,rows+1)); return(NULL) 
  }
  rows <- min(max(1,rows),nR); cols <- min(max(1,cols),nC)
  cN <- colnames(largeMat); rN <- rownames(largeMat)
  if(is.null(cN)) { cN <- paste(1:ncol(largeMat)); clab <- "col#" }
  if(is.null(rN)) { rN <- paste(1:nrow(largeMat)); rlab <- "row#"; rownums=F }
  rD <- spc(min(2,max(nchar(paste(nR)))),".")
  rnD <- spc(min(4,max(nchar(rN[c(1:rows,nR)]))),".")
  linez <- vector("list",rows+3) #row,col =number of rows,cols to print
  rown <- max(nchar(paste(nR)),nchar(rL))*as.numeric(rownums)
  hdr <- (nchar(cN[c(1:cols,nC)]))
  if(is.numeric(largeMat[1,])) {
    ## assess what the largest numbers are likely to be to adjust header spacing if necessary
    long.nums <- max(max(abs(largeMat[1,]),na.rm=T),max(abs(largeMat[,1]),na.rm=T))
    max.before.dp <- nchar(round(long.nums))+3
  } else { max.before.dp <- 6 }
  hdr[hdr<7] <- 7; hdr[hdr<(digits+max.before.dp)] <- (digits+max.before.dp)
  idln <- max(nchar(rlab),nchar(rN[c(1:rows,nR)]))
  long.text.filt <- T
  if(long.text.filt) {
    # look ahead for lengths of textdata > nchar(colnames()) which wouldn't get picked up elsewhere
    dat.cols <- c(1:cols,nC); varln <- a1 <- b1 <- numeric(length(dat.cols))
    for(cc in 1:length(dat.cols)) {
      a1[cc] <- nchar(cN[dat.cols[cc]]); b1[cc] <- max(nchar(paste(RND(largeMat[c(1:rows,nR),dat.cols[cc]],digits))))
      varln[cc] <- max(a1[cc],b1[cc])
    }
    #prv(a1,b1,dat.cols,hdr,varln)
    hdr[varln>hdr] <- varln[varln>hdr]
  }
  if(!ret) { cat("\n"); cat(spc(rown),spc(idln),clab,"\n") }
  dotz <- "  ...  "; dotzh <- " ..... "; dotzn <- "..."
  # make adjustments if matrix is small enough to display all rows/cols
  if(nC<=cols) { dotz <- dotzh <- "" ; cols <- cols-1 }
  if(nR<=rows) { lstln <- 1 } else {  lstln <- 3 }
  ## make adjustments if not displaying rownumbers
  if(!rownums) {
    lstR <- "" ; rD <- ""; jstr <- rep("",times=rows); rL=""
  } else {
    lstR <- nR; jstr <- paste(1:rows)
  }
  linez[[1]] <- c(pad(rL,rown),pad(rlab,idln),pad(cN[c(1:cols)],hdr[1:cols]),
                  dotzh,pad(cN[nC],tail(hdr,1)))
  for (j in 1:rows) { 
    linez[[j+1]] <- c(pad(jstr[j],rown),pad(rN[j],idln),
                      pad(RND(unlist(largeMat[j,1:cols]),digits),hdr[1:cols]),dotz,
                      pad(RND(largeMat[j,nC],digits),tail(hdr,1)))
  }
  linez[[rows+2]] <- c(pad(rD,rown),pad(rnD,idln),pad(rep(dotzn,times=cols),
                                                      hdr[1:cols]),dotz,pad(dotzn,tail(hdr,1)))
  linez[[rows+3]] <- c(pad(lstR,rown),pad(rN[nR],idln),
                       pad(RND(unlist(largeMat[nR,1:cols]),digits),hdr[1:cols]),
                       dotz,pad(RND(largeMat[nR,nC],digits),tail(hdr,1)))
  if(!ret) {
    for (j in 1:(rows+lstln)) {
      cat(paste(linez[[j]],collapse=" "),"\n")
    }
  } else {
    # remove last two lines if all rows are displayed
    if(lstln==1) { for(ii in 1:2) { linez[[length(linez)]] <- NULL }  }
    return(linez)
  }
}




#' Convert objects as arguments to object names
#' 
#' Equivalent to the base function substitute() but can do any length of arguments instead
#' of just one. Converts the objects in parentheses into text arguments as if they
#' had been entered with double quote strings. The objects must exist and be accessible in
#' the environment the function is called from for the function to work (same as for substitute()).
#' One application for this is to be able to create functions where object arguments can be
#' entered without quotation marks (simpler), or where you want to use the name of the object
#' as well as the data in the object.
#'
#' @param x compulsory, simply the first object in the list, no difference to any further objects
#' @param ... any further objects to return string names for.
#' @return character list of x,... object names
#' @export
#' @seealso \code{\link{prv}}, \code{\link{preview}} 
#' @author Nicholas Cooper 
#' @examples
#' myvar <- list(test=c(1,2,3)); var2 <- "testme"; var3 <- 10:14
#' print(myvar)
#' # single variable case, equivalent to base::substitute()
#' print(substitute(myvar))
#' print(Substitute(myvar))
#' # multi variable case, substitute won't work
#' Substitute(myvar,var2,var3)
#' # prv() is a wrapper for preview() allowing arguments without parentheses
#' # which is achieved internally by passing the arguments to Substitute()
#' preview(c("myvar","var2","var3"))
#' prv(myvar,var2,var3)
Substitute <- function(x=NULL,...) {
  varlist <- list(...); out <- character(1)
  if(length(varlist)>0) { 
    extr <- Substitute(...)
  } else {
    extr <- NULL
  }
  if(!is.null(x)) { out[1] <- paste(substitute(x)) }
  out <- c(out,extr)
  return(out[out!=""])
}


#' Output variable states within functions/loops during testing/debugging
#'
#' Same as preview but no labels command, and input is without quotes
#' and should be plain variable names of existing variables (no indices, args, etc)
#' A versatile function to compactly display most common R objects. Will
#' return the object name, type, dimension, and a compact representation of
#' object contents, for instance using prv.large() to display matrices,
#' so as to not overload the console for large objects. Useful for debugging,
#' can be placed inside loops and functions to track values, dimensions, and data types.
#' Particularly when debugging complex code, the automatic display of the variable name
#' prevents confusion versus using regular print statements.
#' By listing variables to track as character(), provides 'cat()' output 
#' of compact and informative variable state information, e.g, variable name, value,
#' datatype and dimension. Can also specify array or list elements, or custom labels.
#' prv() is the same as preview() except it can take objects without using double quotes
#' and has no 'labels' command (and doesn't need one). If expressions are entered rather
#' than variable names, then prv() will attempt to pass the arguments to preview().
#' prv() assumes that the variable(s) to report originate from the environment calling
#' prv(), and if not found there, then it will search through all accessible environments
#' starting with the global environment, and then will report the first instance found,
#' which in exceptional circumstances (be warned) may not be the instance you intended
#' to retrieve.
#' @param ... series of variable(s) to report, separated by commas, which will trigger
#'  automatic labelling of the variable name
#' @param counts a list of array index values; so if calling during a counting loop, the
#'  value can be reported each iteration, also printing the count index; if the list is
#'  named the name will also appear, e.g, variable[count=1]. This list must be the same
#'  length as the variable list ... , and each element [[i]] must contain as many values
#'  as the original corresponding variable list[i] has dimensions
#' @seealso \code{\link{Dim}}
#' @export
#' @examples
#' # create variables of different types to show output styles #
#' testvar1 <- 193
#' testvar2 <- "Atol"
#' testvar3 <- c(1:10)
#' testvar4 <- matrix(rnorm(100),nrow=25)
#' testvar5 <- list(first="test",second=testvar4,third=100:110)
#' preview("testvar1"); prv(testvar1)
#' prv(testvar1,testvar2,testvar3,testvar4)
#' prv(matrix(rnorm(100),nrow=25)) # expression sent to preview() with no label
#' prv(193) # fails as there are no object names involved
prv <- function(...,counts=NULL) {
  options(warn=2)
  txt <- (tryCatch(varlist <- Substitute(...), error = function(e) e))
  options(warn=0)
  if(is(txt)[1]=="simpleError") { 
    #warning("not a function name")
    varlist <- list(...)
    sapply(varlist,preview,prv.call=TRUE)
    return(NULL)
  }
  return(preview(varlist,labels=NULL,counts=counts,prv.call=TRUE))
}


# internal function for print.large
display.var <- function(val,label,cnts=NULL) {
  if(is(cnts)[1]=="list") {
    ## if vars to debug have a counter, update the value and label with count(s)
    if(is(val)[1]=="list") { 
      for (dd in 1:length(cnts)) {
        val <- val[[ cnts[[dd]] ]] 
        if(!is.null(names(cnts))) { 
          label <- paste(label,"[[",names(cnts)[dd],"=",cnts[[dd]],"]]",sep="") 
        } else {
          label <- paste(label,"[[",cnts[[dd]],"]]",sep="")
        }
      }
    } else {
      #val <- val[cnts[[dd]] ]
      #preview(c("val","cnts"))
      if(length(Dim(val))!=length(cnts)) {
        val <- val ; warning("counts did not match dimensions")
      } else {
        arg.list <- vector("list",1+length(cnts)); arg.list[[1]] <- val
        arg.list[2:(1+length(cnts))] <- cnts
        val <- do.call("[",args=arg.list)
        if(!is.null(names(cnts))) { 
          label <- paste(label,"[",
                         paste(paste(names(cnts),"=",cnts,sep=""),collapse=","),"]",sep="") 
        } else {
          label <- paste(label,"[",paste(cnts,collapse=","),"]",sep="")
        }
      }
    }
  } else {
    #counts not a list
  }
  ## display appropriately according to datatype ##
  typ <- is(val)[1]
  if(is.function(val)) {
    cat(label,": function",sep=""); return(invisible())
  }
  if(packages.loaded("bigmemory",cran.check=F)) {
    if(typ=="big.matrix") {
      if(exists("prv.big.matrix",mode="function")) {
        do.call("prv.big.matrix",args=list(val,name=label))
      } else {
        warning("preview() needs the package bigpca to display a big.matrix object")
      }
      return(invisible())
    }
  }
  dv <- Dim(val)
  if(is.numeric(dv)) { if(all(dv==1)) {
    if(is.vector(val)) {
      cat(label,": ",val," (",typ,", ",paste(Dim(val),collapse="*"),")",sep=""); return(invisible())
    }
  } }
  if(is(val)[1]=="list") {
    cat(label," (",typ,", ",paste(Dim(val),collapse="*"),")\n",sep=""); print(headl(val)); return(invisible())
  } else {
    #print(Dim(val))
    if(!is.null(dim(val))) {
      cat(label," (",typ,", ",paste(Dim(val),collapse="*"),")\n",sep="");
      if(length(dim(val))==2) {
        if(ncol(val)>=2 & nrow(val)>=3) {
          prv.large(val,warn=F)
        } else {
          print(head(val))
          if(nrow(val)>6) {
            # if any part not displayed, then indicate using ...
            cat(if(!is.null(rownames(val))) { "  ...    " } else { "" },rep("  ..  ",ncol(val)),"\n")
          }
        }
      } else {
        print(c(head(val),if(length(val)>6) { (" ...") } else { NULL }))  # e.g, for a table
      }
      return(invisible())
    } else {
      cat(label," (",typ,", ",paste(Dim(val),collapse="*"),") [head]:\n",sep="")
      print(head(val))
      return(invisible())
    }
  }
}


#' Output variable states within functions during testing/debugging
#'
#' A versatile function to compactly display most common R objects. Will
#' return the object name, type, dimension, and a compact representation of
#' object contents, for instance using prv.large() to display matrices,
#' so as to not overload the console for large objects. Useful for debugging,
#' can be placed inside loops and functions to track values, dimensions, and data types.
#' Particularly when debugging complex code, the automatic display of the variable name
#' prevents confusion versus using regular print statements.
#' By listing variables to track as character(), provides 'cat()' output 
#' of compact and informative variable state information, e.g, variable name, value,
#' datatype and dimension. Can also specify array or list elements, or custom labels.
#' prv() is the same as preview() except it can take objects without using double quotes
#' and has no 'labels' command (and doesn't need one).
#' @param varlist character vector, the list of variable(s) to report, which will trigger
#'  automatic labelling of the variable name, otherwise if entered as the variable value (ie.
#'  without quotes, then will by default be displayed as 'unknown variable')
#' @param labels, will label 'unknown variables' (see above) if entered as variables without quotes
#' @param counts a list of array index values; so if calling during a counting loop, the
#'  value can be reported each iteration, also printing the count index; if the list is
#'  named the name will also appear, e.g, variable[count=1]. This list must be the same
#'  length as varlist (and labels if not NULL), and each element [[i]] must contain as many values
#'  as the original corresponding varlist[i] has dimensions. The dimensions must result in a 1x1 scalar
#' @param assume.char usually 'varlist' is a character vector of variable names, but in the case
#'  that it is actually a character variable, using assume.char=TRUE will ensure that it will be assumed
#'  the character variable is the object to preview, rather than a list of variable names. So long
#'  as none of the values are found to be variable names in the global environment. preview() can
#'  also find variables in local environments, and if this is where the target variable lies, it is
#'  best to use assume.char=FALSE, otherwise the search for alternative environments might not happen.
#'  Note that in most cases the automatic detection of the input should understand what you want, regardless
#'  of the value of assume.char.
#' @param prv.call It is recommended to always leave this argument as FALSE when calling preview()
#' directly. If set to TRUE, it will first search 2 generations back for the parent frame, instead 
#' of one, as it will assume that the variable(s) to preview are not directly called by preview(),
#' but through a wrapper for preview, such as prv().
#' @seealso \code{\link{Dim}} 
#' @export
#' @examples
#' # create variables of different types to show output styles #
#' testvar1 <- 193
#' testvar2 <- "Atol"
#' testvar3 <- c(1:10)
#' testvar4 <- matrix(rnorm(100),nrow=25)
#' testvar5 <- list(first="test",second=testvar4,third=100:110)
#' preview("testvar1")
#' preview("testvar4")
#' preview(paste("testvar",1:5,sep=""))
#' preview(testvar1,"myvarname")
#' preview(testvar1)
#' # examples with loops and multiple dimensions / lists
#' for (cc in 1:4) {
#'  for (dd in 1:4) { preview("testvar4",counts=list(cc,dd)) }}
#'
#' for (dd in 1:3) { preview("testvar5",counts=list(dd=dd)) }
preview <- function(varlist,labels=NULL,counts=NULL,assume.char=FALSE, prv.call=FALSE) {
  ## for debugging, simplify code to print vars you are checking
  lab <- varlist
  if(is.character(varlist) & (length(labels)<length(varlist))) {
    if(assume.char | length(varlist)>10) {
      if(!any(varlist %in% ls())) {
        unknown.variable <- varlist
        preview("unknown.variable",labels=labels,counts=counts)
        return(invisible(NULL))
      }
    }
  } 
  if(is.character(varlist)) {
    t1 <- grep("[",varlist,fixed=T)
    t2 <- grep("(",varlist,fixed=T)
    if(length(t1)>0 | length(t2)>0) {
      warning("preview() only works with plain variable names, cannot use an index or function",
              "call containing ['s or ('s. To access object indices use the 'counts' argument")    
      return()
    }
  }
  # test whether 'counts' sublists are all of the same length as varlist, else ignore 'counts'
  if(is.list(counts)) {  if(!all(sapply(counts,length)==length(varlist))) { 
    counts <- NULL } } else { if(length(counts)==length(varlist)) { counts <- list(counts) } else { counts <- NULL } }
  #val <- vector("list",length(lab))
  
  ## if data not entered with a label, or as a string (not including prv() converted calls)
  if(!is.character(varlist) | !is.null(labels)) {
    if(is.null(labels) | ((length(labels)!=1) & (length(varlist)!=length(labels)))) {
      display.var(varlist,"unknown variable"); cat("\n")
    } else { 
      for(cc in 1:length(labels)){
        if(is.list(counts)) { cnts <- lapply(counts,"[",cc) } else { cnts <- NULL }
        if(is.list(varlist)) {
          display.var(varlist[[cc]],labels[cc],cnts=cnts)
        } else {
          display.var(varlist[cc],labels[cc],cnts=cnts)
        }
        cat("\n") 
      }
      return(invisible())
    }
    return(invisible())
  }
  if(prv.call) { gens <- 2 } else { gens <- 1 }
  ENVIR <- parent.frame(n=gens)
  for(cc in 1:length(lab)) {
    label <- lab[cc]
    #print(sys.parent())
    #print(sys.nframe())
    #print(sys.frame(-1))#
    mymode <- "any"
    if(exists(label,mode="function")) { if(exists.not.function(label)) { 
      mymode <- exists.not.function(label,T) } } # if object is also a function, what type is the other type?
    #if(mymode=="") { mymode <- "any" }
    val <- NULL
    try(val <- get(label,envir=ENVIR, mode=mymode),silent=T)
    sf <- sys.frames(); cc <- 1
    while(is.null(val) & cc<=length(sf)) { (try(val <- get(label,envir=sf[[cc]],mode=mymode),silent=T)); cc <- cc + 1 }
    if(!is.null(val)) {
      if(is.list(counts)) { cnts <- lapply(counts,"[",cc) } else { cnts <- NULL }
      display.var(val,label,cnts=cnts)
      cat("\n") 
    } else {
      cat("preview() couldn't find variable '",label,"'\n",sep="")
    }
  }
  return(invisible())
}




#' Simulate a dataset with correlated measures
#'
#' Simulate a dataset with correlated measures (normal simulation with e.g, rnorm() usually
#'  only gives small randomly distributed correlations between variables). This is a quick
#'  and unsophisticated method, but should be able to provide a dataset with slightly more
#'  realistic structure than simple rnorm() type functions. Varying the last three parameters
#'  gives some control on the way the data is generated. It starts with a seed random variable,
#'  then creates 'k' random variables with an expected correlation of r=genr() with that seed 
#'  variable. Then after this, one of the variables in the set (including the seed) is randomly
#'  selected to run through the same process of generating 'k' new variables; this is repeated
#'  until columns are full up. 'mix.order' then randomizes the column order destroying the
#'  relationship between column number and correlation structure, although in some cases,
#'  such relationships might be desired as representative of some real life datasets. 
#' @param nrow integer, number of rows to simulate
#' @param ncol integer, number of columns to simulate
#' @param genx the generating function for data, e.g rnorm(), runif(), etc
#' @param genr the generating function for desired correlation, e.g, runif()
#' @param k number of steps generating from the same seed before choosing a new seed
#' @param mix.order whether to randomize the column order after simulating
#' @export
#' @seealso \code{\link{cor.with}}
#' @author Nicholas Cooper 
#' @examples
#' corDat <- sim.cor(200,5)
#' prv(corDat) # preview of simulated normal data with r uniformly varying
#' cor(corDat) # correlation matrix
#' corDat <- sim.cor(500,4,genx=runif,genr=function(x) { 0.5 },mix.order=FALSE)
#' prv(corDat) # preview of simulated uniform data with r fixed at 0.5
#' cor(corDat) # correlation matrix
sim.cor <- function(nrow=100,ncol=100,genx=rnorm,genr=runif,k=3,mix.order=TRUE) {
  #ncol <- 100
  #nrow <- 100
  new.mat <- matrix(numeric(ncol*nrow),nrow=nrow)
  X <- genx(nrow)
  new.mat[,1] <- X
  cnt <- 0
  for (cc in 2:ncol) {
    dd <- cor.with(X,r=genr(1))
    new.mat[,cc] <- dd
    cnt <- cnt+1
    if(cnt>=k) { X <- new.mat[,sample(cc,1)]; cnt <- 0 }
  }
  if(mix.order) {
    new.mat <- new.mat[,sample(ncol(new.mat))]
  }
  return(new.mat)
}


#' Simulate a correlated variable
#'
#' Simulate a variable correlated at level 'r' with cector x (of the same length). Can
#' either 'preserve' the mean and standard-deviation, leave standardizeed, 
#' or select new mean 'mn' and standard deviation 'st'.
#' @param x existing variable, to which you want to simulate a new correlated variable
#' @param r the 'expected' correlation you want to target (randomness 
#'  will mean that the actual correlation will vary around this value)
#' @param preserve logical, whether to preserve the same mean and standard deviation(SD)
#'  as x, for the new variable
#' @param mn optional, set the mean for the new simulated variable [must also set st if using this]
#' @param st optional, set the SD for the new simulated variable [must also set mn if using this]
#' @return return the new variable with an expected correlation of 'r' with x
#' @references http://www.uvm.edu/~dhowell/StatPages/More_Stuff/CorrGen.html
#' @export
#' @seealso \code{\link{sim.cor}}
#' @author Nicholas Cooper 
#' @examples
#' X <- rnorm(10,100,14)
#' cor.with(X,r=.5) # create a variable correlated .5 with X
#' cor(X,cor.with(X)) # check the actual correlation
#' # some variability in the actual correlation, so run 1000 times:
#' print(mean(replicate(1000,{cor(X,cor.with(X))})))
#' cor.with(X,preserve=TRUE) # preserve original mean and standard deviation
#' X[c(4,10)] <- NA # works fine with NAs, but new var will have same missing
#' cor.with(X,mn=50,st=2) # specify new mean and standard deviation
cor.with <- function(x,r=.5,preserve=FALSE,mn=NA,st=NA) {
  # inspired by David C. Howell
  # http://www.uvm.edu/~dhowell/StatPages/More_Stuff/CorrGen.html
  X <- standardize(x)
  L <- length(X)
  y <- rnorm(L)
  a <- r/(sqrt(1-(r^2)))
  Z = a*X + y
  z <- standardize(Z)
  if(preserve) {
    mn <- mean(x,na.rm=T)
    st <- sd(x,na.rm=T)
  }
  if(preserve | (!is.na(mn) & !is.na(st))) {
    z <- (z*st)+mn
  }
  return(z)
}




#' Summarise the dimensions and type of available R example datasets
#' 
#' This function will parse the current workspace to see what R datasets
#' are available. Using the toHTML function from the 'tools' package to interpret
#' the data() call, each dataset is examined in turn for type and dimensionality.
#' Can also use a filter for dataset types, to only show, for instance, matrix 
#' datasets. Also you can specify whether to only look for base datasets, or to
#' search for datasets in all available packages. Result is a printout to the
#' console of the available datasets and their characteristics.
#'
#' @param filter logical, whether to filter datasets by 'types'
#' @param types if filter=TRUE, which data types to include in the result
#' @param all logical, if all=TRUE, look for datasets in all available packages, else just base
#' @param ... if all is false, further arguments to the data() function to search datasets
#' @export
#' @author Nicholas Cooper 
#' @examples
#' summarise.r.datasets()
#' summarise.r.datasets(filter=TRUE,"matrix")
## create a summary of R datasets you could use
summarise.r.datasets <- function(filter=FALSE,types=c("data.frame","matrix"),all=FALSE,...) { 
  # eg., package = .packages(all.available = TRUE)
  if(all) {
    ll <- unlist(strsplit((toHTML(data(package = .packages(all.available = TRUE), envir = environment()))),"\n"))
  } else {
    ll <- unlist(strsplit((toHTML(data(..., envir = environment()))),"\n"))
  }
  ll <- ll[-grep("<",ll,fixed=T)]
  datasets <- ll[-grep(" ",ll,fixed=T)]
  
  for (cc in 1:length(datasets)) { 
    if(exists(datasets[cc])) {
      dd <- NULL
      try(dd <- get(datasets[cc]))
      if(is.null(dd)) { ddd <- ""; iz <- "" } else { ddd <- Dim(dd); iz <- is(dd)[1] }
      if(filter) { if(any(types %in% is(dd))) { disp <- T } else { disp <- F } } else { disp <- T }
      if(disp) {
        cat(paste(datasets[cc])," [",paste(ddd,collapse=","),"] (",iz,")\n",sep="")
      }
    }
  }
}



#' Does object exist ignoring functions
#'   
#' The exists() function can tell you whether an object exists
#' at all, or whether an object exists with a certain type, but
#' it can be useful to know whether an object exists as genuine 
#' data (and not a function) which can be important when a variable
#' or object is accidently or intentionally given the same name as
#' a function. This function usually returns a logical value as to
#' the existence of the object (ignoring functions) but can also
#' be set to return the non-function type if the object exists.
#' @param x the object name to search for
#' @param ret.type logical, if TRUE then will return the objects' type (if it exists) rather
#' than TRUE or FALSE. If the object doesn't exist the empty string will be returned as the type.
#' @return logical, whether non-function object exists, or else the type if ret.type=TRUE
#' @export
#' @author Nicholas Cooper 
#' @examples
#' x <- "test"
#' # the standard exists function, for all modes, correct mode, and other modes:
#' exists("x")
#' exists("x",mode="character")
#' exists("x",mode="numeric")
#' # standard case for a non-function variable
#' exists.not.function("x",TRUE)
#' # compare results for a non-existent variable
#' exists("aVarNotSeen")
#' exists.not.function("aVarNotSeen")
#' # compare results for variable that is a function
#' exists("mean")
#' exists.not.function("mean")
#' # define a variable with same name as a function
#' mean <- 1.4
#' # exists.not.function returns the type of the variable ignoring the function of the same name
#' exists.not.function("mean",TRUE)
#' exists("mean",mode="function")
#' exists("mean",mode="numeric")
exists.not.function <- function(x,ret.type=FALSE) {
  if(!is.character(x)) {
    stop("x should be the name of an object [as character type]")
  }
  other.modes <- c("logical", "integer", "list", "double", "character", "raw", "complex", "NULL")
  ex <- F; type <- ""
  for(cc in 1:length(other.modes)) {
    if(exists(x,mode=other.modes[cc])) { ex <- T ; type <- other.modes[cc] }
  }
  if(ret.type) {
    return(type)
  } else {
    return(ex)
  }
}


#' A more general dimension function
#'
#' A more general 'dim' function. For arrays simply calls the dim() function, but for other data types, tries to
#' provide an equivalent, for instance will call length(x) for vectors, and will
#' recursively report dims for lists, and will attempt something sensible for other datatypes.
#' 
#' @param x the object to find the dimension for
#' @param cat.lists logical, for lists, TRUE will concatenate the dimesions to a single string,
#'  or FALSE will return the sizes as a list of the same structure as the original.
#' @seealso \code{\link{prv}}, \code{\link{preview}}
#' @return dimension(s) of the object
#' @export
#' @examples
#' # create variables of different types to show output styles #
#' Dim(193)
#' Dim(1:10)
#' testvar <- matrix(rnorm(100),nrow=25)
#' Dim(matrix(rnorm(100),nrow=25))
#' Dim(list(first="test",second=testvar,third=100:110))
#' Dim(list(first="test",second=testvar,third=100:110),FALSE)
Dim <- function(x,cat.lists=TRUE) {
  max.dims <- 100
  rez <- NULL
  try(rez <- dim(x))
  if(!is.null(rez)) { return(dim(x)) }
  if(is(x)[1]=="list") { 
    out <- lapply(x,Dim) 
    if(cat.lists) {
      if(length(out)>max.dims) { suf <- paste("... + ",length(out)-max.dims," more",sep="") } else { suf <- "" }
      out <- paste(out[1:min(max.dims,length(out))],collapse="; ")
      out <- paste(out,suf)
    }
  } else { out <- length(x) }
  return(out)  
}


#' Force argument to be a numeric type with length one
#'
#' Sometimes arguments must be numeric, scalar and within a certain range.
#' Rather than using many if statements, this will do everything possible to 
#' coerce input to a scalar, failing that will replace with a default value.
#' Can also provide a maximum and minimum range that the result must lie within.
#' 
#' @param x the object to ensure is a scalar
#' @param default the value to revert to if the format of x is illegal
#' @param min a lower bound for the output, anything below this is set to min
#' @param max an upper bound for the output, anything above this is set to max
#' @seealso \code{\link{force.percentage}}
#' @return the object x if already legal, first element if a vector, the min or
#'  max value if x is outside the specified bounds, or the value of default otherwise
#' @export
#' @examples
#' force.scalar(1.5)
#' force.scalar(NULL,default=.5)
#' force.scalar(NA,default=.4,min=5,max=10) # default is outside range!
#' force.scalar(rnorm(1000))
#' force.scalar(101,max=50)
#' force.scalar(list(0.4,1,2,3,4,"test"))
#' force.scalar(data.frame(test=c(1,2,3),name=c("test","me","few")))
#' force.scalar(Inf)
force.scalar <- function(x,default=1, min=-10^10, max=10^10) {
  if(is.list(x)) { x <- unlist(x) }
  if(!is.numeric(default)) { default <- mean(c(min[1],max[1]),na.rm=T) ; warning("bad default, reverting to max,min mean") }
  if(length(Dim(default))!=1) { default <- mean(c(min[1],max[1]),na.rm=T) ; warning("bad default, reverting to max,min mean") }
  if(length(x)<1) { x <- default }
  if(length(Dim(x))!=1) { x <- default }
  x <- suppressWarnings(as.numeric(x[1]))
  if(!is.numeric(x)) { x <- default }
  if(is.na(x)) { x <- default }
  if(x<min) { x <- min }
  if(x>max) { x <- max }
  return(x)
}


#' Force argument to be a percentage with length one
#'
#' Sometimes it is nice to be able to take a percentage as an argument and not
#' have to specify whether it should be entered as a number between 0 and 100, 
#' e.g, 50 = 50%, or as a decimal between 0 and 1, e.g, 0.5 = 50%. Anything greater
#' than 1 and less than 100 will be divided by 100. Anything outside 0,100 will be
#' set to 0,100 respectively.
#' 
#' @param x the object to ensure is a oercentage
#' @param default the value to revert to if the format of x is illegal
#' @seealso \code{\link{force.scalar}}
#' @return the object x if already legal, first element if a vector, the min or
#'  max value if x is outside the specified bounds, or the value of default otherwise
#' @export
#' @examples
#' # create variables of different types to show output styles #
#' force.percentage(45)
#' force.percentage(450)
#' force.percentage(.45)
#' force.percentage(-45)
#' force.percentage("twenty")
#' force.percentage(NA,default=0.25)
force.percentage <- function(x,default=.5) {
  x <- force.scalar(x,default=default, min=0,max=100)
  while(x>=1) { x <- x/100 }
  return(x)
}


#' Create fake text for testing purposes
#' 
#' Returns randomized input as if reading lines from a file, like 'readLines()'
#' Can be used to test i/o functions, robustness.
#' 
#' @param max.lines maxmimum number of fake lines to read
#' @param max.chars maximum number of characters per line
#' @param pc.space percentage of randomly generated characters that should be a delimiter
#' @param delim what should the simulated delimiter be, e.g, a space, comma etc. If you wish not
#'  to include such either set the delimiter as "", or set pc.space=0.
#' @param can.null whether with probability 1/max.lines to return NULL instead of any lines of text,
#'  which simulates an empty file, which for testing purposes you may want to be able to handle
#' @return a vector of character entries up 'max.chars' long, or sometimes only NULL if can.null=TRUE
#' @export
#' @author Nicholas Cooper
#' @examples
#' fakeLines() # should produce between zero and ten lines of random text, 35% of which are spaces
fakeLines <- function(max.lines=10,max.chars=100,pc.space=.35,delim=" ",can.null=TRUE) {
  all.char <- "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!$%&*()_+-=;'#,./<>?:@~{}[]| "
  ncs <- nchar(all.char)
  pct <- force.percentage(pc.space,default=.35)
  if(pct<=0) { wantx <- 0 } else { wantx <- ncs/((1/pct)-1) }
  all.chars <- c(strsplit(all.char,"")[[1]],rep(delim,times=round(wantx)))
  nlines <- sample(max.lines,1)
  txt <- NULL
  if(!can.null | nlines!=max.lines){
    for(cc in 1:nlines) {
      nc <- sample(max.chars,1)
      lineind <- sample(length(all.chars),size=nc,replace=T)
      charz <- paste(all.chars[lineind],collapse="",sep="")
      txt <- c(txt,charz)
    }
  } else { txt <- c() }
  return(txt)
}



#' Monitor CPU, RAM and Processes
#' 
#' This function runs the unix 'top' command and returns the overall CPU and RAM usage,
#' and optionally the table of processes and resource use for each. Works only with
#' unix-based systems such as Mac OS X and Linux, where 'top' is installed. Default
#' is to return CPU and RAM overall stats, to get detailed stats instead, set Table=TRUE.
#'
#' @param CPU logical, whether to return overall CPU usage information
#' @param RAM logical, whether to return overall RAM usage information
#' @param Table logical, whether to return system information for separate processes. This
#'  is returned as table with all of the same columns as a command line 'top' command. If
#'  'Table=TRUE' is set, then the default becomes not to return the overall CPU/RAM usage stats.
#'  The dataframe returned will have been sorted by descending memory usage.
#' @param procs integer, if Table=TRUE, then the maximum number of processes to return (default 20)
#' @param mem.key character, default for Linux is 'mem' and for Mac OS X, 'physmem', but if the 'top'
#'  command on your system displays memory usage using a different label, then enter it here
#'  (case insensitive) to override defaults.
#' @param cpu.key character, default for Linux and Mac OS X is 'cpu', but if the top
#'  command on your system displays CPU usage using a different label, then enter it here.
#' @return a list containing CPU and RAM usage, or with alternate parameters can return stats for each process
#' @export
#' @author Nicholas Cooper
#' @examples
#' # not run #  top()
#' # not run #  top(Table=TRUE,proc=5)
top <- function(CPU=!Table,RAM=!Table,Table=FALSE,procs=20,mem.key=NULL,cpu.key=NULL) {
  if(!RAM & !CPU & !Table) { warning("Deselected all options, null will be returned"); return(NULL) }
  if(!check.linux.install("top")) {
    warning("'top' command only works on Mac OS X and Linux")
    return(NULL)
  }
  if(toupper(Sys.info()["sysname"])=="DARWIN") { macos <- T } else { macos <- F }
  if(macos) {
    # MAC OS X
    txt <- tryCatch(system("top -l 1",intern=T), error = function(e) e)
    if(length(txt)==0) { warning("command failed"); return(NULL) }
    if(!is.character(txt)) { warning("command failed"); return(NULL) }
    dtt <- divide.top.txt(txt)
    if(!is.list(dtt) | any(is.na(dtt)) ) { warning("unexpected output from top command"); return(NULL) }
    parz <- dtt$table; headr <- dtt$header
    if(!is.character(mem.key)) { mem.key <- "physmem" }
    if(RAM) { ram.gb.list <- suck.mem(headr,key=mem.key) }
  }
  if(!macos) {
    ## LINUX
    txt <- tryCatch(system("top -n 1 -b",intern=T), error = function(e) e)
    if(length(txt)==0) { warning("command failed"); return(NULL) }
    if(!is.character(txt)) { warning("command failed"); return(NULL) }
    dtt <- divide.top.txt(txt)
    if(!is.list(dtt) | any(is.na(dtt)) ) { warning("unexpected output from top command"); return(NULL) }
    parz <- dtt$table; headr <- dtt$header
    if(!is.character(mem.key)) { mem.key <- "mem" }
    if(RAM) { ram.gb.list <- suck.mem(headr,key=mem.key) }
  }
  if(!is.character(cpu.key)) { cpu.key <- "cpu" }
  if(CPU) { cpu.pc.list <- suck.cpu(headr,key=cpu.key) }
  if(Table) {
    tab <- make.top.tab(parz); if(all(is.null(tab))) { return(NULL) }
    mem.col <- grep("mem",colnames(tab),ignore.case=T)[1]
    if(is.na(mem.col)) { mem.col <- grep("RSIZE",colnames(tab),ignore.case=T)[1] }
    cpu.col <- grep("cpu",colnames(tab),ignore.case=T)[1]
    if(is.na(mem.col) | (is.na(mem.col))) {
      warning("did not find 'mem', 'RSIZE' or 'CPU' entries in 'top' output")
    } else {
      tab <- tab[rev(order(tab[,mem.col])),]
      tab <- tab[rev(order(tab[,cpu.col])),]
      tab <- tab[rev(order(tab[,mem.col])),]
    }
    if(is.na(as.numeric(procs))) { procs <- nrow(tab) } else { procs <- round(procs) }
    procs <- min(c(procs,nrow(tab)),na.rm=T)
    if(is.null(dim(tab))) { Table <- FALSE }
  }
  outlist <- NULL; outnms <- NULL
  if(CPU) { outlist <- c(outlist,list(cpu.pc.list)); outnms <- c(outnms,"CPU") }
  if(RAM) { outlist <- c(outlist,list(ram.gb.list)); outnms <- c(outnms,"RAM") }
  if(Table) { outlist <- c(outlist,list(tab[1:procs,])); outnms <- c(outnms,"Table") }
  names(outlist) <- outnms
  return(outlist)
}


# internal function to support top() function
make.top.tab <- function(parz) {
  if(!is.list(parz)) { warning("unexpected input"); return(NULL) }
  cnts <- sapply(parz,length)
  exp.lines <- Mode(cnts)
  shortz <- which(cnts<exp.lines)
  longz <- which(cnts>exp.lines)
  if(length(longz)>0) {  parz[longz] <- lapply(parz[longz],function(X) { X[1:exp.lines] }) }
  if(length(shortz)>0) { parz <- parz[-shortz] }
  LL <- length(parz[[1]]); if(LL < 1) { warning("unexpected input"); return(NULL) }
  df <- as.data.frame(matrix(ncol=LL,nrow=length(parz)))
  if(nrow(df)<1) { warning("unexpected input"); return(NULL)  }
  for(cc in 1:length(parz[[1]])) { df[,cc] <- sapply(parz,"[",cc) }
  if(nrow(df)<2) { warning("expected header row and at least 1 data row"); return(NULL)  }
  tab <- df[-1, ,drop=FALSE]; colnames(tab) <- df[1,]; rownames(tab) <- NULL
  return(tab)
}

# internal function to support top() function
divide.top.txt <- function(txt) {
  parz <- strsplit(txt," +|\t")
  parz <- lapply(parz,function(X) { X <- X[!is.na(X)] ; X[X!=""] } ) 
  headline <- which(sapply(parz,function(X) { all(c("PID","USER") %in% toupper(X)) }))
  if(length(headline)<1) { warning("expected PID and USER column - at least 1 not found"); return(NA) }
  parz <- parz[headline[1]:length(parz)]
  headr <- txt[1:(headline[1]-1)]
  return(list(header=headr,table=parz))
}

# internal function to support top() function
suck.num.from.txt <- function(txt) {
  if(is.na(txt)) { return(NA) }
  if(length(txt)<1) { return(NA) }
  splt <- strsplit(txt,"")
  nmall <- numeric()
  anm <- function(X) { suppressWarnings(as.numeric(X)) }
  for(cc in 1:length(splt)) {
    nm <- sapply(splt[[cc]],function(X) {
      if(!is.na(anm(X))) { anm(X) } else { if(X==".") { X } else { NA } } } )
    nmall[cc] <- anm(paste(narm(nm),collapse="",sep=""))
  }
  return(nmall)
}

# internal function to support top() function
suck.cpu <- function(headr,key="cpu") {
  cpz <- grep(key,headr,ignore.case=T)
  if(length(cpz)>0) {
    cpuline <- headr[cpz[1]]
    ms <- strsplit(cpuline,",")[[1]]
    ms <- gsub("cpu","",ms,ignore.case=T)
    user <- ms[grep("us",ms,ignore.case=T)]
    sys <- ms[grep("sy",ms,ignore.case=T)]
    idle <- ms[grep("id",ms,ignore.case=T)]
    if(length(user)>0) {
      user1 <- rmv.spc(gsub("us","",gsub("user","",user,ignore.case=T)))
      user.gb <- suck.num.from.txt(user1)
    } else { user.gb <- NA }
    if(length(sys)>0) {
      sys1 <- rmv.spc(gsub("sy","",gsub("sys","",sys,ignore.case=T)))
      sys.gb <- suck.num.from.txt(sys1)
    } else { sys.gb <- NA }
    if(length(idle)>0) {
      idle1 <- rmv.spc(gsub("id","",gsub("idle","",idle,ignore.case=T)))
      idle.gb <- suck.num.from.txt(idle1)
    } else { idle.gb <- NA }
    if(is.na(idle.gb) & !is.na(sys.gb) & !is.na(user.gb)) { idle.gb <- 100-user.gb-sys.gb }
    if(is.na(sys.gb) & !is.na(idle.gb) & !is.na(user.gb)) { sys.gb <- 100-user.gb-idle.gb }
    if(is.na(user.gb) & !is.na(sys.gb) & !is.na(idle.gb)) { user.gb <- 100-idle.gb-sys.gb }
  } else { 
    cat("no CPU usage information found\n")
    return(NULL)
  }
  return(list(total=user.gb,idle=idle.gb,sys=sys.gb,unit="%"))
}

# internal function to support top() function
suck.mem <- function(headr,key="Mem") {
  memz <- grep(key,headr,ignore.case=T)
  if(length(memz)>0) {
    memline <- headr[memz[1]]
    ms <- strsplit(memline,",")[[1]]
    ms <- gsub("mem","",ms,ignore.case=T)
    tot <- ms[grep("total",ms,ignore.case=T)]
    free <- ms[grep("free",ms,ignore.case=T)]
    used <- ms[grep("used",ms,ignore.case=T)]
    if(length(tot)>0) {
      tot1 <- rmv.spc(gsub("total","",tot,ignore.case=T))
      tot.gb <- suck.bytes(tot1)
    } else { tot.gb <- NA }
    if(length(free)>0) {
      free1 <- rmv.spc(gsub("free","",free,ignore.case=T))
      free.gb <- suck.bytes(free1)
    } else { free.gb <- NA }
    if(length(used)>0) {
      used1 <- rmv.spc(gsub("used","",used,ignore.case=T))
      used.gb <- suck.bytes(used1)
    } else { used.gb <- NA }
    if(is.na(used.gb) & !is.na(free.gb) & !is.na(tot.gb)) { used.gb <- tot.gb-free.gb }
    if(is.na(free.gb) & !is.na(used.gb) & !is.na(tot.gb)) { free.gb <- tot.gb-used.gb }
    if(is.na(tot.gb) & !is.na(free.gb) & !is.na(used.gb)) { tot.gb <- used.gb+free.gb }
  } else { 
    cat("no RAM usage information found\n")
    return(NULL)
  }
  return(list(total=tot.gb,used=used.gb,free=free.gb,unit="Gb"))
}

# internal function to support top() function  
suck.bytes <- function(tot1,GB=TRUE) {
  mult <- 0
  if(length(grep("k",tot1,ignore.case=T))>0) { mult <- 1000 }
  if(length(grep("m",tot1,ignore.case=T))>0) { mult <- 10^6 }
  if(length(grep("g",tot1,ignore.case=T))>0) { mult <- 10^9 }
  if(mult==0 & length(grep("b",tot1,ignore.case=T))>0) { mult <- 1 }
  if(mult==0) { warning("expected symbol indicating units, defaulting to bytes"); mult <- 1 }
  lst <- c("kb","gb","mb","b","g","m","k")
  tot1 <- suck.num.from.txt(tot1)
  if(is.na(tot1)) { tot1 <- 0 ; warning("no numbers found in text, setting to zero") }
  tot2 <- (as.numeric(tot1)*mult)/10^9 ; 
  if(!GB) { tot2 <- tot2/10^3 }
  return(tot2)
}




#' Check whether a set of packages has been loaded
#' 
#' Returns TRUE if the whole set of packages entered has been loaded, or FALSE
#' otherwise. This can be useful when developing a package where there is optional
#' functionality depending if another package is in use (but the other package is
#' not part of 'depends' because it is not essential). Because 'require' cannot
#' be used within functions submitted as part of a CRAN package.
#' @param pcks character, a package name, or vector of names, if left blank will return all loaded
#' @param ... further package names as character (same as entering via pcks, 
#'  but avoids need for c() in pcks)
#' @param cran.check logical, in the case at least one package is not found, whether
#'  to search CRAN and see whether the package(s) even exist on CRAN.
#' @param repos repository to use if package is not loaded and cran.check=TRUE,
#'  if NULL, will attempt to use the repository in getOptions("repos") or will
#'  default to the imperial.ac.uk mirror. Otherwise the default is to use
#'  all available repositories from getRepositories()
#' @return logical TRUE or FALSE whether the whole list of packages are available
#' @export
#' @author Nicholas Cooper 
#' @examples
#' require(BiocInstaller)
#' packages.loaded("NCmisc","reader")
#' packages.loaded(c("bigpca","nonsenseFailTxt")) # both not found, as second not real
#' packages.loaded(c("bigpca","nonsenseFailTxt"),cran.check=FALSE) # hide warning
#' packages.loaded() # no argument means all loaded packages are listed
#' packages.loaded("snpStats",repos=getRepositories(1)) # doesn't find the bioconductor package on CRAN
#' packages.loaded("snpStats",repos=getRepositories()) # now it can find it by using all repositories
packages.loaded <- function(pcks="",...,cran.check=TRUE,repos=getRepositories()) {
  more <- c(...); if(length(more)>0) { pcks <- c(pcks,paste(more)) }
  if(!is.character(pcks)) { stop("must enter package names as character strings") }
  pt <- "package:"; pkgset <- gsub(pt,"",search()[grep(pt,search(),fixed=TRUE)])
  if(all(pcks=="")) { return(pkgset) }
  answer <- (all(pcks %in% pkgset))
  if(is.null(repos)) { try(repos <- getOption("repos") ) }
  if(is.null(repos)) { repos <- "http://cran.ma.imperial.ac.uk/" }
  #print(repos)
  if(!answer & cran.check) {
    check.exist <- search.cran(pcks,repos=repos)
    for(cc in 1:length(check.exist)) {
      if(!pcks[cc] %in% check.exist[[cc]]) { cat("NB: package",pcks[cc],"is not on CRAN\n") }
    }
  }
  return(answer)
}


#' Split a text file into multiple parts
#' 
#' Wrapper for the bash command 'split' that can separate a text file into multiple 
#' roughly equal sized parts. This function removes the need to remember syntax and
#' suffixes of the bash command
#' @param fn character, file name of the text file to split, if the file is an imcompatible format
#'  the linux command should return an error message to the console
#' @param size integer, the maximum number of lines for the split parts of the file produced
#' @param same.dir logical, whether the resulting files should be moved to the same
#'  directory as the original file, or simply left in the working directory [getwd()]
#' @param verbose logical, whether to report the resulting file names to the console
#' @param suf character, suffix for the split files, default is 'part', the original file
#'  extension will be appended after this suffix
#' @param win logical, set to FALSE if running a standard windows setup (cmd.ext), and the file
#' split will run natively in R. Set to TRUE if you have a unix-alike command system, such as
#' CygWin, sh.exe, csh.exe, tsh.exe, running, and this will then check to see whether the POSIX
#' 'split' command is present (this provides a speed advantage). If in doubt, windows users
#' can always set win=TRUE; the only case where this will cause an issue is if there is a
#' different command installed with the same name (i.e, 'split').
#' @export
#' @return returns the list of file names produced (including path)
#' @author Nicholas Cooper 
#' @examples
#' orig.dir <- getwd(); setwd(tempdir()); # move to temporary dir
#' file.name <- "myfile.txt"
#' writeLines(fakeLines(max.lines=1000),con=file.name)
#' new.files <- file.split(file.name,size=50)
#' unlink(new.files); unlink(file.name)
#' setwd(orig.dir) # reset working dir to original
file.split <- function(fn,size=50000,same.dir=FALSE,verbose=TRUE,suf="part",win=TRUE) {
  if(!file.exists(fn)) { stop("file",fn,"did not exist")}
  if(!is.numeric(size)) { stop("size must be numeric") }
  size <- as.integer(round(size))
  FN <- basename(fn)
  EXT <- get.ext(fn)
  DIR <- dirname(fn)
  if(!same.dir) { DIR <- getwd() }
  file.out <- paste(rmv.ext(FN),ext=suf,sep="_")
  if(!all(check.linux.install(c("split","mv")))) { 
    new.fnz <- suppressWarnings(file.split.windows(fn,size,file.out,DIR,EXT,verbose))
    if(length(new.fnz)<1) { stop("no split part files produced, operation failed") }
  } else {
    cmd <- paste0("split -l ",size," ",fn," ",file.out)
    st <- proc.time()[3]
    jj <- suppressWarnings(suppressMessages(system(cmd,intern = TRUE, ignore.stderr = TRUE)))
    tot <- proc.time()[3]-st
    if(tot>3) {
      cat("command '",cmd,"' was run using system()\n",sep="")
    }
    new.filez <- list.files(pattern=file.out)
    if(length(new.filez)<1) { stop("no split part files produced, operation failed") }
    tt <-  new.filez %in% cat.path("",new.filez,ext=EXT)
    if(any(tt)) {
      # files from a previous run may already be in the directory already with an extension
      new.filez <- new.filez[-which(tt)]
    }
    new.fnz <- cat.path(DIR,new.filez,ext=EXT)
    for (dd in 1:length(new.filez)) {
      system(paste0("mv ",new.filez[dd]," ",new.fnz[dd]))
    }
    if(verbose) {
      cat("split ",fn," into ",length(new.filez)," parts:\n  ",paste(new.fnz,collapse="\n  "),"\n",sep="")
    }
  }
  return(new.fnz)
}


# internal alternative to the split command for windows
file.split.windows <- function(fn,size,file.out,DIR,EXT,verbose=TRUE) {
  read.file <- file(fn)
  open(con=read.file,open="r")
  write.file <- file(cat.path(DIR,file.out,suf=1,ext=EXT))
  open(con=write.file,open="w")
  eof <- FALSE; cc <- -1; dd <- 0
  filenum <- 1
  while(!eof) {
    cc <- cc + 1; dd <- dd + 1
    next.line <- readLines(con=read.file,n=1)
    eof <- length(next.line)==0
    writeLines(next.line,con=write.file)
    if(dd>(size-1) & !eof) {
      filenum <- filenum + 1
      close(con=write.file)
      write.file <- file(cat.path(DIR,file.out,suf=filenum,ext=EXT))
      open(con=write.file,open="w")
      dd <- 0
    }
  }
  close(con=read.file)
  close(con=write.file)
  if(verbose) { cat("split",cc,"lines of",fn,"into",filenum,"files:\n") }
  return(cat.path(DIR,file.out,suf=1:filenum,ext=EXT))
}


#INTERNAL
## COPY FROM READER, SO NCMISC DOESN'T DEPEND ON READER
rmv.ext <- function(fn=NULL,only.known=TRUE,more.known=NULL,print.known=FALSE) {
  # remove file extension from a filename character string
  known.ext <- c("TXT","RDATA","TAB","DAT","CSV","VCF","GCM","BIM","MAP","FAM",
                 "PFB","SH","R","CPP","H","DOC","DOCX","XLS","XLSX","PDF","JPG",
                 "BMP","PNG","TAR","GZ","CNV","PL","PY","ZIP","ORG","RDA","DSC","BCK",
                 "ABW","HTM","HTML",toupper(more.known))
  if(is.null(fn)) { 
    if(print.known) {
      return(known.ext)
    } else {
      warning("couldn't remove extension, not a character()"); return(fn) 
    }
  } else {
    if (all(is.na(fn))) { warning("couldn't remove extension, all values were NA"); return(fn) }
  }
  if(print.known) { cat("known file extensions:\n"); print(known.ext) }
  if(!is.character(fn)) { warning("couldn't remove extension, not a character()"); return(fn) }
  rmv.one <- function(X,known.ext) {
    file.segs <- strsplit(paste(X),".",fixed=TRUE)[[1]]
    lss <- length(file.segs)
    if (lss>1) { 
      if(only.known){
        if(toupper(file.segs[lss]) %in% known.ext) {
          out <- paste(file.segs[-lss],collapse=".") 
        } else { 
          out <- X
        }
      } else {
        out <- paste(file.segs[-lss],collapse=".") 
      }
    } else {
      out <- X 
    }
  }
  return(sapply(fn,rmv.one,known.ext=known.ext))
}

#INTERNAL
## COPY FROM READER, SO NCMISC DOESN'T DEPEND ON READER
cat.path <- function(dir="",fn,pref="",suf="",ext="",must.exist=FALSE) 
{
  dir.ch <- .Platform$file.sep
  if(is.list(fn) & is.ch(fn)) { fn <- unlist(fn) } #; 
  if(length(dir)>1) { dir <- dir[1]; cat("only first dir was used\n") }
  if(length(ext)>1) { ext <- ext[1]; cat("only first extension was used\n") }
  if(length(grep(dir.ch,fn))>0) {
    dir <- dirname(fn)  #split into dir and fn if fn has /'s
    fn <- basename(fn)
  }
  dir <- dir.force.slash(dir)
  if(ext!="") {
    #make sure ext includes the dot
    if(substr(ext,1,1)!=".")   { ext <- paste(".",ext,sep="") }
    #if ext is already built into suffix or filename, remove it from there
    fn <- rmv.ext(paste(fn))
    suf <- rmv.ext(paste(suf))
  }
  location <- paste(dir,pref,fn,suf,ext,sep="")
  if(any(!file.exists(location)) & must.exist) {
    warn <- paste("required file",location,"not found!")
    stop(warn)
  }
  return(location)
}

#INTERNAL
## COPY FROM READER, SO NCMISC DOESN'T DEPEND ON READER
#' Internal function used by cat.path
dir.force.slash <- function(dir) {
  # make sure 'dir' directory specification ends in a / character
  if(!is.null(dim(dir))) { stop("dir should be a vector") }
  dir <- paste(dir)
  dir.ch <- .Platform$file.sep
  the.test <- (dir!="" & substr(dir,nchar(dir),nchar(dir))!=dir.ch)
  dir[the.test] <- paste(dir[the.test],dir.ch,sep="")
  return(dir)
}


#INTERNAL
## COPY FROM READER, SO NCMISC DOESN'T DEPEND ON READER
#' Internal function to assess whether data is a character or list of characters
is.ch <- function(x) { 
  # is function for character() or list of characters
  if(is.null(x)) { return(FALSE) }
  pt1 <- is.character(x)
  if(!pt1 & is.list(x)) { pt2 <- all(sapply(x,is.ch)) } else { pt2 <- pt1 }
  return(as.logical(pt1 | pt2))
}



