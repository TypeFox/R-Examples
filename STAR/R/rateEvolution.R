rateEvolution <- function(x,
                          bw,
                          kernel = c("gaussian", 
                            "epanechnikov", "rectangular", "triangular", "biweight", 
                            "cosine", "optcosine"), 
                          n = 512,
                          from, to,
                          na.rm = FALSE, 
                          ...) {

  data.name <- deparse(substitute(x))
  ## check if x is a spikeTrain object
  if (!is.spikeTrain(x)) x <- as.spikeTrain(x)

  if (length(x) > 1) {
    if (missing(bw)) bw <- 10*median(diff(x))

    class(x) <- NULL
    result <- density(x,bw,kernel=kernel,n=n,from=from,to=to,na.rm=na.rm,...)
    result$y <- result$y*result$n
    result$data.name <- data.name
  } else {
    if (!missing(from) && !missing(to)) {
      xx <- seq(from,to,length.out=n)
      if (missing(bw)) bw <- 3*diff(xx)[1]
      class(x) <- NULL
      result <- density(x,bw,kernel=kernel,n=n,from=from,to=to,na.rm=na.rm,...)
      result$y <- result$y*result$n
      result$data.name <- data.name
    } else {
      result <- list(x=NA,y=NA,bw=NA,n=NA,call=NA,data.name,has.na=NA)
    }
  }
  
  class(result) <- c("rateEvolution","density")
  result
}

plot.rateEvolution <- function(x,
                               main = NULL,
                               xlab = NULL,
                               ylab = "Rate (Hz)",
                               type = "l", 
                               zero.line = TRUE,
                               ...) {
  
  if (is.null(main)) main <- paste(x$data.name,"Rate Evolution")
  if (is.null(xlab)) xlab <- "Time (s)"
  NextMethod("plot",,main=main,xlab=xlab,ylab=ylab,type=type,zero.line=zero.line,...)
}


mkREdf <- function(x,
                   longitudinal,
                   across,
                   bw,
                   kernel=c("gaussian", 
                     "epanechnikov", "rectangular", "triangular", "biweight", 
                     "cosine", "optcosine"), 
                   n=512,
                   from, to,
                   na.rm=FALSE,
                   minusMean=FALSE
                   ) {

  ## check argument x
  if (!inherits(x,"list")) x <- list(x=as.spikeTrain(x))
  if (!is.spikeTrain(x[[1]]) && !is.repeatedTrain(x[[1]]))
    stop("Wrong x.")

  xST <- is.spikeTrain(x[[1]])
  if (missing(longitudinal)) {
    if (xST) longitudinal <- names(x)
    else longitudinal <- names(x[[1]])
  }
  if (missing(across)) {
    if (xST) across <- deparse(substitute(x))
    else across <- names(x)
  } 
  allDF <- expand.grid(longitudinal,across)
  nbT <- dim(allDF)[1]

  theList <- list(kernel=kernel,n=n,na.rm=na.rm)
  if (!missing(bw)) {
    bwGiven <- TRUE
    bw <- rep(bw,length.out=nbT)
  } else {
    bwGiven <- FALSE
  }
  if (!missing(from)) {
    theList$from <- from
  } else {
    if (xST) theList$from <- floor(min(sapply(x,min)))
    else theList$from <- floor(min(sapply(x,function(l) min(sapply(l,min)))))
  }
  if (!missing(to)) {
    theList$to <- to
  } else {
    if (xST) theList$to <- ceiling(max(sapply(x,max)))
    else theList$to <- ceiling(max(sapply(x,function(l) max(sapply(l,max)))))
  }
  
  myRE <- function(idx) {
    if (xST) list4call <- c(list(x=x[[idx]]),theList)
    else list4call <- c(list(x=x[[allDF[idx,2]]][[allDF[idx,1]]]),
                        theList)
    if (bwGiven) list4call$bw <- bw[idx]
    re <- do.call(rateEvolution,list4call)[c("x","y")]
    re$longitudinal <- rep(allDF[idx,1],length(re$x))
    re$across <- rep(allDF[idx,2],length(re$x))
    re
  }
                        
  result <- lapply(1:nbT, myRE)
  result <- data.frame(time=unlist(lapply(result,function(l) l$x)),
                       rate=unlist(lapply(result,function(l) l$y)),
                       longitudinal=unlist(lapply(result,function(l) l$longitudinal)),
                       across=unlist(lapply(result,function(l) l$across))
                       )

  if (minusMean) {
    meanR <- sapply(across,
                    function(n)
                    apply(sapply(longitudinal,
                                 function(l)
                                 result$rate[result$across==n & result$longitudinal==l]
                                 ),1,mean)
                    )
    colnames(meanR) <- across
    result$rate <- as.vector(sapply(across,
                                    function(n)
                                    sapply(longitudinal,
                                           function(l)
                                           result$rate[result$across==n & result$longitudinal==l]-
                                           meanR[,n]
                                           )
                                    )
                             )
  }
  result
  
}
