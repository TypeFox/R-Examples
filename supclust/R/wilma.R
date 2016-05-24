## Wilma - an algorithm for supervised grouping of predictor variables
wilma <- function(x, y, noc, genes = NULL, flip = TRUE,
                  once.per.clust = FALSE, trace = 0)
{
    ## Checking the input
    y <- as.integer(y)
    if(any(y < 0 | y > 1))
        stop("Labels y have to be 0 or 1!")
    n <- length(y)
    if(!is.matrix(x) || !is.numeric(x) || nrow(x) != n)
        stop("'x' must be a numeric matrix with 'n' (= length(y)) observations")
    if((noc <- as.integer(noc)) < 1)
        stop("'noc' must be a positive integer")
    if((trace <- as.integer(trace)) < 0)
        stop("'trace' must be integer >= 0 (or logical)")
    ## C output (Cverb > 0) only for trace >= 2 :
    Cverb <- as.integer(if(trace) Cverb <- trace - 1 else 0)

    ## Customizing the problem and sign-flipping
    iy    <- sort.list(y)# i.e., first the 0's, then the 1's
    io    <- order(iy)
    y     <- y[iy]
    x     <- x[iy,]
    signs <- NULL
    if (flip)
      {
        res   <- sign.flip(x,y)
        x     <- res$flipped.matrix
        signs <- res$signs
      }

    ## Definitions
    n1 <- sum(y == 0)
    n2 <- n - n1 ## = sum(y == 1)
    p  <- ncol(x) ## = length(x)/n


    ## Looking for (optional) starting genes
    used <- rep(FALSE, p)
    if(!is.null(genes)) {
        if(!is.list(genes) || length(genes) != noc)
            stop("starting 'genes' must be a list of length 'noc' (=", noc,")")
        sgl <- unlist(genes)
        ##was sgl <- NULL ; for (i in 1:noc) sgl   <- c(sgl, genes[[i]])
        used[sgl] <- TRUE
    }

    mn.x <- gList <- vector("list", noc)
    steps <- integer(noc)

    ## The loop for supervised clustering
    for (i in 1:noc)
    {
        ## Initial value if (optional) starting clusters are provided
        size <- length(gic <- genes[[i]])
        clMean <- rep(0,n)

        if(trace) cat("\n", "\nCluster ", i, "\n----------\n", sep="")

        ## Forward search and cleaning stages
        nFwd <- 1
        repeat { ## Search forward and backward once :

            if (size > 0)
                clMean <- rowMeans(x[, gic, drop=FALSE])

            if(trace && nFwd > 1) {
                cat("used[]:", which(as.logical(used)),"\n")
                if(trace >= 2 && size > 0)
                    cat("entry gic[]:", gic, "\n")
            }

            res <- .C(R_multicluster,
                      as.double(x), as.integer(y),# 2
                      as.integer(n), as.integer(n1), as.integer(n2),# 5
                      as.integer(p),# 6
                      used = as.integer(used),# 7
                      as.double(clMean),
                      glsize = size,# 9
                      gic    = c(as.integer(gic), integer(p)),
                      scores = integer(size+p),
                      margins=  double(size+p),
                      as.logical(once.per.clust),
                      Cverb)[
                      c("used", "glsize", "gic", "scores", "margins")]

            if(trace) {
                if(nFwd > 1 && trace >= 2 && size > 0)
                    cat("exit gic[1:res$glsize]:", res$gic[1:res$glsize], "\n")
                if(res$glsize > size) {
                    cat("\nAccepted",if(size > 0) "(additionally)","\n")
                    for(j in (size+1):res$glsize)
                        .p.1gen(j,  res$ gic [j],
                                s = res$ scores[j],
                                m = res$margins[j], digits = 3)
                    cat(" gic size changed from ", size," to ",res$glsize,"\n")
                }
                else ## res$glsize == size
                    cat(" gic[] *unchanged*\n")
            }
            gic.old <- gic
            size <- res$glsize
            gic  <- res$gic[1:size]
            used <- res$used

            if(trace) cat("\nEliminating ")
            gic.red  <- back.search(gic, x, y, verbose = trace > 0)
            if(length(gic.red) == size) {
                if(trace)
                    cat(" -- no reduction --> end{repeat} after ",
                        nFwd, if(nFwd == 1)"step" else "steps","\n")
                break
            }

            ## else : *have* reduced

            w <- gic[!(gic %in% gic.red)]
            if(trace) cat(" w= (", paste(w,collapse=','),")", sep="")
            used[w] <- FALSE
            gic <- gic.red
            size <- length(gic)

            if(trace)
                cat(" -- end {repeat} nr. ", nFwd,
                    ". gic enlarged and reduced from ", length(gic.old),
                    " to ", size, "\n", sep="")

            if (nFwd > 1 && length(gic.old) == size && all(gic.old == gic))
                break
            nFwd <- nFwd + 1

        }## end{repeat}

        gList[[i]] <- gic
        steps [i] <- nFwd
        mn.x [[i]] <-
            if(size) sapply(1:size,
                            function(j) rowMeans(x[, gic[1:j], drop=FALSE]))
            else integer(0)

	if(trace)
	    p.1clust(i, gic = gic, x.mean = mn.x[[i]], y = y, nFwd = nFwd)

    }## end for i = 1:noc

    ## Restoring the original order of output and response
    y <- y[io]
    for (i in 1:length(mn.x)) mn.x[[i]] <- mn.x[[i]][io,,drop=FALSE]

    ## Output
    r <- list(clist = gList,
              steps = steps, y = y, x.means = mn.x, noc = noc, signs=signs)
    class(r) <- "wilma" # cluster List
    r
}

## Short overview of what has been done by Wilma
print.wilma <- function(x, ...)
  {
    ## The number of clusters which were fitted
    noc <- x$noc

    ## Final score and margin values
    fvals <- fitted(x)
    s     <- apply(fvals, 2, score, x$y)
    m     <- apply(fvals, 2, margin, x$y)

    ## Formatting the numbers
    nGenes <- unlist(lapply(x$clist,length))
    cScore <- format(s)
    cMargi <- format(round(m,2))
    cG     <- format(nGenes)
    cI     <- format(1:noc)

    ## General overview
    cat("\nWilma called to fit ", x$noc, " cluster", if(noc>1) "s",
        "\n\n", sep="")

    ## Information for each of the predictors
    for (i in 1:noc)
      {
	gic <- x$clist[[i]]
        ng <- length(gic)
        cat("Cluster", cI[i], ": Contains ")
        cat(cG[i], " gene", if(ng != 1)"s", ", final score ", cScore[i],
            ", final margin ", cMargi[i], "\n", sep="")
      }
  }


## The next 3 functions provide a detailed overview about Wilma's clustering
summary.wilma <- function(object, ...)
  {
    cat("'Wilma' object: ")
    printClist(object, ...)
  }

## Auxiliary function for summary.wilma()
printClist <- function(x, ...)
  {
    ## Purpose: print (Method) "gene list" objects, called from print.multiclust
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 21 Jun 2003, 21:38

    if(is.null(gic <- x$clist))
        stop("invalid 'x' argument")
    noc <- length(gic)
    cat("number of clusters 'noc' =", noc,"\n")
    for(i in 1:noc)
        p.1clust(i, gic = gic[[i]],
                 x.mean = x$x.means[[i]], y = x$y, nFwd = x$steps[i])
    invisible(x)
  }

## Another auxiliary function for summary.wilma()
p.1clust <- function(i, gic, x.mean, y, nFwd)
  {
    cat("\nFinal Cluster", i,
        "\n----------------\n\n")
    size <- length(gic)
    if(size >= 1)
        for (j in 1:size) {
            s <- score (x.mean[,j], y)
            m <- margin(x.mean[,j], y)
            .p.1gen(j=NULL, id = gic[j], s, m, digits = 3)
        }
    else cat(" __empty__ \n")
    invisible()
  }

.p.1gen <- function(j, id, s, m, digits = 3)
  {
    cat("Gen",
        if(is.null(j)) "" else paste("[", j, "]",if(j <= 9)" ",sep=""),
        ":", formatC(id, wid=6),
        "  Score:", formatC(s, wid=4),
        "  Margin:",
        formatC(format(c(round(m, digits=digits),0.735))[1], wid=7, dig=3),
        "\n", sep="")
  }



## Returns the fitted values
fitted.wilma <- function(object, ...)
  {
    ## Fitted values
    out <- NULL
    for (i in 1:object$noc)
      {
        out <- cbind(out, (object$x.means[[i]])[,ncol(object$x.means[[i]])])
      }

    if (is.null(vNames <- colnames(out))) ## Naming the predictors
        vNames <- paste("Predictor", 1:object$noc)
    dimnames(out) <- list(1:length(object$y), vNames)
    out
}

## 2-dimensional projection of Wilma's output
plot.wilma <- function(x, xlab = NULL, ylab = NULL, col = seq(x$yvals),
                       main = "2-Dimensional Projection of Wilma's Output", ...)
  {
    ## Fitted values
    fvals <- fitted(x)

    ## If only 1 cluster was fitted
    if(x$noc <= 1)
      {
        if(is.null(xlab)) xlab <- "Mean Expression of Wilma's Cluster 1"
        if(is.null(ylab)) ylab <- "Class label"
        plot(fvals[,1], x$y, type = "n", xlab=xlab, ylab=ylab, main=main, ...)
        text(fvals[,1], x$y, x$y, col = col[1 + x$y])
        return()
      }

    ## For 2 clusters and more
    if(is.null(xlab)) xlab <- "Mean Expression of Wilma's Cluster 1"
    if(is.null(ylab)) ylab <- "Mean Expression of Wilma's Cluster 2"
    plot(fvals[,1], fvals[,2], type="n", xlab=xlab, ylab=ylab, main=main, ...)
    text(fvals[,1], fvals[,2], x$y, col=col[1 + x$y])
    invisible()
}


predict.wilma <- function(object, newdata = NULL, type = c("fitted", "class"),
                          classifier = c("nnr", "dlda", "logreg", "aggtrees"),
                          noc = object$noc, ...)
  {
    ## Checking the input
    type       <- match.arg(type)
    classifier <- match.arg(classifier)

    ## Return fitted values for the training data
    if (is.null(newdata))
      {
        if (type == "fitted") {
          if (length(noc)==1 && noc>object$noc)
          stop("You cannot predict with more predictors than you have fitted")

          if (length(noc)==1 && noc<=object$noc)
          return(fitted(object))

          if (length(noc)>1 & max(noc)>object$noc)
          stop("You cannot predict with more predictors than you have fitted")

          if (length(noc)>1 & max(noc)<=object$noc)
          return(fitted(object)[, noc, drop = FALSE])
        }

        if (type == "class") {
          if (length(noc)==1 && noc>object$noc)
          stop("You cannot predict with more predictors than you have fitted")

          if (length(noc)==1 && noc<=object$noc)
            {
              xlearn <- (fitted(object))[,1:noc, drop = FALSE]
              xtest  <- (fitted(object))[,1:noc, drop = FALSE]
              return(switch(classifier,
                            nnr      = nnr(xlearn, xtest, object$y),
                            dlda     = dlda(xlearn, xtest, object$y),
                            logreg   = logreg(xlearn, xtest, object$y),
                            aggtrees = aggtrees(xlearn, xtest, object$y)))
            }

          if (length(noc)>1 & max(noc)>object$noc)
          stop("You cannot predict with more predictors than you have fitted")

          if (length(noc)>1 && max(noc)<=object$noc)
            {
              clmt <- NULL
              for (i in 1:length(noc))
                {
                  xl <- (fitted(object))[, 1:(noc[i]), drop = FALSE]
                  xt <- (fitted(object))[, 1:(noc[i]), drop = FALSE]
                  cl <- switch(classifier,
                               nnr      = nnr(xl, xt, object$y),
                               dlda     = dlda(xl, xt, object$y),
                               logreg   = logreg(xl, xt, object$y),
                               aggtrees = aggtrees(xl, xt, object$y))
                  cl           <- matrix(cl, ncol=1)
                  dimnames(cl) <- list(1:nrow(xl), paste(noc[i],"Predictors"))
                  clmt         <- cbind(clmt, cl)
                }
              return(clmt)
            }
        }
      }

    ## Returning fitted values or class labels for test data
    else
      {
        ## Check the dimensions of the new data
        if (!is.null(object$signs) && ncol(newdata)!=length(object$signs))
          {
            stop(paste("The new data need to have the same number of",
                       "predictor variables as the training data"))
          }

        ## Flip the signs of the new data
        if (!is.null(object$signs))
          {
            newdata <- t(t(newdata)*object$signs)
          }

        ## Fitted values for the new data
        xtest <- NULL
        for (j in 1:object$noc)
          {
            xtest <- cbind(xtest, rowMeans(newdata[,object$clist[[j]],
                                                   drop = FALSE]))
          }

        ## Naming the predictors and samples
        sampnames <- 1:nrow(newdata)
        clusnames <- character(object$noc)
        for(i in 1:object$noc)      clusnames[i] <- paste("Predictor", i)
        dimnames(xtest) <- list(sampnames, clusnames)

        ## Return the fitted values for the new data if requested
        if (type == "fitted") {
          if (length(noc)==1 && noc>object$noc)
          stop("You cannot predict with more predictors than you have fitted")

          if (length(noc)==1 && noc<=object$noc)
          return(xtest[,1:noc, drop = FALSE])

          if (length(noc)>1 & max(noc)>object$noc)
          stop("You cannot predict with more predictors than you have fitted")

          if (length(noc)>1 & max(noc)<=object$noc)
          return(xtest[,noc, drop = FALSE])
        }

        ## Do 0/1-classification
        if (length(noc)==1 && noc>object$noc)
            stop("You cannot predict with more predictors than you have fitted")

        if (length(noc)==1 && noc==object$noc)
          {
            xlearn <- fitted(object)
            return(switch(classifier,
                          nnr      = nnr(xlearn, xtest, object$y),
                          dlda     = dlda(xlearn, xtest, object$y),
                          logreg   = logreg(xlearn, xtest, object$y),
                          aggtrees = aggtrees(xlearn, xtest, object$y)))
          }

        if (length(noc)==1 && noc<object$noc)
          {
            xlearn <- (fitted(object))[, 1:noc, drop = FALSE]
            xtest  <- xtest[, 1:noc, drop = FALSE]
            return(switch(classifier,
                          nnr      = nnr(xlearn, xtest, object$y),
                          dlda     = dlda(xlearn, xtest, object$y),
                          logreg   = logreg(xlearn, xtest, object$y),
                          aggtrees = aggtrees(xlearn, xtest, object$y)))
          }

        if (length(noc)>1 && max(noc)>object$noc)
            stop("You cannot predict with more predictors than you have fitted")

        if (length(noc)>1 && max(noc)<=object$noc)
          {
            clmt <- NULL
            for (i in 1:length(noc))
              {
                xl <- (fitted(object))[, 1:(noc[i]), drop = FALSE]
                xt <- xtest[, 1:(noc[i]), drop = FALSE]
                cl <- switch(classifier,
                             nnr      = nnr(xl, xt, object$y),
                             dlda     = dlda(xl, xt, object$y),
                             logreg   = logreg(xl, xt, object$y),
                             aggtrees = aggtrees(xl, xt, object$y))
                cl           <- matrix(cl, ncol=1)
                dimnames(cl) <- list(1:nrow(newdata),paste(noc[i],"Predictors"))
                clmt         <- cbind(clmt, cl)
              }
            return(clmt)
          }
      }
  }


