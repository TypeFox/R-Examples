plot.rdacv <- function(x, type=c("both", "error", "gene"), 
                       nice=FALSE, ...){

  if(class(x) != 'rdacv') {
    stop("You must supply a cross-validation object.")
  }
  else{
    minerr <- min(x$cv.err)/x$n
    one.se <- ceiling((sqrt(minerr*(1-minerr)/x$n)*1.645+minerr)*x$n)
    pos <- which(x$cv.err == one.se, TRUE)
    minpos <- which(x$cv.err == min(x$cv.err), TRUE)

    type <- match.arg(type)
    switch(type, 
    both={
      tmperr <- x$cv.err
      dimnames(tmperr) <- list(x$alpha, x$delta)
      tmpgene <- x$ngene
      dimnames(tmpgene) <- list(x$alpha, x$delta)

      par(ask=TRUE)
      rda.plotmat(tmperr, se=TRUE,
                  main=paste("Heatmap of ", x$nfold,
                             "-fold CV Error"),
                  pos=pos, minpos=minpos, nice=nice)
      rda.plotmat(tmpgene, se=TRUE,
                  main=paste("Heatmap of Number of Genes Remained"),
                  pos=pos, minpos=minpos, nice=nice)
      par(ask=FALSE)
      return(list(one.se.pos=pos, min.cv.pos=minpos))
    },
    error={
      tmperr <- x$cv.err
      dimnames(tmperr) <- list(x$alpha, x$delta)
      rda.plotmat(tmperr, se=TRUE,
                  main=paste("Heatmap of ", x$nfold, "-fold CV Error"),
                  pos=pos, minpos=minpos, nice=nice)
      return(list(one.se.pos=pos, min.cv.pos=minpos))
    },
    gene={
      tmpgene <- x$ngene
      dimnames(tmpgene) <- list(x$alpha, x$delta)
      rda.plotmat(tmpgene, se=TRUE,
                  main=paste("Heatmap of Number of Genes Remained"),
                  pos=pos, minpos=minpos, nice=nice)
      return(list(one.se.pos=pos, min.cv.pos=minpos))
    })
  }
}

