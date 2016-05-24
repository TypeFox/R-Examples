bayesx.reformate <- function(x)
{
  if(!is.null(x$fixed.effects) && length(x$fixed.effects) > 0L) {
    e <- x$fixed.effects
    type <- if(colnames(e)[1] == "pmode") "REML" else "MCMC"
    if(type == "REML") {
      tvalues <- e[,1L]/e[,4L]
      if(!is.null(x$model.fit$N)) {
        pvalues <- 2 * pt(-abs(tvalues), df = x$model.fit$N - 1)
        e <- cbind(e[,1L], e[,4L], tvalues, pvalues)
        if(!is.matrix(e))
          e <- matrix(e, nrow = 1L)
        colnames(e) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
      } else {
        e <- cbind(e[,1L], e[,4L], tvalues)
        if(!is.matrix(e))
          e <- matrix(e, nrow = 1L)
        colnames(e) <- c("Estimate", "Std. Error", "t value")
      }
    }
    if(type == "MCMC") {
      e <- e[, c(1L, 2L, 3L, 5L, 7L)]
      if(!is.matrix(e))
        e <- matrix(e, nrow = 1L)
      colnames(e) <- c("Mean", "Sd", "2.5%", "50%", "97.5%")
    }
    eattrn <- names(attributes(x$fixed.effects))
    for(i in 1L:length(eattrn))
      if(eattrn[i] != "dim" && eattrn[i] != "dimnames")
        attr(e, eattrn[i]) <- attr(x$fixed.effects, eattrn[i])
    rownames(e) <- rownames(x$fixed.effects)
    if(!is.null(attr(e, "sample")))
      if(ncol(attr(e, "sample")) == length(rownames(e)))
        colnames(attr(e, "sample")) <- rownames(e)
    x$fixed.effects <- e
  }
  if(!is.null(x$effects) && length(x$effects) > 0L) {
    for(i in 1L:length(x$effects)) {
      if(!is.null(x$effects[[i]])) {
        if(is.list(x$effects[[i]])) {
          if(length(x$effects[[i]]) > 0L)
            for(j in 1L:length(x$effects[[i]]))
              x$effects[[i]][[j]] <- chacol(x$effects[[i]][[j]])
        } else x$effects[[i]] <- chacol(x$effects[[i]])
      }
    }
  }
  if(!is.null(x$smooth.hyp) && length(x$smooth.hyp) > 0L)
    x$smooth.hyp <- recol(chacol(x$smooth.hyp))
  if(!is.null(x$random.hyp) && length(x$random.hyp) > 0L)
    x$random.hyp <- recol(chacol(x$random.hyp))
  if(!is.null(x$variance) && length(x$variance) > 0L)
    x$variance <- recol(chacol(x$variance))

  return(x)
}

