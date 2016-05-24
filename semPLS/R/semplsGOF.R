dgrho <- function(object, ...){
  UseMethod("dgrho", object)
}

# Dillon-Goldstein's rho (Composite Reliability in SmartPLS)
# requires: outer loadings (factor scores), model
dgrho.sempls <- function(object, ...){
    dgr <- matrix(NA, nrow=length(object$model$latent), ncol=2)
    rownames(dgr) <- object$model$latent
    colnames(dgr) <- c("Dillon-Goldstein's rho", "reflective MVs")
    for(i in object$model$latent){
        if(attr(object$model$blocks[[i]], "mode")=="B"){
            next
        }
        x <- object$outer_loadings[, i]
        ind <- which(x!=0)
        if(length(ind)==1){
            dgr[i,2] <- 1
            next
        }
        else {
            x <- x[ind]
            dgr[i,1] <- sum(x)^2 / (sum(x)^2 + sum(1-x^2))
            dgr[i,2] <- length(ind)
        }
    }
    class(dgr) <- c("dgrho", class(dgr))
    return(dgr)
}

print.dgrho <- function(x, na.print=".", digits=2, ...){
  xChar <- format(as.data.frame(x), digits=digits, ...)
  xChar[is.na(x)] <- na.print
  print(xChar, ...)
  invisible(x) 
}

communality <- function(object, ...){
  UseMethod("communality", object)
}

# requires: outer loadings (factor scores), model
communality.sempls <- function(object, ...){
    com <- matrix(NA, nrow=length(object$model$latent), ncol=2)
    rownames(com) <- object$model$latent
    colnames(com) <- c("communality", "reflective MVs")
    for(i in object$model$latent){
        if(attr(object$model$blocks[[i]], "mode")=="B"){
            next
        }
        x <- object$outer_loadings[, i]
        ind <- which(x!=0)
        if(length(ind)==1){
            com[i,2] <- 1
            next
        }
        else {
            x <- x[ind]
            com[i,1] <- 1/length(x)*sum(x^2)
            com[i,2] <- length(ind)
        }
    }
    class(com) <- c("communality", class(com))
    return(com)
}

print.communality <- function(x, na.print=".", digits=2, ...){
  #xChar <- format(as.data.frame(unclass(x)), digits=digits, ...)
  xChar <- format(as.data.frame(x), digits=digits, ...)
  xChar[is.na(x)] <- na.print
  print(xChar, ...)
  aveCom <- sum(x[!is.na(x[,1] ),2], na.rm=TRUE)^-1 * sum(x[,1] * x[,2], na.rm=TRUE)
  cat(paste("\n\tAverage communality:", signif(aveCom, digits=digits), "\n"))
  invisible(x) 
}


redundancy <- function(object, ...){
  UseMethod("redundancy", object)
}
# Redundancy Example:
# requires: rSquared(predict, factor scores), communality(outer loadings (factor scores), model)
redundancy.sempls <- function(object, ...){
    red <- as.matrix(communality(object)[,1] * rSquared(object)[,1])
    colnames(red) <- "redundancy"
    class(red) <- c("redundancy", class(red))
    return(red)
}

print.redundancy <- function(x, na.print=".", digits=2, ...){
  print.table(x, na.print=na.print, digits=digits, ...)
  ## aveRed <- nrow(x)^-1 * sum(x[,1], na.rm=TRUE)
  aveRed <- mean(x[,1], na.rm=TRUE)
  cat(paste("\n\tAverage redundancy:", signif(aveRed, digits=digits), "\n"))
  invisible(x)
}


rSquared2 <- function(object, ...){
  UseMethod("rSquared2", object)
}
# requires: rSquared(predict, factor scores), communality(outer loadings (factor scores), model)
rSquared2.sempls <- function(object, na.rm=FALSE, ...){
  Y_hat <- predict(object, what="LVs", ...)
  if(sum(is.na(Y_hat)) > 0 & !na.rm) stop("Use argument 'na.rm=TRUE'!")
  R_squared <- apply(Y_hat, 2, var, na.rm=na.rm) / apply(object$factor_scores, 2, var, na.rm=na.rm)
  R_squared <- as.matrix(R_squared)
  R_squared <- cbind(R_squared, NA,colSums(object$model$D))
  colnames(R_squared) <- c("R-squared", "R-squared-corrected", "predecessors")
  R_squared[R_squared[,"predecessors"]==0, "R-squared"] <- NA
  ## correction
  ## correct <- function(rSqrd, J, N) {rSqrd - J*(1-rSqrd)/(N-J-1)}
  ## Fixed(2012-09-21): since there is no intercept, see summary.lm()
  correct <- function(rSqrd, J, N) {rSqrd - J*(1-rSqrd)/(N-J)}
  N <- object$N
  J <- R_squared[, "predecessors"]
  R_squared[, "R-squared-corrected"] <- correct(R_squared[, "R-squared"], J, N)
  R_squared <- as.data.frame(R_squared)
  class(R_squared) <- c("rSquared2", class(R_squared))
  return(R_squared)
}

print.rSquared2 <- function(x, na.print=".", digits=2, ...){
  ## xChar <- format(as.data.frame(unclass(x)), digits=digits, ...)
  xChar <- format(as.data.frame(x), digits=digits, ...)
  xChar[is.na(x)] <- na.print
  print(xChar)
  ## aveRsquared <- nrow(x)^-1 * sum(x[,1], na.rm=TRUE)
  ## see: PLS-Handbook, p. 58 (... where J is the total number of
  ## endogenous latent variables in the model.) 
  aveRsquared <- mean(x[,1], na.rm=TRUE)
  cat(paste("\n\tAverage R-squared:", signif(aveRsquared, digits=digits), "\n"))
  invisible(x)
}


gof <- function(object, ...){
  UseMethod("gof", object)
}

# requires: rSquared, communality
gof.sempls <- function(object, ...){
    rSq <- rSquared(object)
    ## aveRsq <- nrow(rSq)^-1 * sum(rSq[,1], na.rm=TRUE)
    aveRsq <- mean(rSq[,1], na.rm=TRUE)
    com <- communality(object)
    aveCom <- sum(com[!is.na(com[,1]), 2], na.rm=TRUE)^-1 *
      sum(com[,1] * com[,2], na.rm=TRUE)
    gof <- sqrt(aveCom * aveRsq)
    gof <- matrix(c(aveRsq, aveCom, gof), nrow=3, ncol=1)
    rownames(gof) <- c("Average R-squared", "Average Communality", "GoF")
    colnames(gof) <- c("Value")
    class(gof) <- c("gof", class(gof))
    return(gof)
}

print.gof <- function(x, na.print=".", digits=2, ...){
  print.table(as.matrix(x), na.print=na.print, digits=digits, ...)
  invisible(x)
}
  
