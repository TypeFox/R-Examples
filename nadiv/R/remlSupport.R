###########################################################

# Code for the pin function is taken and modified from R code
# with the same name that was created and posted online by 
# Ian White at the University of Edinburgh

##########################################################
pin <- function (object, transform){
  pframe <- as.list(object$gammas)
  names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
  tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)), pframe)
  X <- as.vector(attr(tvalue, "gradient"))
  X[object$gammas.type == 1] <- 0
  tname <- if (length(transform) == 3) transform[[2]] else ""
  n <- length(pframe)
  i <- rep(1:n, 1:n)
  j <- sequence(1:n)
  k <- 1 + (i > j)
  Vmat <- object$ai
  se <- sqrt(sum(Vmat * X[i] * X[j] * k))
 return(data.frame(row.names = tname, Estimate = tvalue, SE = se))
}


pcc <- function(object, traces = NULL, tol = 0.01, silent = FALSE){
    if(is.null(object) & is.null(traces)){
       stop("one of object or traces must be non-NULL")
    }

    if(!is.null(object)){
       if(!silent) message("trimming the asreml monitor matrix")
       rc <- dim(object$monitor)
       traces <- object$monitor[4:rc[1], 1:(rc[2]-1)]
       if(object$loglik == 0.00 | object$converge == FALSE) return(FALSE)
    } 
    rc <- dim(traces)
    penultimate <- rc[2] - 1
    pchange <- abs(traces[, rc[2]] - traces[, penultimate]) / c(apply(traces[, penultimate:rc[2]], MARGIN = 1, FUN = max))
     return(all(pchange < tol))
}

   
