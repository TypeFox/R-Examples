# Sobol' indices estimation (Martinez 2011)
# Plus: Theoretical confidence intervals from correlation coefficient-based confidence interval
#
# J-M. Martinez, Analyse de sensibilite globale par decomposition de la variance,
# Presentation a la journee des GdR Ondes et MASCOT-NUM, 13 janvier 2011, 
# Institut Henri Poincare, Paris, France.
#
# Bertrand Iooss (2015)
# Modified by Frank Weber (2016)


sobolmartinez <- function(model = NULL, X1, X2, nboot = 0, conf = 0.95, ...) {
  if ((ncol(X1) != ncol(X2)) | (nrow(X1) != nrow(X2)))
    stop("The samples X1 and X2 must have the same dimensions")
  p <- ncol(X1)
  
  X <- rbind(X1,X2)
  for (i in 1:p) {
    Xb <- X1
    Xb[,i] <- X2[,i]
    X <- rbind(X, Xb) 
  }
  
  x <- list(model = model, X1 = X1, X2 = X2, nboot = nboot, conf = conf, X = X,
            call = match.call())
  class(x) <- "sobolmartinez"
  
  if (!is.null(x$model)) {
    response(x, other_types_allowed = TRUE, ...)
    tell(x)
  }
  
  return(x)
}


estim.sobolmartinez <- function(data, i = NULL, estimStd = FALSE, conf = 0){
  if(class(data) == "matrix"){
    # This means x$y is a numeric vector.
    if(is.null(i)) i <- 1:nrow(data)
    d <- as.matrix(data[i, ]) # as.matrix for colSums
    n <- nrow(d)
    p <- ncol(d) - 2
    
    V <- var(d[, 1])
    ecor <- sapply(1:p, function(ii){
      cor(d[, 2], d[, ii + 2], use = "pairwise.complete.obs")
    })
    ecorcompl <- sapply(1:p, function(ii){
      cor(d[, 1], d[, ii + 2], use = "pairwise.complete.obs")
    })
    if(estimStd){
      VV <- matrix(V, nrow = 1, ncol = 3, 
                   dimnames = list(1, c("estim", "CIinf", "CIsup")))
      est_matrix <- sapply(1:p, function(ii){
        confcor <- cor.test(d[, 2],d[, ii + 2], conf.level = conf)
        estcor <- c(ecor[ii], 
                    confcor$conf.int[1], 
                    confcor$conf.int[2])
        confcor <- cor.test(d[, 1], d[, ii + 2], conf.level = conf)
        estcorcompl <- c(ecorcompl[ii], 
                         confcor$conf.int[2], 
                         confcor$conf.int[1]) # on intervertit car apres on prend l'oppose
        return(c(estcor, estcorcompl))
      })
      estcor <- t(est_matrix[1:3, ])
      estcorcompl <- t(est_matrix[4:6, ])
      dimnames(estcor) <- list(2:(p + 1), c("estim", "CIinf", "CIsup"))
      dimnames(estcorcompl) <- list((p + 2):(2*p + 1), 
                                    c("estim", "CIinf", "CIsup"))
      return(rbind(VV, estcor, estcorcompl))
    } else{
      return(c(V, ecor, ecorcompl))
    }
  } else if(class(data) == "array"){
    if(estimStd){
      stop("Confidence intervals not supported if \"data\" is an array")
    }
    if(is.null(i)) i <- 1:dim(data)[1]
    n <- length(i)
    p <- dim(data)[2] - 2
    
    # Define a helper function:
    one_dim3 <- function(d_array){
      V <- apply(d_array, 3, function(d_matrix){
        var(d_matrix[, 1])
      })
      ematrix <- apply(d_array, 3, function(d_matrix){
        ecor <- sapply(1:p, function(ii){
          cor(d_matrix[, 2], d_matrix[, ii + 2], use = "pairwise.complete.obs")
        })
        ecorcompl <- sapply(1:p, function(ii){
          cor(d_matrix[, 1], d_matrix[, ii + 2], use = "pairwise.complete.obs")
        })
        c(ecor, ecorcompl)
      })
      return(rbind(V, ematrix, deparse.level = 0))
    }
    
    if(length(dim(data)) == 3){
      # This means x$y is a matrix.
      d <- data[i, , , drop = FALSE]
      return(one_dim3(d))
    } else if(length(dim(data)) == 4){
      # This means x$y is a 3-dimensional array.
      d <- data[i, , , , drop = FALSE]
      all_dim3 <- sapply(1:dim(data)[4], function(i){
        one_dim3(array(data[ , , , i], 
                       dim = dim(data)[1:3], 
                       dimnames = dimnames(data)[1:3]))
      }, simplify = "array")
      dimnames(all_dim3)[[3]] <- dimnames(data)[[4]]
      return(all_dim3)
    }
  }
}


tell.sobolmartinez <- function(x, y = NULL, return.var = NULL, ...) {
  id <- deparse(substitute(x))
  
  if (! is.null(y)) {
    x$y <- y
  } else if (is.null(x$y)) {
    stop("y not found")
  }
  
  p <- ncol(x$X1)
  n <- nrow(x$X1)
  
  if(class(x$y) == "numeric"){
    data <- matrix(x$y, nrow = n)
    
    # estimation of the partial variances (V, D1 and Dt)
    if (x$nboot == 0){
      V <- data.frame(original = estim.sobolmartinez(data, 1:n, TRUE, x$conf))
      colnames(V) <- c("original", "min. c.i.", "max. c.i.")
    }
    else{
      V.boot <- boot(data, estim.sobolmartinez, R = x$nboot)
      V <- bootstats(V.boot, x$conf, "basic")
      rownames(V) <- c("global", 
                       colnames(x$X1), 
                       paste("-", colnames(x$X1), sep = ""))
    }
    
    # estimation of the Sobol' indices (S1 and St)
    if (x$nboot == 0) {
      S <- V[2:(p + 1), 1:3, drop = FALSE]
      T <- 1 - V[(p + 2):(2 * p + 1), 1:3, drop = FALSE]
      
    } else {
      S.boot <- V.boot
      S.boot$t0 <- V.boot$t0[2:(p + 1)]
      S.boot$t <- V.boot$t[,2:(p + 1)]
      S <- bootstats(S.boot, x$conf, "basic")
      
      T.boot <- V.boot
      T.boot$t0 <- 1 - V.boot$t0[(p + 2):(2 * p + 1)]
      T.boot$t <- 1 - V.boot$t[,(p + 2):(2 * p + 1)]
      T <- bootstats(T.boot, x$conf, "basic")
    }
    rownames(S) <- colnames(x$X1)
    rownames(T) <- colnames(x$X1)
  } else if(class(x$y) == "matrix"){
    data <- array(x$y, dim = c(n, nrow(x$y) / n, ncol(x$y)), 
                  dimnames = list(NULL, NULL, colnames(x$y)))
    if(x$nboot == 0){
      V <- estim.sobolmartinez(data, 1:n, estimStd = FALSE)
      rownames(V) <- c("global", 
                       colnames(x$X1), 
                       paste("-", colnames(x$X1), sep = ""))
      S <- V[2:(p + 1), , drop = FALSE]
      T <- 1 - V[(p + 2):(2 * p + 1), , drop = FALSE]
      rownames(T) <- colnames(x$X1)
    } else{
      V.boot <- lapply(1:ncol(x$y), function(col_idx){
        boot(as.matrix(data[, , col_idx]), estim.sobolmartinez, R = x$nboot)
      })
      V <- sapply(1:length(V.boot), function(col_idx){
        as.matrix(bootstats(V.boot[[col_idx]], x$conf, "basic"))
      }, simplify = "array")
      dimnames(V) <- list(
        c("global", colnames(x$X1), paste("-", colnames(x$X1), sep = "")),
        dimnames(V)[[2]],
        colnames(x$y))
      S <- sapply(1:length(V.boot), function(col_idx){
        S.boot_col <- V.boot[[col_idx]]
        S.boot_col$t0 <- V.boot[[col_idx]]$t0[2:(p + 1)]
        S.boot_col$t <- V.boot[[col_idx]]$t[, 2:(p + 1)]
        as.matrix(bootstats(S.boot_col, x$conf, "basic"))
      }, simplify = "array")
      T <- sapply(1:length(V.boot), function(col_idx){
        T.boot_col <- V.boot[[col_idx]]
        T.boot_col$t0 <- 1 - V.boot[[col_idx]]$t0[(p + 2):(2 * p + 1)]
        T.boot_col$t <- 1 - V.boot[[col_idx]]$t[, (p + 2):(2 * p + 1)]
        as.matrix(bootstats(T.boot_col, x$conf, "basic"))
      }, simplify = "array")
      dimnames(S) <- dimnames(T) <- list(colnames(x$X1),
                                         dimnames(V)[[2]],
                                         colnames(x$y))
    }
  } else if(class(x$y) == "array"){
    data <- array(x$y, dim = c(n, dim(x$y)[1] / n, dim(x$y)[2:3]), 
                  dimnames = list(NULL, NULL, 
                                  dimnames(x$y)[[2]], dimnames(x$y)[[3]]))
    if(x$nboot == 0){
      V <- estim.sobolmartinez(data, 1:n, estimStd = FALSE)
      dimnames(V)[[1]] <- c("global", 
                            colnames(x$X1), 
                            paste("-", colnames(x$X1), sep = ""))
      S <- V[2:(p + 1), , , drop = FALSE]
      T <- 1 - V[(p + 2):(2 * p + 1), , , drop = FALSE]
      dimnames(T)[[1]] <- colnames(x$X1)
    } else{
      V.boot <- lapply(1:dim(x$y)[[3]], function(dim3_idx){
        lapply(1:dim(x$y)[[2]], function(dim2_idx){
          boot(as.matrix(data[, , dim2_idx, dim3_idx]), estim.sobolmartinez, R = x$nboot)
        })
      })
      V <- sapply(1:dim(x$y)[[3]], function(dim3_idx){
        sapply(1:dim(x$y)[[2]], function(dim2_idx){
          as.matrix(bootstats(V.boot[[dim3_idx]][[dim2_idx]], x$conf, "basic"))
        }, simplify = "array")
      }, simplify = "array")
      dimnames(V) <- list(c("global", 
                            colnames(x$X1), 
                            paste("-", colnames(x$X1), sep = "")),
                          dimnames(V)[[2]],
                          dimnames(x$y)[[2]],
                          dimnames(x$y)[[3]])
      S <- sapply(1:dim(x$y)[[3]], function(dim3_idx){
        sapply(1:dim(x$y)[[2]], function(dim2_idx){
          S.boot_dim2 <- V.boot[[dim3_idx]][[dim2_idx]]
          S.boot_dim2$t0 <- 
            V.boot[[dim3_idx]][[dim2_idx]]$t0[2:(p + 1)]
          S.boot_dim2$t <- 
            V.boot[[dim3_idx]][[dim2_idx]]$t[, 2:(p + 1)]
          as.matrix(bootstats(S.boot_dim2, x$conf, "basic"))
        }, simplify = "array")
      }, simplify = "array")
      T <- sapply(1:dim(x$y)[[3]], function(dim3_idx){
        sapply(1:dim(x$y)[[2]], function(dim2_idx){
          T.boot_dim2 <- V.boot[[dim3_idx]][[dim2_idx]]
          T.boot_dim2$t0 <- 
            1 - V.boot[[dim3_idx]][[dim2_idx]]$t0[(p + 2):(2 * p + 1)]
          T.boot_dim2$t <- 
            1 - V.boot[[dim3_idx]][[dim2_idx]]$t[, (p + 2):(2 * p + 1)]
          as.matrix(bootstats(T.boot_dim2, x$conf, "basic"))
        }, simplify = "array")
      }, simplify = "array")
      dimnames(S) <- dimnames(T) <- list(colnames(x$X1),
                                         dimnames(V)[[2]],
                                         dimnames(x$y)[[2]],
                                         dimnames(x$y)[[3]])
    }
  }
  
  # return
  x$V <- V
  x$S <- S
  x$T <- T
  
  for (i in return.var) {
    x[[i]] <- get(i)
  }
  
  assign(id, x, parent.frame())
}


print.sobolmartinez <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (!is.null(x$y)) {
    if (class(x$y) == "numeric") {
      cat("\nModel runs:", length(x$y), "\n")
    } else if (class(x$y) == "matrix") {
      cat("\nModel runs:", nrow(x$y), "\n")
    } else if (class(x$y) == "array") {
      cat("\nModel runs:", dim(x$y)[1], "\n")
    }
    cat("\nFirst order indices:\n")
    print(x$S)
    cat("\nTotal indices:\n")
    print(x$T)
  } else {
    cat("\n(empty)\n")
  }
}


plot.sobolmartinez <- function(x, ylim = c(0, 1), 
                               y_col = NULL, y_dim3 = NULL, ...) {
  if (!is.null(x$y)) {
    p <- ncol(x$X1)
    pch = c(21, 24)
    if(class(x$y) == "numeric"){
      nodeplot(x$S, xlim = c(1, p + 1), ylim = ylim, pch = pch[1])
      nodeplot(x$T, xlim = c(1, p + 1), ylim = ylim, labels = FALSE,
               pch = pch[2], at = (1:p)+.3, add = TRUE)
    } else if(class(x$y) %in% c("matrix", "array")){
      if(is.null(y_col)) y_col <- 1
      if(class(x$y) == "matrix" && !is.null(y_dim3)){
        y_dim3 <- NULL
        warning("Argument \"y_dim3\" is ignored since the model output is ",
                "a matrix")
      }
      if(class(x$y) == "array" && is.null(y_dim3)) y_dim3 <- 1
      nodeplot(x$S, xlim = c(1, p + 1), ylim = ylim, pch = pch[1], 
               y_col = y_col, y_dim3 = y_dim3)
      nodeplot(x$T, xlim = c(1, p + 1), ylim = ylim, labels = FALSE,
               pch = pch[2], at = (1:p)+.3, add = TRUE, 
               y_col = y_col, y_dim3 = y_dim3)
    }
    legend(x = "topright", legend = c("main effect", "total effect"), pch = pch)
  }
}
