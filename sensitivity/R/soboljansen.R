# Sobol' indices estimation (Jansen 1999 - Saltelli 2010)
#
# Author: Bertrand Iooss 2012
# Modified by Frank Weber (2016)


soboljansen <- function(model = NULL, X1, X2, nboot = 0, conf = 0.95, ...) {
  if ((ncol(X1) != ncol(X2)) | (nrow(X1) != nrow(X2)))
    stop("The samples X1 and X2 must have the same dimensions")
  p <- ncol(X1)
  
  X <- rbind(X1, X2)
  for (i in 1:p) {
    Xb <- X1
    Xb[,i] <- X2[,i]
    X <- rbind(X, Xb) 
  }
  
  x <- list(model = model, X1 = X1, X2 = X2, nboot = nboot, conf = conf, X = X,
            call = match.call())
  class(x) <- "soboljansen"
  
  if (!is.null(x$model)) {
    response(x, other_types_allowed = TRUE, ...)
    tell(x)
  }
  
  return(x)
}


estim.soboljansen <- function(data, i = NULL) {
  if(class(data) == "matrix"){
    # This means x$y is a numeric vector.
    if(is.null(i)) i <- 1:nrow(data)
    d <- as.matrix(data[i, ]) # as.matrix for colSums
    n <- nrow(d)
    V <- var(d[, 1])
    VCE <- V - (colSums((d[,2] - d[, - c(1, 2)])^2) / (2 * n - 1))
    VCE.compl <- (colSums((d[,1] - d[, - c(1, 2)])^2) / (2 * n - 1))
    return(c(V, VCE, VCE.compl))
  } else if(class(data) == "array"){
    if(is.null(i)) i <- 1:dim(data)[1]
    n <- length(i)
    p <- dim(data)[2] - 2
    
    # Define a helper function:
    one_dim3 <- function(d_array){
      V <- apply(d_array, 3, function(d_matrix){
        var(d_matrix[, 1])
      })
      SumSq <- apply(d_array, 3, function(d_matrix){
        matrix(c(
          colSums((d_matrix[, 2] - d_matrix[, -c(1, 2), drop = FALSE])^2) / 
            (2 * n - 1),
          colSums((d_matrix[, 1] - d_matrix[, -c(1, 2), drop = FALSE])^2) / 
            (2 * n - 1)), 
          ncol = 2)
      })
      V_rep <- matrix(rep(V, each = p), ncol = dim(d_array)[3], 
                      dimnames = list(NULL, dimnames(d_array)[[3]]))
      VCE <- V_rep - SumSq[1:p, , drop = FALSE]
      VCE.compl <- SumSq[(p + 1):(2 * p), , drop = FALSE]
      return(rbind(V, VCE, VCE.compl, deparse.level = 0))
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


tell.soboljansen <- function(x, y = NULL, return.var = NULL, ...) {
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
      V <- data.frame(original = estim.soboljansen(data))
    } else{
      V.boot <- boot(data, estim.soboljansen, R = x$nboot)
      V <- bootstats(V.boot, x$conf, "basic")
    }
    rownames(V) <- c("global", 
                     colnames(x$X1), 
                     paste("-", colnames(x$X1), sep = ""))
    
    # estimation of the Sobol' indices (S1 and St)
    if (x$nboot == 0) {
      S <- V[2:(p + 1), 1, drop = FALSE] / V[1,1]
      T <- V[(p + 2):(2 * p + 1), 1, drop = FALSE] / V[1,1]
      rownames(T) <- colnames(x$X1)
    } else {
      S.boot <- V.boot
      S.boot$t0 <- V.boot$t0[2:(p + 1)] / V.boot$t0[1]
      S.boot$t <- V.boot$t[,2:(p + 1)] / V.boot$t[,1]
      S <- bootstats(S.boot, x$conf, "basic")
      
      T.boot <- V.boot
      T.boot$t0 <- V.boot$t0[(p + 2):(2 * p + 1)] / V.boot$t0[1]
      T.boot$t <- V.boot$t[,(p + 2):(2 * p + 1)] / V.boot$t[,1]
      T <- bootstats(T.boot, x$conf, "basic")
      rownames(S) <- colnames(x$X1)
      rownames(T) <- colnames(x$X1)
    }
  } else if(class(x$y) == "matrix"){
    data <- array(x$y, dim = c(n, nrow(x$y) / n, ncol(x$y)), 
                  dimnames = list(NULL, NULL, colnames(x$y)))
    if(x$nboot == 0){
      V <- estim.soboljansen(data)
      rownames(V) <- c("global", 
                       colnames(x$X1), 
                       paste("-", colnames(x$X1), sep = ""))
      V_global <- matrix(rep(V[1, ], p), nrow = p, byrow = TRUE)
      S <- V[2:(p + 1), , drop = FALSE] / V_global
      T <- V[(p + 2):(2 * p + 1), , drop = FALSE] / V_global
      rownames(T) <- colnames(x$X1)
    } else{
      V.boot <- lapply(1:ncol(x$y), function(col_idx){
        boot(as.matrix(data[, , col_idx]), estim.soboljansen, R = x$nboot)
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
        S.boot_col$t0 <- V.boot[[col_idx]]$t0[2:(p + 1)] / V.boot[[col_idx]]$t0[1]
        S.boot_col$t <- V.boot[[col_idx]]$t[, 2:(p + 1)] / V.boot[[col_idx]]$t[, 1]
        as.matrix(bootstats(S.boot_col, x$conf, "basic"))
      }, simplify = "array")
      T <- sapply(1:length(V.boot), function(col_idx){
        T.boot_col <- V.boot[[col_idx]]
        T.boot_col$t0 <- V.boot[[col_idx]]$t0[(p + 2):(2 * p + 1)] / V.boot[[col_idx]]$t0[1]
        T.boot_col$t <- V.boot[[col_idx]]$t[, (p + 2):(2 * p + 1)] / V.boot[[col_idx]]$t[, 1]
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
      V <- estim.soboljansen(data)
      dimnames(V)[[1]] <- c("global", 
                            colnames(x$X1), 
                            paste("-", colnames(x$X1), sep = ""))
      V_global <- array(rep(V[1, , ], each = p), dim = c(p, dim(x$y)[2:3]))
      S <- V[2:(p + 1), , , drop = FALSE] / V_global
      T <- V[(p + 2):(2 * p + 1), , , drop = FALSE] / V_global
      dimnames(T)[[1]] <- colnames(x$X1)
    } else{
      V.boot <- lapply(1:dim(x$y)[[3]], function(dim3_idx){
        lapply(1:dim(x$y)[[2]], function(dim2_idx){
          boot(as.matrix(data[, , dim2_idx, dim3_idx]), estim.soboljansen, R = x$nboot)
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
            V.boot[[dim3_idx]][[dim2_idx]]$t0[2:(p + 1)] / 
            V.boot[[dim3_idx]][[dim2_idx]]$t0[1]
          S.boot_dim2$t <- 
            V.boot[[dim3_idx]][[dim2_idx]]$t[, 2:(p + 1)] / 
            V.boot[[dim3_idx]][[dim2_idx]]$t[, 1]
          as.matrix(bootstats(S.boot_dim2, x$conf, "basic"))
        }, simplify = "array")
      }, simplify = "array")
      T <- sapply(1:dim(x$y)[[3]], function(dim3_idx){
        sapply(1:dim(x$y)[[2]], function(dim2_idx){
          T.boot_dim2 <- V.boot[[dim3_idx]][[dim2_idx]]
          T.boot_dim2$t0 <- 
            V.boot[[dim3_idx]][[dim2_idx]]$t0[(p + 2):(2 * p + 1)] / 
            V.boot[[dim3_idx]][[dim2_idx]]$t0[1]
          T.boot_dim2$t <- 
            V.boot[[dim3_idx]][[dim2_idx]]$t[, (p + 2):(2 * p + 1)] / 
            V.boot[[dim3_idx]][[dim2_idx]]$t[, 1]
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


print.soboljansen <- function(x, ...) {
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


plot.soboljansen <- function(x, ylim = c(0, 1), 
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
