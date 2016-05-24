# This provides a projection on 2d.

project.on.2d <- function(x, emobj = NULL, pi = NULL, Mu = NULL,
    LTSigma = NULL, class = NULL, method = c("PP", "SVD")){
  if(is.null(emobj)){
    emobj <- list(pi = pi, Mu = Mu, LTSigma = LTSigma, class = class)
  }
  var <- LTSigma2variance(emobj$LTSigma)

  if(method[1] == "PP"){
    ### Convert original S to spheric condition.
    x.pp <- x
    Sigma.pp <- var
    for(k.var in 1:dim(var)[3]){
      tmp <- eigen(var[,, k.var])
      Sigma.k.inv <- tmp$vector %*% diag(sqrt(1/tmp$values)) %*% t(tmp$vector)

      tmp <- x[emobj$class == k.var,]
      tmp.mu <- colMeans(tmp)
      tmp <- t(t(tmp) - tmp.mu) %*% Sigma.k.inv
      tmp <- t(t(tmp) + as.vector(tmp.mu))
      x.pp[emobj$class == k.var,] <- tmp
    }

    ### Use PPtree to get a better view and a projection matrix.
    tmp.pp <- PPtree::PP.optimize.random("LDA", 2, x.pp, emobj$class,
                                         std = FALSE)
    proj.mat <- tmp.pp$proj.best
  } else if(method[1] == "SVD"){
    ### Obtain a projection matrix based on the largest two components.
    x.svd <- svd(x)
    proj.mat <- diag(sqrt(x.svd$d)) %*% x.svd$v[, 1:2]
  } else{
    stop("method is not found.")
  }

  ### Project to 2D and convert S back to elipsoid condition.
  x.new <- x %*% proj.mat
  mu.new <- emobj$Mu %*% proj.mat
  var.new <- array(0, dim = c(2, 2, dim(var)[3]))
  for(k.var in 1:dim(var)[3]){
    tmp <- t(proj.mat) %*% Sigma.pp[,, k.var] %*% proj.mat
    var.new[,, k.var] <- tmp
  }

  ### Return.
  ret <- list(da = x.new,
              Pi = emobj$pi, Mu = mu.new, S = var.new,
              class = emobj$class,
              proj.mat = proj.mat)
  ret
} # End of project.on.2d().

