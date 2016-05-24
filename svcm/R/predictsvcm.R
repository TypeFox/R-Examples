predictsvcm <- function(lambda, type = type){

  if (type == "TP") {
    ##build penalty matrix
    if (length(lambda) == 1) {
      P <- lambda * P.whole
    } else if (length(lambda) == 2 & ndims == 2){
      P <- as(lambda[1] * P.x + lambda[2] * P.y, "dsCMatrix")
    } else if (length(lambda) == 3 & ndims == 3){
      P <- as(lambda[1] * P.x + lambda[2] * P.y + lambda[3] * P.z, "dsCMatrix")
    } else {
      stop(paste("\nUse either global lambda or dimension-specific lambda",
                 "of length", ndims, "!\n"))
    }
    
    ##compute sandwich matrix
    S <- as(BB + P, "dsCMatrix")
    assign("S", S, env = svcm.env)
    ##compute amplitudes (=coefficients)
    A <- solve(S, RHS)
    assign("A", as.vector(A), env = svcm.env)
    ##compute effects: could be replaced by vec(XAB') = (B %x% X)vec(A),
    ##                 exploiting tensor product structure.
    eta <- B %*% A
    assign("eta", eta, env = svcm.env)
    ##compute trace of hat matrix, i.e. effective dimension
    EDofTP()
    
  } else if (type == "SEQ") {

    if (length(lambda) == 1) lambda <- rep(lambda, ndims)
    
    ##compute amplitudes (=coefficients)
    ##compute trace of hat matrix = kronecker of all B_i S_i^{-1} B_i',i=x,y,z;
    ##Note that the trace of Kronecker products is the product of the traces.
    A <- SEQpsi(as.matrix(LS.X), Y, 1)
    ED <- sum(diag(X %*% LS.X))
    LS <- solve(as.matrix(BB.x + lambda[1] * P.x)) %*% t(B.x)
    A <- SEQpsi(as.matrix(LS), A, 2)
    ED <- ED * sum(diag(B.x %*% LS))
    LS <- solve(as.matrix(BB.y + lambda[2] * P.y)) %*% t(B.y)
    A <- SEQpsi(as.matrix(LS), A, 3)
    ED <- ED * sum(diag(B.y %*% LS))
    if (ndims == 3) {
      LS <- solve(as.matrix(BB.z + lambda[3] * P.z)) %*% t(B.z)
      A <- SEQpsi(as.matrix(LS), A, 4)
      ED <- ED * sum(diag(B.z %*% LS))
    }
    assign("A", A, env = svcm.env)
    assign("ED", ED, env = svcm.env)
#    cat(paste("Effective dimension =", ED, ".\n"))
    
    ##compute effects
    eta <- SEQpsi(as.matrix(X), A, 1)
    eta <- SEQpsi(as.matrix(B.x), eta, 2)
    eta <- SEQpsi(as.matrix(B.y), eta, 3)
    if (ndims == 3) {
      eta <- SEQpsi(as.matrix(B.z), eta, 4)
    }
    assign("eta", eta, env = svcm.env)

  }
  
}
