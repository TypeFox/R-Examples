enet.selection.A.b.from.enet <-
function(enet.fit){
  # unpack the required info
  X = enet.fit$X
  n = nrow(X)
  beta = enet.fit$beta
  lambda = enet.fit$lambda
  alpha = enet.fit$alpha
  mu = enet.fit$mu
  
  E = which (abs(beta) > 10^-8) # selected set
  #print(E)
  
  if (length (E) > 0){
    X.E = X[,E, drop=FALSE]
    X.Emin = X[,-E, drop=FALSE]
    s.E = sign (beta[E])
    
    if (length(E) == 1){
      D.s.E = matrix (s.E, nrow=1, ncol=1)
    }else{
      D.s.E = diag(s.E)
    }
    
    # compute some of the important matrices
    Q = solve (t(X.E)%*%X.E + lambda*(1-alpha)*diag(length(E)), t(X.E))
    P = X.E%*%Q
    
    # compute A
    A = rbind ( (t(X.Emin) - t(X.Emin)%*%P)/lambda/alpha,
                -(t(X.Emin) - t(X.Emin)%*%P)/lambda/alpha,
                -D.s.E%*%Q
    )
    
    # compute b
    b = c( 1 - t(X.Emin)%*%t(Q)%*%s.E,
           1 + t(X.Emin)%*%t(Q)%*%s.E,
           -(alpha*lambda)*D.s.E%*%solve(t(X.E)%*%X.E + lambda*(1-alpha)*diag(length(E)), s.E)
    )
    
    # check whether these matrices came from specified or unspecified mu and adjust accordingly
    if (is.null(mu)){ # had to estimate mu, so adjust A
      A = A - apply (A, 1, mean)%*%t(rep(1, ncol(A)))
    }else{ # first subtracted mu from y, so adjust b vector
      b = b + mu*apply(A, 1, sum)
    }
    
    return(list (which.col=E, A=A, b=b))
  }else{ # no variables selected
    X.Emin = X
    
    # compute A
    A = rbind ( t(X.Emin)/lambda/alpha,
                -t(X.Emin)/lambda/alpha
    )
    
    # compute b
    b = c( rep(1, ncol(X.Emin)),
           rep(1, ncol(X.Emin))
    )
    
    # check whether these matrices came from specified or unspecified mu and adjust accordingly
    if (is.null(mu)){ # had to estimate mu, so adjust A
      A = A - apply (A, 1, sum)%*%t(rep(1, ncol(A)))/n
    }else{ # first subtracted mu from y, so adjust b vector
      b = b + mu*apply(A, 1, sum)
    }
    
    return(list (which.col=E, A=A, b=b))
  }
}
