PolarDecomp<-function(A)
  {
    ### A =U %*% P  is the polar decomposition
    ###      W         S           V
    ###  A = E$u %*% diag(E$d)%*%  t(E$v)
    ###  U is orthogonal rotation matrix
    ###  P is the stretch tensor
    E = svd(A)
    P = E$v%*% diag(E$d)%*%  t(E$v)
    U =  E$u %*% t(E$v)

    return(list(P=P, U=U))

  }
