`B_function` <-
function(m,C,Gmatx,Wtx,mx,Cx)
  {
    Rx <- Gmatx %*% C %*% t(Gmatx) + Wtx# Gmatx = G_{t+1}
    B  <- C %*% t(Gmatx) %*% solve( Rx )
    ms <- t(m) + B%*%(t(mx) - Gmatx %*%t(m)) # mx=ms_{t+1}
    Cs <- C + B%*%(Cx-Rx)%*%t(B)
    list(ms=ms,Cs=Cs,B=B)
  }

