new.lambda <- function(penden.env,lambda0) {
  eps <- 0.001*lambda0
  epsdf <- 0.01
  M <- get("M",penden.env)
  N <- get("N",penden.env)
  help <- c(((M-1)*N+1):(M*N))
  p <- get("p",penden.env)
  m <- get("m",penden.env)
  beta.val.help <- c(t(get("beta.val",penden.env)))[-help]
  
  calc <- TRUE
  lambda <- lambda0
  u <- t(beta.val.help)%*%get("Dm",penden.env)%*%beta.val.help
  hh <-1
  while(calc) {
    if(hh==11) break
    
    Derv2.tmp <- Derv2(penden.env,lambda[hh])
    df <- sum(diag(my.positive.definite.solve(Derv2.tmp$Derv2.pen)%*%Derv2.tmp$Derv2.cal))

    if(df-(p*(m-1)) < epsdf) break
      #{print("df kleiner m-1+epsdf")
     
    
    df<- df-(p*(m-1))

    help <- abs(df/u - lambda[hh])
    if((df/u)<0) break
    if(help<eps) calc <- FALSE else {
      lambda[hh+1] <- df/u
      hh <- hh+1
    }
  }
  return(lambda[hh])
}
