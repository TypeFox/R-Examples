
sim.data <-
  function(model=c("ar1","ar2"),time=time,n.obs=n.obs, n.var=n.var,prob0=NULL,
  network=c("random","scale-free","hub","user_defined"),
  prec=NULL,gamma1=NULL,gamma2=NULL){
  model=match.arg(model)
  network=match.arg(network)
  t=time
  n=n.obs
  d=n.var
  r= round(runif(1),4)*10000
  if(model=="ar1") {
    if(network=="random") {
      L = sugm.generator(n=n,d=d,graph="random", prob=prob0, seed=1234+r, vis = FALSE)
      LL = sugm.generator(n=n,d=d,graph="random", prob=prob0, seed=4567+r, vis = FALSE)
      LLL = sugm.generator(n=n,d=d,graph="random", prob=prob0, seed=1564+r, vis = FALSE)

    }
    else if(network=="scale-free") {
      L = sugm.generator(n=n,d=d,graph="scale-free", prob=prob0, seed=1234+r, vis = FALSE)
      LL = sugm.generator(n=n,d=d,graph="scale-free", prob=prob0, seed=4567+r, vis = FALSE)
      LLL = sugm.generator(n=n,d=d,graph="scale-free", prob=prob0, seed=1564+r, vis = FALSE)

    }
    else if(network=="hub") {
     L = sugm.generator(n=n,d=d,graph="hub", prob=prob0, seed=1234+r, vis = FALSE)
     LL = sugm.generator(n=n,d=d,graph="hub", prob=prob0, seed=4567+r, vis = FALSE)
     LLL = sugm.generator(n=n,d=d,graph="hub", prob=prob0, seed=1564+r, vis = FALSE)

   }
   if(network=="user_defined"){
       mu <- rep(0,d)
       true_theta <- prec
       sigma1 <- solve(prec)
       true_gamma <- gamma1
       B <- true_gamma
   }
  else{
  true_theta = as.matrix(L$omega*L$theta)
  diag(true_theta)=1     ##theta is the precision matrix
  sigma1 <- L$sigma
  mu <- rep(0,d)
  thetaL = as.matrix(LL$omega*LL$theta)
  thetaLL = as.matrix(LLL$omega*LLL$theta)
  lwt = thetaL*(1*lower.tri(thetaL, diag =FALSE))
  upt = thetaLL*(1*upper.tri(thetaLL, diag = FALSE))
  B=upt+lwt
  ua=rbinom(d,1,0.3)
  uu=runif(d,0,1)
  uau = ua*uu
  diag(B)= uau
   B11 = B
 for(i in 1:d){
   for(j in 1:d){
    if(B[i,j] != 0)
       { m = rbinom(1,1,0.6)
         if(m ==0) B11[i,j] = -B[i,j]
       }

  }}
  true_gamma=B11

  }        #Gamma is the autoregressive coefficient matrix
  ##data generation
  xtn <-array(NA,c(t,d,n))
  xtt <- array(NA,c(t,d,1))
  for(i in 1:n){
    x0 <- rmvnorm(1, mu, sigma1, method="svd")
    for(j in 1:t){
      et <- rmvnorm(1, mu, sigma1, method="svd")
      xt <-  x0 %*% B11 + et
      xtt[j,,] <- xt
      x0 <- xt
   }
   xtn[,,i] <- round(xtt,3)
 }
 xy=matrix(aperm(xtn, c(3,1,2)), ncol=d)
 data1 <- as.longitudinal(xy, repeats=n)
 return(list(data1=data1,theta=true_theta, gamma=true_gamma))
 }
 if(model=="ar2") {
  if(network=="random") {
    L = sugm.generator(n=n,d=d,graph="random", prob=prob0, seed=12346+r, vis = FALSE)
    LL = sugm.generator(n=n,d=d,graph="random", prob=prob0, seed=45678+r, vis = FALSE)
    LLL = sugm.generator(n=n,d=d,graph="random", prob=prob0, seed=43219+r, vis = FALSE)
    LL1 = sugm.generator(n=n,d=d,graph="random", prob=prob0, seed=14578+r, vis = FALSE)
    LLL1 = sugm.generator(n=n,d=d,graph="random", prob=prob0, seed=96879+r, vis = FALSE)
  }
  else if(network=="scale-free") {
   L = sugm.generator(n=n,d=d,graph="scale-free", prob=prob0, seed=12346+r, vis = FALSE)
   LL = sugm.generator(n=n,d=d,graph="scale-free", prob=prob0, seed=45678+r, vis = FALSE)
   LLL = sugm.generator(n=n,d=d,graph="scale-free", prob=prob0, seed=43219+r, vis = FALSE)
   LL1 = sugm.generator(n=n,d=d,graph="scale-free", prob=prob0, seed=14578+r, vis = FALSE)
   LLL1 = sugm.generator(n=n,d=d,graph="scale-free", prob=prob0, seed=96879+r, vis = FALSE)
  }
 else if(network=="hub") {
  L = sugm.generator(n=n,d=d,graph="hub", prob=prob0, seed=12346+r, vis = FALSE)
  LL = sugm.generator(n=n,d=d,graph="hub", prob=prob0, seed=45678+r, vis = FALSE)
  LLL = sugm.generator(n=n,d=d,graph="hub", prob=prob0, seed=43219+r, vis = FALSE)
  LL1 = sugm.generator(n=n,d=d,graph="hub", prob=prob0, seed=14578+r, vis = FALSE)
  LLL1 = sugm.generator(n=n,d=d,graph="hub", prob=prob0, seed=96879+r, vis = FALSE)
 }
 if(network=="user_defined"){
        mu <- rep(0,d)
       true_theta <- prec
       sigma1 <- solve(prec)
       B1 <- gamma2
       B2 <- gamma1
      }
  else{
  true_theta = as.matrix(L$omega*L$theta)
  diag(true_theta)=1     ##theta is the precision matrix
  sigma1 <- L$sigma
  mu <- rep(0,d)
  thetaL = as.matrix(LL$omega*LL$theta)
  thetaLL = as.matrix(LLL$omega*LLL$theta)
  lwt = thetaL*(1*lower.tri(thetaL, diag =FALSE))
  upt = thetaLL*(1*upper.tri(thetaLL, diag = FALSE))
  B1c=upt+lwt
  ua=rbinom(d,1,0.02)
  uu=runif(d,0,1)
  uau = ua*uu
  diag(B1c)= uau
   B11 = B1c
 for(i in 1:d){
   for(j in 1:d){
    if(B1c[i,j] != 0)
       { m = rbinom(1,1,0.4)
         if(m ==0) B11[i,j] = -B1c[i,j]
       }

  }}
  B1=B11
  thetaL1 = as.matrix(LL1$omega*LL1$theta)
  thetaLL1 = as.matrix(LLL1$omega*LLL1$theta)
  lwt = thetaL1*(1*lower.tri(thetaL1, diag =FALSE))
  upt = thetaLL1*(1*upper.tri(thetaLL1, diag = FALSE))
  B2c=upt+lwt
  ua=rbinom(d,1,0.02)
  uu=runif(d,0,1)
  uau = ua*uu
  diag(B2c)= uau
   B22 = B2c
 for(i in 1:d){
   for(j in 1:d){
    if(B2c[i,j] != 0)
       { m = rbinom(1,1,0.4)
         if(m ==0) B22[i,j] = -B2c[i,j]
       }

  }}
  B2=B22
 
 }
#B1and B2 are the autoregressive coefficient matrices
##data generation
 xtn <-array(NA,c(t,d,n))
 xtt <- array(NA,c(t,d,1))
 for(i in 1:n){
  x0 <- rmvnorm(1, mu, sigma1, method="svd")
  x1 <- rmvnorm(1, mu, sigma1, method="svd")
  for(j in 1:t){
    et <- rmvnorm(1, mu, sigma1, method="svd")
    xt <-   x1 %*% B2 +  x0 %*% B1 + et
    xtt[j,,] <- xt
    x0 <- x1
    x1 <- xt
  }
   xtn[,,i] <- round(xtt,3)
  }
  true_gamma <- rbind(B2, B1)
  xy=matrix(aperm(xtn, c(3,1,2)), ncol=d)
  data1 <- as.longitudinal(xy, repeats=n)

 return(list(data1=data1,theta=true_theta, gamma=true_gamma))
 }
}


 