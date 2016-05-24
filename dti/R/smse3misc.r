betagamma <- function(g){
  dg <- dim(g)
  ngrad <- if(!is.null(dg)) dg[2] else 1
  bghat <- .Fortran("bgstats",
                    as.double(g),
                    as.integer(ngrad),
                    double(2*ngrad),
                    bghat = double(2*ngrad*ngrad),
                    PACKAGE = "dti")$bghat
  dim(bghat) <- c(2, ngrad, ngrad)
  ## sphaerische Coordinaten fuer Gradienten-Paare
  bghat
}

matrm <- function(b, g){
  matrix(c(cos(b), 0, sin(b), 
           sin(b)*sin(g), cos(g), -cos(b)*sin(g), 
           -cos(g)*sin(b), sin(g), cos(b)*cos(g)),
         3, 3)
}

#
#
#

getkappas <- function(grad, trace = 0, dist = 1){
  #
  #  dist = 1: k4^2+k5^2+|k6|
  #  dist = 2: k4^2+k5^2+k6^2
  #  dist = 3: acos(g_i%*%g_j)
  # new version with analytic expm(par[1]*mx)
  #
  krit <- function(par, matm, beta){
    ## sum((matm-expm(par[1]*m4)%*%expm(par[2]*m5)%*%expm(par[3]*m6))^2)
    .Fortran("k456krb",
             as.double(par),
             as.double(beta),
             as.double(matm),
             erg = double(1),
             PACKAGE = "dti")$erg
  }
  krit5 <- function(x,p,pk4,matm,beta){
    # for line search with respect to k5 to get second solution
    p1 <- p+c(pk4/2,x,pi)
    krit(p1,matm,beta)
  }
  ngrad <- dim(grad)[2]
  if(dist<3){
    prta <- Sys.time()
    cat("Start computing spherical distances", format(Sys.time()), "\n")
    kappa456 <- kappa456a <- array(0, c(3, ngrad, ngrad))
    bghat <- betagamma(grad)
    for (i in 1:ngrad) for (j in 1:ngrad) {
      bg <- bghat[, i, j]
      #   fix for discontinuity
      if(abs(cos(bg[1])) < 1.e-6) bg[1] = pi/2 - 1e-6*sign(cos(bg[1]))
      matm <- matrm(bg[1], bg[2])
      cbg1 <- cos(bg[1])
      k456 <- runif(3, -.1, .1)
      maxit <- 10000
      fnscale <- 1
      z <- optim(k456, krit, method = "BFGS", matm = matm, beta=bg[1],
                 control = list(trace = trace, fnscale=fnscale, maxit=maxit, reltol = 1e-12, abstol = 1e-16))
      count <- 5
      while (z$value > 1e-14&count>0) {
        ## cat("i",i,"j",j,"value",z$value,"par",z$par,"\n")
        maxit <- maxit*1.5
        fnscale <- fnscale/2
        k456 <- runif(3, -1, 1)
        z <- optim(k456, krit, method = "BFGS", matm = matm, beta=bg[1],
                   control = list(trace = trace, fnscale=fnscale, maxit=maxit, reltol = 1e-12, abstol = 1e-16))
        ## cat(" new value",z$value,"par",z$par,"\n")
        count <- count - 1
        if(count==0)
          cat("failed for bg:",bg,"value",z$value,"par",z$par,"\n")
      }
      z$par[1] <- z$par[1]/cbg1
      kappa456[, i, j] <- z$par
      pk4 <- abs(2*pi*cbg1/(2-cbg1^2)^.5)
      if(kappa456[1,i,j] < -pk4/2) kappa456[1,i,j] <- kappa456[1,i,j] - trunc(kappa456[1, i, j]/pk4-1)*pk4
      if(kappa456[1,i,j] >= pk4/2) kappa456[1,i,j] <- kappa456[1,i,j] - trunc(kappa456[1, i, j]/pk4+1)*pk4
      if(kappa456[2,i,j] < -pi) kappa456[2,i,j] <- kappa456[2,i,j] - trunc(kappa456[2, i, j]/2/pi-1)*2*pi
      if(kappa456[2,i,j] > pi) kappa456[2,i,j] <- kappa456[2,i,j] - trunc(kappa456[2, i, j]/2/pi+1)*2*pi
      if(kappa456[3,i,j] < -pi) kappa456[3,i,j] <- kappa456[3,i,j] - trunc(kappa456[3, i, j]/2/pi-1)*2*pi
      if(kappa456[3,i,j] > pi) kappa456[3,i,j] <- kappa456[3,i,j] - trunc(kappa456[3, i, j]/2/pi+1)*2*pi
      kappa456a[,i,j] <- kappa456[,i,j]+c(pk4/2,0,pi)
      kpar <- kappa456[,i,j]*c(cbg1,1,1)
      kappa456a[2,i,j] <- optimize(krit5,c(-pi,pi),p=kpar,pk4=pk4,
                                   matm=matm,beta=bg[1],maximum=FALSE)$minimum
      if(kappa456a[1,i,j] < -pk4/2) kappa456a[1,i,j] <- kappa456a[1,i,j] - trunc(kappa456a[1, i, j]/pk4-1)*pk4
      if(kappa456a[1,i,j] >= pk4/2) kappa456a[1,i,j] <- kappa456a[1,i,j] - trunc(kappa456a[1, i, j]/pk4+1)*pk4
      if(kappa456a[2,i,j] < -pi) kappa456a[2,i,j] <- kappa456a[2,i,j] - trunc(kappa456a[2, i, j]/2/pi-1)*2*pi
      if(kappa456a[2,i,j] > pi) kappa456a[2,i,j] <- kappa456a[2,i,j] - trunc(kappa456a[2, i, j]/2/pi+1)*2*pi
      if(kappa456a[3,i,j] < -pi) kappa456a[3,i,j] <- kappa456a[3,i,j] - trunc(kappa456a[3, i, j]/2/pi-1)*2*pi
      if(kappa456a[3,i,j] > pi) kappa456a[3,i,j] <- kappa456a[3,i,j] - trunc(kappa456a[3, i, j]/2/pi+1)*2*pi
    }
    dka <- switch(dist,kappa456[1,,]^2+kappa456[2,,]^2+abs(kappa456[3,,]),
                  kappa456[1,,]^2+kappa456[2,,]^2+kappa456[3,,]^2)
    dkb <- switch(dist,kappa456a[1,,]^2+kappa456a[2,,]^2+abs(kappa456a[3,,]),
                  kappa456a[1,,]^2+kappa456a[2,,]^2+kappa456a[3,,]^2)
    dim(kappa456) <- dim(kappa456a) <- c(3,ngrad*ngrad)
    kappa456[,dkb<dka] <- kappa456a[,dkb<dka]
    dim(kappa456) <- c(3,ngrad,ngrad)
    prtb <- Sys.time()
    cat("End computing spherical distances", format(Sys.time()), "\n")
  } else {
    kappa456 <- array(0, c(3, ngrad, ngrad))
    bghat <- betagamma(grad)
    for(i in 1:ngrad) kappa456[1,i,] <- acos(pmin(1,abs(grad[,i]%*%grad)))
  }
  list(k456 = kappa456, bghat = bghat, dist=dist)
}




suggestkappa <- function(grad,vred=1,dist=1){
  #
  #  get a kappa value from variance reduction on the sphere
  #
  gstats <- getkappas(grad,dist=dist)
  ngrad <- dim(grad)[2]
  vred <- min(vred,ngrad-1)
  d <- switch(dist,apply(gstats$k456[1:2,,]^2,2:3,sum)+abs(gstats$k456[3,,]),
              apply(gstats$k456^2,2:3,sum),
              apply(gstats$k456^2,2:3,sum))
  kmin <- sqrt(min(d[d>1e-8]))# just to prevent from taking zero
  kappa <- kmin
  vredk <- 1
  while(vredk < vred){
    kappa <- kappa*1.005
    w <- matrix(pmax(1-d/kappa^2,0),ngrad,ngrad)
    vredk <- mean(apply(w,1,sum)^2/apply(w^2,1,sum))
  }
  list(kappa=kappa,vred=vredk)
}
vredsphere <- function(grad,kappa,dist=1){
  #
  #  compute initial variance reduction on the sphere 
  #  for given kappa
  #
  gstats <- getkappas(grad,dist=dist)
  ngrad <- dim(grad)[2]
  d <- switch(dist,apply(gstats$k456[1:2,,]^2,2:3,sum)+abs(gstats$k456[3,,]),
              apply(gstats$k456^2,2:3,sum),
              apply(gstats$k456^2,2:3,sum))
  w <- matrix(pmax(1-d/kappa^2,0),ngrad,ngrad)
  mean(apply(w,1,sum)^2/apply(w^2,1,sum))
}


