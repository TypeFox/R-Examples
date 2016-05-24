#### (invisible function)
#### sub-function "get_mixcomp()" for select_poisson()
#### last version: 2015/03/26
#### This function is not meant to be called directly by the user. 

get_mixcomp <- function(y, mcomp){  
  
  i0  <- which(y==0)    # index of [y=0]
  n0  <- length(i0)     # length of [y=0]
  igz <- which(y > 0)   # index of [y>0]
  ygz <- y[igz]         # y>0
  ngz <- length(igz)
  
  if (all(ygz <= 30000)){
    my <- mcomp$m[ygz,]
    vy <- mcomp$v[ygz,]
    wy <- mcomp$w[ygz,]    
  } else {
    comp.m <- lapply(ygz, function(x){
      m <- compute.mixture(x, type = "log.gamma")$m
      return(as.data.frame(t(m)))
    })
    my <- as.matrix(do.call(rbind.fill, comp.m))
    dimnames(my) <- NULL
    my[is.na(my)] <- 0      
    if (ncol(my) < 10){
      my <- cbind(my, matrix(0, length(y), 10 - ncol(my)))
    }
    
    comp.v <- lapply(ygz, function(x){
      v <- compute.mixture(x, type = "log.gamma")$v
      return(as.data.frame(t(v)))
    })
    vy <- as.matrix(do.call(rbind.fill, comp.v))
    dimnames(vy) <- NULL
    vy[is.na(vy)] <- 0      
    if (ncol(vy) < 10){
      vy <- cbind(vy, matrix(0, length(y), 10 - ncol(vy)))
    }  
    
    comp.w <- lapply(ygz, function(x){
      w <- compute.mixture(x, type = "log.gamma")$p
      return(as.data.frame(t(w)))
    })
    wy <- as.matrix(do.call(rbind.fill, comp.w))
    dimnames(wy) <- NULL
    wy[is.na(wy)] <- 0    
    if (ncol(wy) < 10){
      wy <- cbind(wy, matrix(0, length(y), 10 - ncol(wy)))
    }     
  }
  return(list(i0 = i0, n0 = n0, igz = igz, ygz = ygz, ngz = ngz, 
              my = my, vy = vy, wy = wy))
}



#### (invisible function)
#### sub-function "iams1_poisson()" for select_poisson()
#### last version: 2015/03/26
#### This function is not meant to be called directly by the user. 

## Data augmentation step 1 of the "Improved Auxiliary Mixture Sampling" (IAMS)
## algorithm of Fruehwirth-Schnatter et al. (2009):
## Sample the (latent) inter-arrival times tau of an assumed Poisson process
## in the unit time interval with intensity [offset*lambda] 
## (see Fruehwirth-Schnatter et al., 2009)

iams1_poisson <- function(y, mu, compmix.pois){ 	
	tau         <- matrix(0, length(y), 2)
	tau[, 1]    <- rexp(length(y), mu)
	ta2         <- rbeta(length(y) - compmix.pois$n0, compmix.pois$ygz, 1)
	tau[compmix.pois$igz, 1] <- 1 - ta2 + tau[compmix.pois$igz, 1]
	tau[compmix.pois$igz, 2] <- ta2
	tau[compmix.pois$i0, 1]  <- 1 + tau[compmix.pois$i0, 1]
	t1         <- tau[, 1]
	if (compmix.pois$n0==0){
    t2 <- tau[, 2]
	} else {
    t2 <- tau[-compmix.pois$i0, 2]
	}
	return(list(t1 = t1, t2 = t2))
}



#### (invisible function)
#### sub-function "iams2_poisson()" for select_poisson()
#### last version: 2015/03/26
#### This function is not meant to be called directly by the user. 

## Data augmentation step 2 of the "Improved Auxiliary Mixture Sampling" (IAMS)
## algorithm of Fruehwirth-Schnatter et al. (2009):
## Approximation of the error distributions by Gaussian mixtures, where the
## latent component indicators r are introduced as missing data.

iams2_poisson <- function(n, tau1, tau2, logMu, logMugz, cm1, compmix){
  
  lwy <- log(compmix$wy)
  lwy[which(is.finite(lwy)==FALSE)] <- 0
  lvy <- log(compmix$vy)
  lvy[which(is.finite(lvy)==FALSE)] <- 0
  vy2 <- compmix$vy
  vy2[which(vy2==0)] <- 1
  kill <- matrix(compmix$vy > 0, compmix$ngz, 10)
  c2 <- (lwy - 0.5*lvy) 
  
  rgm   <- cm1$c1 - 0.5*((-log(tau1) - logMu)%*%t(rep(1,10)) - rep(1,n)%*%t(cm1$comp$m))^2/(rep(1,n)%*%t(cm1$comp$v))
  rgm[which(rgm==0)] <- NA
  mx <- apply(rgm, 1, max, na.rm = TRUE)
  e1    <- exp(rgm - mx)
  e1[which(is.na(e1))] <- 0
  rgmod <- e1/rowSums(e1, na.rm = TRUE)
  Fn    <- rgmod%*%upper.tri(matrix(1, 10, 10), diag = TRUE)
  
  # determination of random indicators R1
  u <- runif(n, 0, 1)
  R <- 11 - rowSums(u < Fn)
  
  xx     <- (-log(tau2) - logMugz)%*%t(rep(1, 10))*kill
  rgmx   <- c2 - 0.5*(xx - compmix$my)^2/vy2
  rgmx[which(rgmx==0)] <- NA
  mx <- apply(rgmx, 1, max, na.rm = TRUE)
  e1     <- exp(rgmx - mx)
  e1[which(is.na(e1))] <- 0 
  rgmodx <- e1/rowSums(e1, na.rm = TRUE)
  Fx     <- rgmodx%*%upper.tri(matrix(1, 10, 10), diag = TRUE)
  
  # determination of random indicators R2
  ux <- runif(compmix$ngz, 0, 1)
  R <- c(R, 11 - rowSums(ux < Fx))
  return(R)
}