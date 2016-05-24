# October 14 2009 - Rewriting this....

#library(Hmisc) 
#library(MASS)





distmat <- function(vec){ as.matrix(dist(vec)) }
 
semisupercritmultidim <- function(d,z,alpha,y){
  tots <- 0
  if(is.matrix(z)){
    for(i in 1:ncol(z)){
      tots <- tots+sum((d/sqrt(ncol(z)) - outer(as.numeric(z[,i]), as.numeric(z[,i]), "-"))[y==2,y==1]^2)
     }
  } else {
     tots <- sum((d-outer(as.numeric(z), as.numeric(z), "-"))[y==2, y==1]^2)
  }
  dz <- distmat(z)
  stress <- 0.5*sum((d-dz)^2)
  stress.between <- sum(((d - dz)[y==1,y==2])^2)
  stress.within <- 0.5*(sum(((d - dz)[y==1,y==1])^2) + sum(((d - dz)[y==2,y==2])^2))
  return(list(crits=alpha*tots+(1-alpha)*stress, stress=stress, super=tots, stress.between=stress.between,stress.within=stress.within))
}  


TrainSuperMDSOnce <- function(d=NULL, y, alpha=.5, S=2, x=NULL, z=NULL){
  if(alpha<0 || alpha>1) stop("alpha must be between zero and 1.")
  if(!is.null(x) && !is.null(d)) stop("Please enter either x or d; the other should be NULL.")
  if(!is.null(x)) d <- distmat(x)
  if(is.null(x) && is.null(d)) stop("Either x or d must be non-null.")
  if(is.null(z))  z <- cmdscale(d, k=S)
  n1 <- sum(y==1); n2 <- sum(y==2); n <- length(y)
  if(n1+n2 != n) stop("Each element of y must be a 1 or a 2.")
  stress.between <- stress.within <- stress <- super <- crits <- NULL
  crit.out <- semisupercritmultidim(d,z,alpha,y)
  crits <- c(crits, crit.out$crits)
  stress <- c(stress, crit.out$stress)
  super <- c(super, crit.out$super)
  stress.between <- c(stress.between, crit.out$stress.between)
  stress.within <- c(stress.within, crit.out$stress.within)
  while(length(crits)<2 || abs(crits[length(crits)]-crits[length(crits)-n])/min(crits) > 1e-2){
    for(i in 1:n){
      if(S>1){
        if(y[i]==1){
          znew <- (1-alpha)*apply(z[-i,],2,sum) + alpha*apply(z[y==2,], 2, sum) - (alpha/sqrt(S))*sum(d[i,y==2]) - (1-alpha)*apply(sweep(scale(z[-i,], center=z[i,], scale=FALSE), 1, d[i,-i]/distmat(z)[i,-i], "*"), 2, sum)
          znew <- znew/((n-1)*(1-alpha)+n2*alpha)
        } else if (y[i]==2){
          znew <- (1-alpha)*apply(z[-i,],2,sum) + alpha*apply(z[y==1,], 2, sum) + (alpha/sqrt(S))*sum(d[i,y==1]) - (1-alpha)*apply(sweep(scale(z[-i,], center=z[i,], scale=FALSE), 1, d[i,-i]/distmat(z)[i,-i], "*"), 2, sum)
          znew <- znew/((n-1)*(1-alpha)+n1*alpha)
        }
        z[i,] <- znew
        crit.out <- semisupercritmultidim(d,z,alpha,y)
        crits <- c(crits, crit.out$crits)
        stress <- c(stress, crit.out$stress)
        super <- c(super, crit.out$super)
        stress.between <- c(stress.between, crit.out$stress.between)
        stress.within <- c(stress.within, crit.out$stress.within)
      } else {
        if(y[i]==1){
          znew <- (1-alpha)*sum(z[-i]) + alpha*sum(z[y==2]) - alpha*sum(d[i,y==2]) - (1-alpha)*sum(sign(z[-i]-z[i])*(d[i,-i]))
          znew <- znew/((n-1)*(1-alpha)+n2*alpha)
        } else if (y[i]==2){
          znew <- (1-alpha)*sum(z[-i,]) + alpha*sum(z[y==1]) + alpha*sum(d[i,y==1]) - (1-alpha)*sum(sign(z[-i]-z[i])*(d[i,-i]))
          znew <- znew/((n-1)*(1-alpha)+n1*alpha)
        }
        z[i] <- znew
        crit.out <- semisupercritmultidim(d,z,alpha,y)
        crits <- c(crits, crit.out$crits)
        stress <- c(stress, crit.out$stress)
        super <- c(super, crit.out$super)
        stress.within <- c(stress.within, crit.out$stress.within)
        stress.between <- c(stress.between, crit.out$stress.between)
      }
    }
    if(sum(is.na(crits)>0)) stop("some crits are NA")
  }
  teout <- TestSuperMDS(list(crits=crits,z=z, d=d, x=x, y=y, alpha=alpha), dtetr=d)
  cutpoint <- quantile(teout$crit1s-teout$crit2s, sum(y==1)/length(y))
  return(list(crits=crits,stress=stress[length(stress)], super=super[length(super)], z=z, d=d, x=x, y=y, alpha=alpha, cutpoint=cutpoint,
              stress.between=stress.between[length(stress.between)], stress.within=stress.within[length(stress.within)]))
}

evaltestcrit <- function(dnew, S, z, znew, ynew, alpha, y){
  n <- nrow(z)
  if(S>1) tots <- (1-alpha)*sum((dnew-distmat(rbind(z,znew))[n+1, -(n+1)])^2)
  if(S==1) tots <- (1-alpha)*sum((dnew-abs(z-znew))^2)
  if(ynew==1){
    if(S>1){
      for(s in 1:S){
        tots <- tots + alpha*sum((dnew[y==2]/sqrt(S) - (z[y==2,s]-znew[s]))^2)
      }
    } else if (S==1){
      tots <- tots+alpha*sum((dnew[y==2]-(z[y==2]-znew))^2)
    }
  } else if (ynew==2){
    if(S>1){
      for(s in 1:S){
        tots <- tots + alpha*sum((dnew[y==1]/sqrt(S) - (znew[s]-z[y==1,s]))^2)
      }
    } else if (S==1){
      tots <- tots+alpha*sum((dnew[y==1]-(znew-z[y==1]))^2)
    }
  }
  return(tots)
}


TestSuperMDSSingleObs <- function(dnew, alpha, z, y, S){
  n <- length(y); n1 <- sum(y==1); n2 <- sum(y==2)
  # Case 1: ynew=1
  if(S>1) znew <- apply(z[y==1,], 2, sum)
  if(S==1) znew <- mean(z[y==1])
  crit1 <- c(evaltestcrit(dnew, S, z, znew, 1, alpha,y))
  while(length(crit1)<3 || abs(crit1[length(crit1)]-crit1[length(crit1)-1])/crit1[length(crit1)] > 1e-4){
    if(S>1) znew <- alpha*apply(z[y==2,],2,sum) - (alpha/sqrt(S))*sum(dnew[y==2]) + (1-alpha)*apply(z, 2, sum) - (1-alpha)*apply(sweep(scale(z, center=znew, scale=FALSE), 1, dnew/distmat(rbind(z, znew))[n+1,-(n+1)], "*"), 2, sum, na.rm=TRUE)
    if(S==1) znew <- alpha*sum(z[y==2]) - alpha*sum(dnew[y==2]) + (1-alpha)*sum(z) - (1-alpha)*sum(dnew*sign(z-znew))
    znew <- znew/(n2*alpha+n*(1-alpha))
    crit1 <- c(crit1, evaltestcrit(dnew, S, z, znew, 1, alpha,y))
  }
  if(max(diff(crit1))>1e-8) print("Max(diff(crit1)) too big.")
  znew1 <- znew
  # Case 2: ynew=2
  if(S>1) znew <- apply(z[y==2,], 2, sum)
  if(S==1) znew <- mean(z[y==2])
  crit2 <- c(evaltestcrit(dnew, S, z, znew, 2, alpha, y))
  while(length(crit2)<3 || abs(crit2[length(crit2)]-crit2[length(crit2)-1])/crit2[length(crit2)] > 1e-4){
    if(S>1) znew <- alpha*apply(z[y==1,],2,sum) + (alpha/sqrt(S))*sum(dnew[y==1]) + (1-alpha)*apply(z, 2, sum) - (1-alpha)*apply(sweep(scale(z, center=znew, scale=FALSE), 1, dnew/distmat(rbind(z, znew))[n+1,-(n+1)], "*"), 2, sum, na.rm=TRUE)
    if(S==1) znew <- alpha*sum(z[y==1]) + alpha*sum(dnew[y==1]) + (1-alpha)*sum(z) - (1-alpha)*sum(dnew*sign(z-znew))
    znew <- znew/(n1*alpha+n*(1-alpha))
    crit2 <- c(crit2, evaltestcrit(dnew, S, z, znew, 2, alpha, y))
  }
  if(max(diff(crit2))>1e-8) print("Max(diff(crit2)) too big.")
  znew2 <- znew
  return(list(znew1=znew1, crit1=min(crit1), znew2=znew2, crit2=min(crit2)))
}


