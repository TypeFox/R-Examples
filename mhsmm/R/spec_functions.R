mstep.pois <- function (x, wt) 
{
    k = ncol(wt)
    lambda = numeric(k)
    for (i in 1:k) lambda[i]=weighted.mean(x,wt[,i])
    list(lambda=lambda)
}

dpois.hsmm <- function (x, j, model) dpois(x,model$parms.emission$lambda[j])

rpois.hsmm <- function (j, model)  rpois(1, model$parms.emission$lambda[j])

rmvnorm.hsmm <- function(j,model) 
  rmvnorm(1,mean=model$parms.emission$mu[[j]],sigma=model$parms.emission$sigma[[j]])

mstep.mvnorm <- function(x, wt) {
  idx <-  apply(is.na(x),1,any) # Find rows with NA's (cov.wt does not like them)
  x  <- x[!idx,,drop=FALSE]
  wt <- wt[!idx,,drop=FALSE]
#  print(class(x))
  emission <- list(mu = list(), sigma = list())
   for (i in 1:ncol(wt)) {  ### CHANGE HERE: Must be wt, NOT x
   tmp <- cov.wt(x, wt[, i])
   emission$mu[[i]] <- tmp$center
   emission$sigma[[i]] <- tmp$cov
   }
   emission
} 

dmvnorm.hsmm <- function(x, j, model) {
  ans <- dmvnorm(x, mean = model$parms.emission$mu[[j]],
  sigma = model$parms.emission$sigma[[j]])
  ans[is.na(ans)] <- 1
  ans 
}

mstep.norm <- function(x,wt) {
    k = ncol(wt)
    mu = numeric(k)
    sigma = numeric(k)
    for(i in 1:k) {
      tmp = cov.wt(data.frame(x[!is.na(x)]),wt[!is.na(x),i])
      mu[i] = tmp$center
      sigma[i] = tmp$cov
    }
  list(mu=mu,sigma=sigma)
}

dnorm.hsmm <- function(x,j,model) {
  ret = dnorm(x,model$parms.emission$mu[j],sqrt(model$parms.emission$sigma[j]))
  ret[is.na(ret)] = 1
  ret           
}

rnorm.hsmm <- function(j,model)  rnorm(1,model$parms.emission$mu[j],sqrt(model$parms.emission$sigma[j]))
