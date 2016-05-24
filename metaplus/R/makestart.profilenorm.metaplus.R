makestart.profilenorm.metaplus <- function(yi,sei,mods=NULL,fixed=NULL) {
  
  isreg <- !is.null(mods)
  
  ll.profilenorm <- function(par,yi,sei,mods) {
    isreg <- !missing(mods)
    muhat <- par[1]
    tau2 <- par[2]
    if (isreg) xcoef <- par[3:length(par)]
    w <- 1.0/(tau2+sei^2)
    if (isreg) negll <- 0.5*sum(log(2*pi)+log(1/w)+w*(yi-muhat-as.vector(mods %*% xcoef))^2)
    else negll <- 0.5*sum(log(2*pi)+log(1/w)+w*(yi-muhat)^2)
    if (is.nan(negll)) negll <- NA
    if (!is.finite(negll)) negll <- NA
    return(negll)
  }

  infixed <- fixed
  fixed <- unlist(infixed)
  names(fixed) <- names(infixed)
  
  # obtain starting values
  
  if (isreg) {
    start.meta <- rma(yi=yi, sei=sei, mods=as.data.frame(mods), method="DL")
    start.val <- c(start.meta$b[1,1],start.meta$tau2,start.meta$b[2:dim(start.meta$b)[1],1])
    lower.val <- c(-Inf,0.0,rep(-Inf,dim(mods)[2]))
  } else {
    start.meta <- rma(yi=yi, sei=sei, method="DL")
     start.val <- c(start.meta$b[1,1],start.meta$tau2)
    lower.val <- c(-Inf,0.0)
  }
  if (isreg) names(start.val) <- c("muhat","tau2",dimnames(mods)[[2]])
  else names(start.val) <- c("muhat","tau2")

  parnames(ll.profilenorm) <- names(start.val)
  
#  thefixed <- unlist(fixed)
#  if (length(thefixed)>0) {
#    fixedparms <- (1:length(start.val))[names(start.val)==names(thefixed)]
#    start.val <- start.val[-fixedparms]
#    lower.val <- lower.val[-fixedparms]
#  }
  
  names(lower.val) <- names(start.val)
  
  start.val <- start.val+0.001
  
  if (isreg) profilenorm.fit <- mymle(ll.profilenorm,start=start.val,vecpar=TRUE,
                                   optimizer="user",optimfun=myoptim,data=list(yi=yi,sei=sei,mods=as.matrix(mods)),
                                   skip.hessian=TRUE,
                                   lower=lower.val,
                                   fixed=fixed)
  else profilenorm.fit <- mymle(ll.profilenorm,start=start.val,vecpar=TRUE,
                            optimizer="user",optimfun=myoptim,data=list(yi=yi,sei=sei),
                             skip.hessian=TRUE,
                             lower=lower.val,
                             fixed=fixed)
  
  results <- c(as.list(profilenorm.fit@coef),fixed)
 if (!is.null(fixed)) {
   if (isreg) thenames <- c("muhat","tau2",dimnames(mods)[[2]])
   else thenames <- c("muhat","tau2")
   results <- unlist(results)
   results <- results[thenames]
   names(results) <- thenames
   results <- as.list(results)
 }
 return(params=results)
}
