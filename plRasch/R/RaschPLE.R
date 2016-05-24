### 'RaschPLE'  Fit Rasch Family Models using Pseudolikelihood
### Esitmation, capable of dealing with polytomous items, and
### multidimensional latent variables.
###
### Syntax:
###   RaschPLE(data, item.mtx, trait.mtx)
###
### Arguments:
###   data:
###   item.mtx:
###   trait.mtx:
###
### Value:
###
###
###
###
### Details:
###   The model is
###                       exp( w[i,h]' theta[v] + beta[i,h] )
###   P(X[v,i] = h) = ---------------------------------------------
###                    sum_l{ exp( w[i,l]' theta[v] + beta[i,l] ) }
###   where
###     X[v,i] is the response of vth individual to ith item;
###     w[i,h] is a vector of know category weight or score for
### response h of ith item;
###     theta[v] is a vector of latent traits for vth individual;
###     beta[i,h] is the item difficult parameter for ith item;
### associated with response h.
###
###   The function  only returns the item parameter beta.
###
###   Essentially, it is a wrapper function: the equvialent llla model
### are fitted.
###
### Reference:
###   Carolyn Anderson, Zhushan Li, and Jeroen Vermunt, 2006,
### Esitmation of Models in the Rasch Family for Polytomous Items and
### Multidimensional Latent Variables, Journal of Statistical
### Software.
###
### Author: Zhushan "Mandy" Li
### Last Updated: Aug 23, 2006

raschRap <- function(data, item.mtx, trait.mtx){
  RaschPLE(data, item.mtx, trait.mtx)
}

RaschPLE <- function(data, item.mtx, trait.mtx){
  ## fit the corresponding llla model
  lllafit <- llla(data, item.mtx, trait.mtx)
  llla2Rasch(lllafit,data)
}

### llla2Rasch converts the result from function 'llla',
### makes correction on item parameters, so that they are
### consistent with the Rasch model
llla2Rasch <- function(lllafit, data){
  stopifnot(inherits(lllafit, "llla"))

  ncat <- lllafit$ncat
  nitem <- lllafit$nitem
  nexaminee <- lllafit$nexaminee

  item.mtx <- lllafit$item.mtx
  trait.mtx <- lllafit$trait.mtx

  ##
  coef.estimates <-lllafit$coefficients;

  phi.names <- names(coef.estimates)[grep("phi", names(coef.estimates))]

  phi.coef.idx <- match(phi.names, names(lllafit$coefficients))
  phi.estimates <- coef.estimates[phi.coef.idx]
  nphi <- length(phi.estimates)

  item.coef.idx <- (1:length(coef.estimates))[-phi.coef.idx]
  item.estimates <- coef.estimates[item.coef.idx]

  ## The item parameter estimates given by llla need to be
  ## subtracted by a certain number to get the estimates
  ## of beta in Rasch Model

  pldata <- plStackData(data, item.mtx, trait.mtx)

  trait.low <- trait.mtx & lower.tri(trait.mtx,diag=TRUE)
  which.low <- which(trait.low, arr.ind=TRUE)

  ## contruct a matrix to keep track the paramters's original positions in the trait matrix
  trait.num <- matrix(0, nrow=NROW(trait.mtx), ncol=NCOL(trait.mtx))
  for(i in 1:nphi){
    trait.num[which.low[i,1], which.low[i,2]] <- i
    trait.num[which.low[i,2], which.low[i,1]] <- i
  }

  ##
  phi2trait.mtx <- matrix(0, nrow=nphi, ncol=NCOL(trait.mtx))
  for(i in 1:ncol(phi2trait.mtx) ) {
    phi2trait.mtx[trait.num[i,], i] <-   phi.estimates[trait.num[i,]]
  }
  rownames(phi2trait.mtx) <- phi.names
  colnames(phi2trait.mtx) <- colnames(trait.mtx)

  ##
  plphi <- pldata[,phi.names]
  plskills <- plphi %*% phi2trait.mtx

  plskills.mean <- apply(plskills,2,mean)
  offset <- sum(plskills.mean)

  ## return the estimates of beta
  beta <- item.estimates + rep(offset,ncat-1)*rep(1:(ncat-1), each=nitem)
  beta.se <- lllafit$se[item.coef.idx]
  beta.covb <- lllafit$covb[item.coef.idx, item.coef.idx]
  ##
  return(list(coefficients=beta, se=beta.se, covb=beta.covb))
}


## llla2Rasch.2 <- function(lllafit){
##   stopifnot(inherits(lllafit, "llla"))

##   ncat <- lllafit$ncat
##   nitem <- lllafit$nitem
##   nexaminee <- lllafit$nexaminee

##   ##
##   coef.estimates <-lllafit$coefficients;

##   phi.names <- names(coef.estimates)[grep("phi", names(coef.estimates))]

##   phi.coef.idx <- match(phi.names, names(lllafit$coefficients))
##   phi.estimates <- coef.estimates[phi.coef.idx]

##   item.coef.idx <- (1:length(coef.estimates))[-phi.coef.idx]
##   item.estimates <- coef.estimates[item.coef.idx]

##   ## The item parameter estimates given by llla need to be
##   ## subtracted by a certain number to get the estimates
##   ## of beta in Rasch Model

##   ## Construct the offset values
##   pl.offset <- plStackData(matrix((ncat-1)/2, nrow=1, ncol=nitem), item.mtx, trait.mtx)
##   pl.offset.phi <- as.matrix(pl.offset[, phi.names])
##   offset <- pl.offset.phi %*% phi.estimates

##   ## return the estimates of beta
##   beta <- item.estimates + rep(offset,ncat-1)*rep(1:(ncat-1), each=nitem)
##   beta.se <- lllafit$se[item.coef.idx]
##   beta.covb <- lllafit$covb[item.coef.idx, item.coef.idx]
##   ##
##   return(list(coefficients=beta, se=beta.se, covb=beta.covb))
## }


