#' @title Simulate Response Pattern under Dichotomous and Polytomous Rasch Model
#' @export simra
#' @description function for simulation of response patterns following the dichotomous and/or polytomous Rasch model based on the category probabilities given the model parameters.
#' At default, when just calling \code{simra()} 1 replication of responses to 5 items with difficulties -2, -1, 0, 1, 2 from 100 persons with ability drawn from N(0|1) are sampled.
#' @details no details in the moment.
#' @param itempar a "matrix" with \code{nrow = k} (number of items) and ncol = m (maximum number of thresholds), holding the 'thurstonian' thresholds of the respective item. Some of the rightmost matrix entries may be NA, depending on the number of categories of the respective item.  
#' @param theta either one of the following (1) a numeric vector of length \code{n} providing the values of the person parameter (ability) for \code{n} persons to be used in the simulation. (2) a integer defining the number of values \code{n} to draw from N(0|1).
#' @param pers_obj an object of class \code{"pers"} as a result from function \code{\link{pers}}. If an object of class \code{"pers"} is assigned to this argument the model parameters in it are taken to simulate the responses. At default (\code{pers_obj = NULL}) simulation is done by considering the model parameters given in the other arguments (see above).
#' @param replicate an integer defining how many replicates (data matrices) \code{r} to draw based on the model parameters.
#' @param seed a numeric vector with legnth of number of replications used for \code{\link{set.seed}} prior to each replicate to keep the result repeatable. If \code{seed = NULL} no seed is set.  
#' @param ... arguments passed through.
#' @return an array with \code{dim(n,k,r)} response patterns (\code{k} items in colums \code{n} persons in rows and \code{r} replications in the third dimension).
#' @examples ########
#' simra() # 100 dichotomous probabilistic response pattern
#' ### 100 polytomous response pattern (4 items; each 4 answer categories)
#' v <- c(-1.0,-0.5,0.0,0.5,-0.75,-0.25,0.25,0.75,-0.5,0.0,0.5,1.0)
#' itempar <- matrix(v,nrow = 4,ncol = 3)
#' simra(itempar = itempar)
#' simra(itempar = itempar,replicate = 10) # draw 10 replications
#' 
####################################################

simra <- function(itempar=matrix(seq(-2,2,length=5)), theta=100, pers_obj=NULL, replicate=1, seed=seq(1,replicate,1), ...  ){
  #-----------internal function ------------------ mit fehler 4-3-2015 --> korrigiert am 10-3-2015
  pvx.matrix<-function(theta_v,thres,xm_v=NULL){
    # func. by joerg-henrik heine jhheine(at)googlemail.com
    # ein dimension dazu und zum merken 
    # theta_v: ein vector oder zahl; 
    # thres: thurstonian thresholds eines items
    # xm_v: vector welche kategorie prob jeweils ausgegeben werden soll
    # korrigierte formel aus markus buch seite 330
    s<-0:length(thres)
    thres0<-c(0,thres)
    oben_v <- exp(apply((s%o%theta_v),2,function(x){x-cumsum(thres0)})) # ok - für theta_v als vector oder zahl
    unten_v<- apply( exp( apply((s%o%theta_v),2,function(x){x-cumsum(thres0)}) ) , 2 ,sum) # #sum in cumsum 
    px_v <- mapply(FUN=function(o,u){  o / u }, o=as.list(as.data.frame(oben_v)), u=  as.list(unten_v)     ) # u as list etc
    rownames(px_v)<-paste("cat",0:(length(thres0)-1),sep=".") 
    colnames(px_v)<-theta_v 
    P_v <- apply(px_v,2,sum) # test ok - für theta_v als vector oder zahl
    if(length(xm_v)==0){return( (px_v) )}
    if(length(xm_v)!=0){mapply(function(p,ic){p[ic]}, as.list(as.data.frame(px_v)), xm_v)}
  }     
  #--------------------------------------------------- 
  
  if (length(theta)==1){
    theta <- sort(rnorm(theta))
    Pname<-nchar(paste(length(theta)))
    personnames <- paste("P",formatC(1:length(theta), width = Pname, format = "d", flag = "0"),sep="")  
  }
  
  if (length(pers_obj)!=0){
    itemnames <- rownames(pers_obj$pair$threshold)
    personnames <- as.character(pers_obj$pers$persID)
    itempar <- pers_obj$pair$threshold
    theta <- pers_obj$pers$WLE
  }
  
  if(length(rownames(itempar))==0){
    Iname<-nchar(paste(dim(itempar)[1]))
    itemnames <- paste("I",formatC(1:nrow(itempar), width = Iname, format = "d", flag = "0"),sep="")  
  }else{itemnames <- rownames(itempar)}
  
  thresL <- lapply(1:nrow(itempar), function(i) {na.omit(itempar[i,])})
  names(thresL) <- rownames(itempar)
  
  probs <- lapply(thresL,FUN=function(x){t(pvx.matrix(theta_v=theta,thres=x))})
  names(probs) <- names(thresL)
  
  Rname<-nchar(paste(replicate))
  repnames <- paste("Rep",formatC(1:replicate, width = Rname, format = "d", flag = "0"),sep="")  
    
  k <- length(thresL) # number items
  n <- length(theta) # number persons
  m <- sapply(thresL, length)+1 # number of catgories per item 
  location <- sapply(thresL, mean) # not returned yet
  erg <- array(NA, dim=c(n,k,replicate))
  for (r in 1:replicate){
    if(length(seed)!=0){set.seed=seed[r]}
    X <- matrix(NA, n, k)
    for (j in 1:k) {
      for (i in 1:n) {
        if(length(seed)!=0){set.seed=seed[r]}
        X[i, j] <- sample(m[j], 1, prob = probs[[j]][i, ])
      }
    }
  X <- X-1
  #dimnames(X) <- list(personnames,itemnames)
  erg[ , ,r] <- X
  }
  dimnames(erg) <- list(personnames,itemnames,repnames)
  return(erg)
}
  


