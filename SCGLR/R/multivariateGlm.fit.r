#' @rdname multivariateGlm
#' @export 
#' @importFrom stats as.formula glm 
#' @param Y matrix of dependent variables.
#' @param comp matrix of covariates.
# @param family a vector of character giving the family distribution of each response
# @param offset used for the poisson dependent variables.
# A vector or a matrix of size: number of observations * number of Poisson dependent variables is expected
# @param size a matrix giving the number of trials for each Binomial dependent variable
# ncol(size) must be equal to the number of Binomial variables
# @return the list, each item of which is the glm object associated with each response
multivariateGlm.fit <- function(Y,comp,family,offset,size){
  q <- ncol(Y)
  out <- list()
  if(!is.null(comp)){
    if (is.matrix(comp)|| is.data.frame(comp))
    {
      if(is.null(colnames(comp))){
      colnames(comp) <- paste("c",1:ncol(comp),sep="")
      }
      form <- as.formula(paste("obs~",paste(colnames(comp),collapse="+")))
      ##out <- matrix(0,ncol(comp)+1,ncol(Y))
    }else
    {
      form <- as.formula("obs~comp")
      ##out <- matrix(0,2,ncol(Y))
    }
  }else{
    form <- as.formula("obs~1")
    comp <- rep(1,nrow(Y))
  }
  kpois <- kbern <-  1
  for(i in seq(q)){
    if(family[i]=="poisson"){
      if(is.null(offset)){
        out[[i]] <- glm(form,data=data.frame(obs=Y[,i],comp),family="poisson")
      }else{
        out[[i]] <- glm(form,data=data.frame(obs=Y[,i],comp),family=family[i],
                        offset=log(offset[,kpois]))
        kpois <- kpois+1
      }
    }
    if(family[i]=="bernoulli"){
      out[[i]] <- glm(form,data=data.frame(obs=Y[,i],comp),family="binomial")  
    }
    if(family[i]=="binomial"){
      out[[i]] <- glm(form,data=data.frame(obs=Y[,i],comp),family="binomial",
                      weights=size[,kbern])
      kbern <- kbern+1
    }
    if(family[i]=="gaussian"){
      out[[i]] <- glm(form,data=data.frame(obs=Y[,i],comp),family="gaussian")
    }  
  }
  
  return(out)
}

