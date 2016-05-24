#' @title Function that predicts the responses from the covariates for a new sample
#' @export
#' @importFrom stats as.formula model.matrix
#' @description  Function that predicts the responses from the covariates for a new sample
#' @param Xnew a data frame containing the values of the covariates for the new sample.
#' @param family a vector of character specifying the distributions of the responses.
#' @param beta the matrix of coefficients estimated from the calibration sample.
#' @param offset used for the poisson dependent variables.
#' A vector or a matrix of size: number of observations * number of Poisson dependent variables is expected.
#' @return a matrix of predicted values.
#' @examples \dontrun{
#' library(SCGLR)
#' 
#' # load sample data
#' data(genus)
#' 
#' # get variable names from dataset
#' n <- names(genus)
#' ny <- n[grep("^gen",n)]    # Y <- names that begins with "gen"
#' nx <- n[-grep("^gen",n)]   # X <- remaining names
#' 
#' # remove "geology" and "surface" from nx as surface is
#' # offset and we want to use geology as additional covariate
#' nx <-nx[!nx%in%c("geology","surface")]
#' 
#' # build multivariate formula
#' # we also add "lat*lon" as computed covariate
#' form <- multivariateFormula(ny,c(nx,"I(lat*lon)"),c("geology"))
#' 
#' # split genus dataset
#' sub <- sample(1:nrow(genus),100,replace=FALSE)
#' sub_fit <- (1:nrow(genus))[-sub]
#' 
#' # define family 
#' fam <- rep("poisson",length(ny))
#' 
#' # fit the model
#' genus.scglr <- scglr(formula=form, data=genus, family=fam, K=4, 
#'  offset=genus$surface, subset=sub_fit)
#' 
#' # xnew, the design matrix associated to sub-sample used for prediction
#' # rhs parameters is introduced to take into account that the covariate 
#' # of the formula is composed of two differents sets
#' xnew <- model.matrix(form, data=genus[sub,], rhs=1:2)[,-1]
#' 
#' # prediction based on the scglr approch
#' pred.scglr <- multivariatePredictGlm(xnew,family=fam,
#'  beta=genus.scglr$beta,offset=genus$surface[sub])
#' cor.scglr <-diag(cor(pred.scglr,genus[sub,ny])) 
#' plot(cor.scglr, col="red",ylim=c(-1,1))
#' 
#' # prediction based on classical poisson glm
#' genus.glm <- multivariateGlm(formula=form, data=genus, family=fam, 
#'  offset=genus$surface, subset=sub_fit)
#' coefs <- sapply(genus.glm,coef)
#' 
#' # rhs parameters is introduced to take into account that the covariate
#' # part of the formula is composed of two differents sets
#' pred.glm <- multivariatePredictGlm(xnew,family=fam,beta=coefs,
#'  offset=genus$surface[sub])
#' cor.glm <- diag(cor(pred.glm,genus[sub,ny]))
#' 
#' points(cor.glm, col="blue")
#' }
multivariatePredictGlm <- function(Xnew,family,beta,offset=NULL){
  if(is.matrix(Xnew)|| is.data.frame(Xnew)){
  Xnew <- as.data.frame(Xnew)
  if(is.null(names(Xnew))) names(Xnew) <- paste("x",1:ncol(Xnew),sep="")
  form <- as.formula(paste("~",paste(names(Xnew),collapse="+"),sep="")) 
  x <- model.matrix(form, data=Xnew)
  prediction <- x%*%as.matrix(beta)
  }else{
    prediction <- Xnew%*%matrix(beta,nrow=1)
  }
  if(!is.null(offset)) {
    if(is.vector(offset)) {
      offset <- matrix(offset,nrow(Xnew), sum(family=="poisson"))
    } else {
      offset <- as.matrix(offset)
    }
  } 
  if(sum(family=="bernoulli")>0){
    tmp <- prediction[,which(family=="bernoulli")]
    prediction[,which(family=="bernoulli")] <- exp(tmp)/(1+exp(tmp))
  }
  if(sum(family=="binomial")>0){
    tmp <- prediction[,which(family=="binomial")]
    prediction[,which(family=="binomial")] <- exp(tmp)/(1+exp(tmp))
  }
  if(sum("poisson"%in%family)>0){
    if(is.null(offset)){
      tmp <-  prediction[,which(family=="poisson")]  
    }else{
      tmp <- log(offset)+prediction[,which(family=="poisson")] 
    }
      prediction[,which(family=="poisson")] <- exp(tmp)
  }
  return(prediction)
}
