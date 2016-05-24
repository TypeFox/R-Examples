#' Prediction from fitted Item focussed Trees 
#' 
#' @description
#' The function returns predictions of item parameters 
#' obtained by item focussed recursive partitioning based on the Rasch Model or the 
#' Logistic Regression Approach for DIF detection.  
#' 
#' @param object Object of class \code{\link[DIFtree]{DIFtree}}
#' @param item Number of the item, for which the prediction shall be returned
#' @param newdata New data.frame, for which the prediction shall be returned
#' @param ... Further arguments passed to or from other methods 
#' 
#' @details
#' For \code{"Rasch"} model the function returns the predicted item difficulty. 
#' For \code{"Logistic"} models the function returns the predicted intercept and/or slope.  
#' 
#' @author Moritz Berger <moritz.berger@stat.uni-muenchen.de> \cr \url{http://www.statistik.lmu.de/~mberger/}
#' 
#' @references 
#' Berger, Moritz and Tutz, Gerhard (2015): Detection of Uniform and Non-Uniform Differential Item Functioning 
#' by Item Focussed Trees, Cornell University Library, arXiv:1511.07178
#' 
#' Tutz, Gerhard and Berger, Moritz (2015): Item Focused Trees for the Identification of Items
#' in Differential Item Functioning, Psychometrika, published online, DOI: 10.1007/s11336-015-9488-3 
#' 
#' @seealso \code{\link[DIFtree]{DIFtree}}, \code{\link[DIFtree]{plot.DIFtree}}, \code{\link[DIFtree]{summary.DIFtree}}
#' 
#' @examples 
#' data(data_sim)
#'  
#' Y <- data_sim[,1]
#' X <- data_sim[,-1]
#' 
#' Xnew <- data.frame("x1"=c(0,1),"x2"=c(-1.1,2.5),"x3"=c(1,0),"x4"=c(-0.2,0.7))
#'  
#' \dontrun{
#'  
#' mod <- DIFtree(Y=Y,X=X,model="Logistic",type="udif",alpha=0.05,nperm=1000,trace=TRUE)
#'  
#' predict(mod,item=1,Xnew)
#' }
#' 
#
#' @method predict DIFtree
#' @export

predict.DIFtree <-
function(object, # object of class DIFtree
                            item,
                            newdata,...){
  
  # check input 
  if(!is.data.frame(newdata)){
    stop("X must be of class 'data.frame'")
  }

  X <- prepareX(newdata)

  if(ncol(X)!=ncol(object$X)){
    stop("X does not have the right structure!")
  }
  if(!all(names(X) %in% names(object$X))){
    stop("X does not have the right structure!")
  }
  n_pred <- nrow(X)
  
  model <- which(c("Rasch","Logistic") %in% paste(object$call))
  if(model==1){
    if(is.null(object$splits)){
      params_hat <- rep(object$coefficients$betas_nodif[paste0("beta",item)],n_pred)
    } else{
      info <- object$splits[object$splits[,"item"]==item,]
      if(nrow(info)==0){
       params_hat <- rep(object$coefficients$betas_nodif[paste0("beta",item)],n_pred)
      } else{
        params <- object$coefficients$betas_dif[[paste(item)]]
        params_hat <- whole_prediction(info,item,params,X)
      }
    }
  }
  if(model==2){
    type <- which(c("udif","dif","nudif") %in% paste(object$call))
    if(type==1){
      if(is.null(object$splits)){
        params_hat <- rep(object$coefficients$gammas_nodif[paste0("gamma",item)],n_pred)
      } else{
        info <- object$splits[object$splits[,"item"]==item,]
        if(nrow(info)==0){
          params_hat <- rep(object$coefficients$gammas_nodif[paste0("gamma",item)],n_pred)
        } else{
          params <- object$coefficients$gammas_dif[[paste(item)]]
          params_hat <- whole_prediction(info,item,params,X)
        }
      }
    }
    if(type==2){
      dif_items  <- unique(c(object$splits[[1]][,"item"],object$splits[[2]][,"item"]))
      params_hat <- matrix(NA,nrow=n_pred,ncol=2)
      colnames(params_hat) <- c("intercept","slope")
      if(!item %in% dif_items){
        params_hat[,1] <- rep(object$coefficients$gammas_nodif[paste0("gamma",item)],n_pred)
        params_hat[,2] <- rep(object$coefficients$alphas_nodif[paste0("alpha",item)],n_pred)
      } else{
        info1 <- object$splits[[1]][object$splits[[1]][,"item"]==item,]
        if(is.null(info1) || nrow(info1)==0){
          params_hat[,1] <- rep(object$coefficients$gammas_dif[[paste(item)]],n_pred)
        } else{
          params <- object$coefficients$gammas_dif[[paste(item)]]
          params_hat[,1] <- whole_prediction(info1,item,params,X)
        }
        info2 <- object$splits[[2]][object$splits[[2]][,"item"]==item,]
        if(is.null(info2) || nrow(info2)==0){
          params_hat[,2] <- rep(object$coefficients$alphas_dif[[paste(item)]],n_pred)
        } else{
          params <- object$coefficients$alphas_dif[[paste(item)]]
          params_hat[,2] <- whole_prediction(info2,item,params,X)
        }
      }
    }
    if(type==3){
      params_hat <- matrix(NA,nrow=n_pred,ncol=2)
      colnames(params_hat) <- c("intercept","slope")
      if(is.null(object$splits)){
        params_hat[,1] <- rep(object$coefficients$gammas_nodif[paste0("gamma",item)],n_pred)
        params_hat[,2] <- rep(object$coefficients$alphas_nodif[paste0("alpha",item)],n_pred)
      } else{
        info <- object$splits[object$splits[,"item"]==item,]
        if(nrow(info)==0){
          params_hat[,1] <- rep(object$coefficients$gammas_nodif[paste0("gamma",item)],n_pred)
          params_hat[,2] <- rep(object$coefficients$alphas_nodif[paste0("alpha",item)],n_pred)
        } else{
          params1 <- object$coefficients$gammas_dif[[paste(item)]]
          params_hat[,1] <- whole_prediction(info,item,params1,X)
          params2 <- object$coefficients$alphas_dif[[paste(item)]]
          params_hat[,2] <- whole_prediction(info,item,params2,X)
        }
      }
    }
  }
  return(params_hat)
  
  invisible(object)
}
