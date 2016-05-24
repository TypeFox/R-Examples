#' Extract and Visualize GWRM Model Residuals
#'
#' residuals is a method which extracts model residuals from \code{"gw"}, commonly returned by \code{gw} function. Optionally, it produces a normal plot with a simulated envelope of the residuals.
#'
#' @param object	object of class \code{"gw"} holding the fitted model
#' @param type type of residuals to be extracted. Default is \code{"pearson"}. \code{"response"} and \code{"deviance"} are also available. Deviance residuals are defined as \eqn{2[ln f(y_i|y_i)-ln f(\widehat{\mu}_i|y_i)]}, so that their sum is the value of the deviance statistic.
#' @param rep	number of replications for envelope construction. Default is 19, that is the smallest 95 percent band that can be built.
#' @param envelope	a logical value to specify if the envelope is required.
#' @param title	a title for the envelope.
#' @param trace	if \code{TRUE} a sort of information is printed during the running time.
#' @param parallel if \code{TRUE} use parallel executation.
#' @param ncores is the number of cores that we use if \code{parallel} is \code{TRUE}.
#' @param ... 	further arguments passed to or from other methods.
#'
#' @details The usual Q-Q plot may show an unsatisfactory pattern of the residuals of a model fitted: then we are led to think that the model is badly specificated. The normal plot with simulated envelope indicates that under the distribution of the response variable the model is OK if only a few points fall off the envelope.
#' @return Residuals values and plot
#'
#' @examples
#' data(goals)
#' set.seed(1)
#' fit0 <- gw(goals ~ position, data = goals[sample(1:nrow(goals), 75), ])
#' residuals(fit0, type = "pearson", rep = 19, envelope = TRUE, trace = FALSE, ncores = 2)
#'
#' @importFrom stats qnorm
#' @importFrom graphics lines plot polygon
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster clusterExport detectCores
#' @importFrom foreach foreach %dopar% %do%
#' @export

residuals.gw <- function(object, type = "pearson", rep = 19, envelope = FALSE, title = "Simulated Envelope of Residuals", trace = FALSE, parallel=TRUE, ncores=2,  ...){

  resid.gw <- function(object, type){
    mu <- object$fitted.values
    k <- object$betaIIpars[1]
    ro <- object$betaIIpars[2]
    if (type == 'pearson'){
      if (ro < 2) stop('Variance is infinite')
      else{
        varianza <- ((k + ro - 1) / (ro - 2)) * (mu + (mu ^ 2) / k)
        residuos <- sort(object$residuals / sqrt(varianza))
        return(residuos)
      }
    }
    if (type == 'deviance'){
      diflogdgw<-function(p,k,ro){
        y<-p[1]
        a<-p[2]
        if(y==0) -lgamma(k+ro)-lgamma(a+ro)+lgamma(ro)+lgamma(a+k+ro)
        else{
          lgamma(y*(ro-1)/k+ro)+lgamma(y*(ro-1)/k+y)-lgamma(y*(ro-1)/k)-lgamma(y*(ro-1)/k+k+ro+y)-lgamma(a+ro)-lgamma(a+y)+lgamma(a)+lgamma(a+k+ro+y)
        }
      }
      y <- object$Y
      a <- mu*(ro-1)/k
      p <- cbind(y,a)
      residuos <- 2*sort(apply(p,1,k,ro,FUN=diflogdgw))
      return(residuos)
    }
    if (type == 'response'){
      residuos <- sort(object$residuals)
      return(residuos)
    }
  }

  if (envelope == FALSE) return(resid.gw(object, type))

  else{
    st <- proc.time()

    k <- object$betaIIpars[1]
    ro <- object$betaIIpars[2]
    mu <- object$fitted.values
    a <- mu * (ro - 1) / k

    n<-sum(object$W)

    if (parallel){
      cl <- makeCluster(ncores)
      registerDoParallel(cl)
      clusterExport(cl, list("object","rgw","type"),
                    envir=environment())

      #residuos.sim <- foreach(j=1:(rep-1), .combine = cbind,.multicombine = TRUE,.inorder=FALSE,.packages=c("GWRM"),.verbose=as.logical(trace)) %dopar%{
      residuos.sim <- foreach(j=1:rep, .combine = cbind,.multicombine = TRUE,.inorder=FALSE,.packages=c("GWRM"),.verbose=as.logical(trace)) %dopar%{
        converged<-FALSE
        varResponse<-getResponse(object$formula)
        datos <- object$data[rep(1:nrow(object$data), object$W),]
        datos[varResponse]<-as.matrix(rgw(n, a, k, ro))
        while(!converged){
          fit <- try(GWRM::gw(object$formula, data = datos, k = object$k), silent = TRUE)
          if(fit$aic>0 && fit$betaIIpars[2]>2)
            converged<-TRUE
          else{ ##Generate new response values
            datos[varResponse]<-as.matrix(rgw(n, a, k, ro))
          }
        }
        as.matrix(resid.gw(fit,type))
      }
      stopCluster(cl)
    }else{
      #residuos.sim <- foreach(j=1:(rep-1), .combine = cbind,.multicombine = TRUE,.inorder=FALSE,.packages=c("GWRM"),.verbose=as.logical(trace)) %do%{
      residuos.sim <- foreach(j=1:rep, .combine = cbind,.multicombine = TRUE,.inorder=FALSE,.packages=c("GWRM"),.verbose=as.logical(trace)) %do%{
        converged<-FALSE
        varResponse<-getResponse(object$formula)
        datos <- object$data[rep(1:nrow(object$data), object$W),]
        datos[varResponse]<-as.matrix(rgw(n, a, k, ro))
        while(!converged){
          fit <- try(GWRM::gw(object$formula, data = datos, k = object$k), silent = TRUE)
          if(fit$aic>0 && fit$betaIIpars[2]>2)
            converged<-TRUE
          else{ ##Generate new response values
            datos[varResponse]<-as.matrix(rgw(n, a, k, ro))
          }
        }
        as.matrix(resid.gw(fit,type))
      }
    }
    residuos <- resid.gw(object, type)
    residuos.sim<-cbind(as.matrix(rep(residuos,object$W)),residuos.sim)


    minimos <- apply(residuos.sim[, 2:rep], 1, min)
    maximos <- apply(residuos.sim[, 2:rep], 1, max)

    t <- 1:n
    normal.score <- qnorm(t / (n + 1))
    xx <- c(normal.score, rev(normal.score))
    yy <- c(minimos, rev(maximos))
    plot(normal.score, residuos, type = "l", xlab = "Standard normal quantiles", ylab = paste("Residuals ","(",type[1],")", sep = ""), main = title)
    polygon(xx, yy, col = "gray", border = NA)
    lines(normal.score, residuos)
#     for (j in 2:rep){
#       lines(normal.score, residuos.sim[,j])
#     }
  }

  et <- proc.time()
  if (trace > 0)
    cat(paste("\nOverall envelope simulation process took", round((et - st)[3], 2), 'seconds'))
  else cat("\n")

  ans <-  list(type = type, residuals = residuos, sim.residuals = residuos.sim)
  class(ans) <- "residuals.gw"
  return(ans)
}
