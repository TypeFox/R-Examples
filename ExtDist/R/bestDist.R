#' @title Finding the best distribution for a (weighted) sample.
#' @description This function chooses the best fitted distribution, based on 
#' a specified criterion.
#' @rdname bestDist
#' @name bestDist
#' @details When comparing models fitted by maximum likelihood to the same data, the smaller the AIC, BIC or MDL, the better the fit.
#' When comparing models using the log-likelihood criterion, the larger the log-likelihood the better the fit.
#' @note The MDL criterion only works for parameter estimation by numerical maximum likelihood.
#' @param X Sample observations.
#' @param w An optional vector of sample weights.
#' @param candDist A vector of candidate distributions.
#' @param criterion The basis on which the best fitted distribution is chosen.
#' 
#' @author Haizhen Wu and A. Jonathan R. Godfrey.
#' 
#' @return An object of class character containing the name of the best distribution and its corresponding parameter estimates.
#' @export bestDist
#'
#' @examples
#' X <- rBeta_ab(30, a = 0, b = 1, shape1 = 2, shape2 = 10)
#' 
#' # Determining the best distribution from the list of candidate distributions for the data X
#' Best.Dist <- bestDist(X, candDist = c("Laplace","Normal","Beta_ab"), criterion = "logLik")
#' 
#' # Printing the parameter estimates of the best distribution
#' attributes(Best.Dist)$best.dist.par

bestDist <- function(X, w = rep(1,length(X))/length(X), candDist = c("Beta_ab","Laplace","Normal"), criterion = c("AICc", "logLik",
                      "AIC", "BIC", "MDL")){
  
    if(!(criterion %in% c("logLik","AIC","AICc","BIC","MDL"))) {return("criterion unknown")}
    
    w <- w/sum(w)*length(X)
    
    est.pars <- lapply(paste0("e",candDist), do.call, args=list(X,w)) 
    
    criterion.value <- unlist(lapply(est.pars, getS3method(criterion, class = "eDist")))
    names(criterion.value) <- candDist
    
    if(criterion %in% c("logLik")){
      # When comparing models fitted by maximum likelihood to the same data,
      # the larger the log-likelihood the better the fit. 
      new.order <- order(criterion.value,decreasing=T)
      best <- candDist[new.order][1]
      est.par <- est.pars[new.order][[1]]
    } else {
      # When comparing models fitted by maximum likelihood to the same data,
      # the smaller the AIC, BIC or MDL the better the fit.
      new.order <- order(criterion.value,decreasing=F)
      best <- candDist[new.order][1]
      est.par <- est.pars[new.order][[1]]
    }
    
    attributes(best)$best.dist.par <- est.par
    attributes(best)$criterion.value <- criterion.value
    
    return(best)
  }

