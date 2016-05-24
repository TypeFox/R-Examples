#' @title Parameter Estimation Evaluation.
#' @description A function to evaluate the parameter estimation function.
#' 
#' @rdname eval.estimation
#' @name eval.estimation
#' 
#' @param rdist Random variable generating function.
#' @param edist Parameter estimation function.
#' @param n Sample size.
#' @param rep.num Number of replicates. 
#' @param params True parameters of the distribution.
#' @param method Estimation method.
#'
#' @return A list containing the mean and sd of the estimated parameters.\cr
#' \cr
#'  na.cont returns the number of "na"s that appeared in the parameter estimation.
#'
#' @author Haizhen Wu and A. Jonathan R. Godfrey.
#'
#' @examples
#' eval.estimation(rdist=rNormal,edist=eNormal,n = 100, rep.num = 50, 
#' params = list(mean = 1,sd = 5))

#' @export eval.estimation
eval.estimation <- function(rdist, edist, n = 20, rep.num = 1e3, params, method = "numerical.MLE"){
    k <- length(params)
    
    est.par <- array(NA, dim = c(rep.num,k))
    start.time <- proc.time()
    for(i in 1:rep.num){
      X <- rdist(n=n,params=params)
      est.par[i,] <- as.numeric(unlist(edist(X,method=method)))
      #       print(paste("i=",i))
    }
    return(list(method = method,
                est.mean = apply(est.par,2, mean, na.rm =T), 
                est.sd = apply(est.par,2, sd, na.rm =T), 
                time = proc.time() - start.time,
                na.cont = sum(is.na(est.par[,1])))
    )
  }
