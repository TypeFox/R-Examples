#' plots \eqn{\beta^{rc}} given hyperpriori parameters
#' 
#' prioriPlot simulates the \eqn{beta_i^rc} values at first level given specific parameters at hyperpriori level
#' @param pars list of parameters for hyperpriori. 
#' if which="gamma" then parameter has to be a list with shape and rate as parameters
#' if which="expo" then parameter has to be a list with only lam
#' @param which specified priori. "gamma" or "expo"
#' @param cols integer specifying how many columns the RxC-Table should have
#' @param alphaSample integer specifying the number of times new alpha-values are drawn
#' @param betaSample integer specifying the number of times betas will be drawn for each alpha-value
#' @param plot logical TRUE/FALSE if histogram should be plotted
#' @param ... additional arguments for "hist" function
#' 
#' @details
#' Calculation is made via the marginal beta distribution
#' 
#' function structure:
#' \itemize{
#'   \item \code{"gamma"} choose one parameter for every alpha_rc-parameter or a two matrices
#'                       of parameters specifying lambda's for every alpha_rc-parameter
#'   \item \code{"expo"} choose one parameter for every alpha_rc-parameter or a one matrix
#'                       of parameters specifying lambda's for every alpha_rc-parameter
#' }
#' 
#' @return
#' nested \code{list} with each element containing another \code{list}.
#' First level are rows and second level are columns per row.
#' 
#' @examples
#' \dontrun{
#' test1 <- prioriPlot(list(shape=4,rate=2), "gamma",cols=4)
#' str(test1)
#' 
#' pars <- list(shape=matrix(1:9,3,3),rate=matrix(9:1,3,3))
#' test2 <- prioriPlot(pars, "gamma",breaks=100)
#' test3 <- prioriPlot(list(shape=8,rate=2),"gamma",breaks=100,cols=3)
#' 
#' pars4 <- list(shape=matrix(c(6,6,6),1,3), rate=matrix(c(4,4,4),1,3))
#' test4 <- prioriPlot(pars4, "gamma",breaks=100)
#' 
#' pars5 <- list(lam=2)
#' test5 <- prioriPlot(pars5, "expo",cols=4, breaks=100)
#' 
#' pars6 <- list(lam=matrix(1:9,3,3)/100)
#' test6 <- prioriPlot(pars6, "expo", breaks=25, col=grey(0.8))
#' 
#' # example for 3x4-table
#' set.seed(568)
#' pars7 <- list(shape=matrix(sample(1:20,12), 3,4), rate=matrix(sample(1:20,12),3,4))
#' test7 <- prioriPlot(pars7, "gamma",breaks=50)
#' }
#' 
#' @export
#' 


prioriPlot <- function(pars, which, cols,alphaSample =10000, betaSample=300,plot=TRUE,...){
  if(!which %in% c("gamma", "expo"))
    stop("wrong priori specified", call.=FALSE)
  #############################################################################
  ## Gamma-distribution
  #############################################################################  
  if(which=="gamma"){
    if( !all(c("shape", "rate") %in% names(pars)) )
      stop("no element called shape/rate in parameters!", call.=FALSE)
    
    if(is.matrix(pars$shape) & is.matrix(pars$rate)){
      ############################ Fall mit  mehreren rates/shapes
      if(!all(pars$shape>0) | !all(pars$rate>0))
        stop("priori parameters have to be positive number!", call.=FALSE)
      if(!all(dim(pars$shape)==dim(pars$rate)))
        stop("shape and rate parameters must have same dimensions!", call.=FALSE)
      
      r <- nrow(pars$shape)
      c <- ncol(pars$shape)
      
      par(mfrow=c(r,c))
      betaList <- vector("list",r)
      for(rr in 1:r){
        tmpList <- vector("list",c)
        for(cc in 1:c){
          cat("shape: ", pars$shape[rr,cc], ", rate: ", pars$rate[rr,cc], "\n")
          betas <- sapply(1:alphaSample, function(as){
              # Ziehen für die ganze reihe
            alpha_rc <- sapply(1:c, function(j) rgamma(1,pars$shape[rr,j],pars$rate[rr,j]))
              # spezifische column wählen
              # restterm wird berechnet
            alpha_rcsum <- sum(alpha_rc) - alpha_rc[cc]
            rbeta(betaSample, alpha_rc[cc], alpha_rcsum)
          })
          tmpList[[cc]] <- c(betas)
          if(plot)
            hist(c(betas),xlim=c(0,1),xlab=expression(beta^rc),
                 main=paste0("Gamma(shape: ",pars$shape[rr,cc]," rate: ",pars$rate[rr,cc],")\n r=",rr," c=",cc),...) 
        }
        betaList[[rr]]<-tmpList
      }
      par(mfrow=c(1,1))
      
      return(betaList)
            
    } else{ 
      ############################ Fall mit nur einer rate/shape
      if(pars$shape<=0 | pars$rate <=0)
        stop("gamma parameters have to be positive numbers!", call.=FALSE)
      
      betas <- sapply(1:alphaSample, function(as){
        alpha_rc <- rgamma(cols, shape=pars$shape, rate=pars$rate)
        alpha_rcsum <- sum(alpha_rc) - alpha_rc[1]
        rbeta(betaSample, alpha_rc[1], alpha_rcsum)
      })
      if(plot){
        hist(c(betas),xlim=c(0,1),xlab=expression(beta^rc),
             main=paste0("Gamma(shape: ", pars$shape, ", rate: ", pars$rate, ")\n columns: ", cols),...)
      }
      return(betas)
    }
  }
  #############################################################################
  ## Exponential-distribution
  #############################################################################
  if(which=="expo"){
    if(!"lam" %in% names(pars))
      stop("no element called lam in parameters!", call.=FALSE)
    
    if(is.matrix(pars$lam)){
      if(!all(pars$lam>0))
        stop("priori parameters have to be positive number!", call.=FALSE)
    
      r <- nrow(pars$lam)
      c <- ncol(pars$lam)
      
      par(mfrow=c(r,c))
      betaList <- vector("list",r)
      for(rr in 1:r){
        tmpList <- vector("list",c)
        for(cc in 1:c){
          cat("lambda: ", pars$lam[rr,cc], "\n")
          betas <- sapply(1:alphaSample, function(as){
            # Ziehen für die ganze reihe
            alpha_rc <- sapply(1:c, function(j) rexp(1,pars$lam[rr,j]))
            # spezifische column wählen
            # restterm wird berechnet
            alpha_rcsum <- sum(alpha_rc) - alpha_rc[cc]
            rbeta(betaSample, alpha_rc[cc], alpha_rcsum)
          })
          tmpList[[cc]] <- c(betas)
          
          if(plot)
            hist(c(betas),xlim=c(0,1),xlab=expression(beta^rc),
                 main=paste0("Exp(lambda: ", pars$lam[rr,cc],")\n r=", rr, " c=", cc),...)
            
        }
        betaList[[rr]]<-tmpList
      }
      par(mfrow=c(1,1))
      
      return(betaList)
      
    } else{
      ############################ Fall mit nur einem lambda
      if(pars$lam<=0)
        stop("exponential parameter has to be positive number!", call.=FALSE)
      
      betas <- sapply(1:alphaSample, function(as){
        alpha_rc <- rexp(cols, rate=pars$lam)
        alpha_rcsum <- sum(alpha_rc) - alpha_rc[1]
        rbeta(betaSample, alpha_rc[1], alpha_rcsum)
      })
      if(plot){
        hist(c(betas),xlim=c(0,1),xlab=expression(beta^rc),
             main=paste0("Exp(lambda: ", pars$lam, ")\n columns: ", cols),...)
      }
      return(betas)
    }
  } 
}



