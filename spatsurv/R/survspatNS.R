##' survspatNS function
##'
##' A function to perform maximun likelihood inference for non-spatial survival data.
##'
##' @param formula the model formula in a format compatible with the function flexsurvreg from the flexsurv package 
##' @param data a SpatialPointsDataFrame object containing the survival data as one of the columns
##' @param dist choice of distribution function for baseline hazard. Current options are: exponentialHaz, weibullHaz, gompertzHaz, makehamHaz, tpowHaz
##' @param control additional control parameters, see ?inference.control
##' @return an object inheriting class 'mcmcspatsurv' for which there exist methods for printing, summarising and making inference from.
##' @seealso \link{tpowHaz}, \link{exponentialHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{weibullHaz}, 
##' \link{covmodel}, link{ExponentialCovFct}, \code{SpikedExponentialCovFct},
##' \link{mcmcpars}, \link{mcmcPriors}, \link{inference.control} 
##' @references 
##' \enumerate{
##'     \item Benjamin M. Taylor. Auxiliary Variable Markov Chain Monte Carlo for Spatial Survival and Geostatistical Models. Benjamin M. Taylor. Submitted. \url{http://arxiv.org/abs/1501.01665}
##' }
##' @export


survspatNS <- function( formula,
                        data,
                        dist,
                        control=inference.control()){
                        
    control$hessian <- TRUE                        
                     
    responsename <- as.character(formula[[2]])
    survivaldata <- data[[responsename]]
    checkSurvivalData(survivaldata) 
                     

    # start timing,
    
    start <- Sys.time()

    control$dist <- dist
           
    ##########
    # This chunk of code borrowed from flexsurvreg    
    ##########  
                      
    call <- match.call()
    indx <- match(c("formula", "data"), names(call), nomatch = 0)
    if (indx[1] == 0){ 
        stop("A \"formula\" argument is required")
    }
    temp <- call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    m <- eval(temp, parent.frame())

    Terms <- attr(m, "terms")
    X <- model.matrix(Terms, m)
    
    ##########
    # End of borrowed code    
    ##########
 
    X <- X[, -1, drop = FALSE]                               
              
    info <- distinfo(dist)()
    
    control$omegatrans <- info$trans
    control$omegaitrans <- info$itrans
    control$omegajacobian <- info$jacobian # used in computing the derivative of the log posterior with respect to the transformed omega (since it is easier to compute with respect to omega) 
    control$omegahessian <- info$hessian

    control$censoringtype <- attr(survivaldata,"type")

    if(control$censoringtype=="left" | control$censoringtype=="right"){
        control$censored <- survivaldata[,"status"]==0
        control$notcensored <- !control$censored

        control$Ctest <- any(control$censored)
        control$Utest <- any(control$notcensored)        
        
    }
    else{
        control$rightcensored <- survivaldata[,"status"] == 0
        control$notcensored <- survivaldata[,"status"] == 1
        control$leftcensored <- survivaldata[,"status"] == 2
        control$intervalcensored <- survivaldata[,"status"] == 3

        control$Rtest <- any(control$rightcensored)        
        control$Utest <- any(control$notcensored) 
        control$Ltest <- any(control$leftcensored)
        control$Itest <- any(control$intervalcensored)
    }

    #######
    
    cat("\n","Maximum likelihood using BFGS ...","\n")
    mlmod <- maxlikparamPHsurv(surv=survivaldata,X=X,control=control)
    estim <- mlmod$par
    print(mlmod)
    cat("Done.\n")

    end <- Sys.time()
    
    #browser() 
 

    retlist <- list()
    retlist$formula <- formula
    retlist$dist <- dist
    retlist$control <- control
    
    retlist$terms <- Terms
    retlist$mlmod <- mlmod

    # construct artificial samples from which the baseline hazard can be obtained with confidence intervals
    ch <- t(chol(solve(mlmod$hessian)))
    samp <- t(mlmod$par+ch%*%matrix(rnorm(1000*ncol(ch)),ncol(ch),1000))

    betasamp <- samp[,1:ncol(X),drop=FALSE]
    
    ####
    #   Back transform for output
    ####
    omegasamp <- samp[,(ncol(X)+1):ncol(samp),drop=FALSE]
    if(ncol(omegasamp)>1){
        omegasamp <- t(apply(omegasamp,1,control$omegaitrans))
    }
    else{
        omegasamp <- t(t(apply(omegasamp,1,control$omegaitrans)))
    }
    colnames(omegasamp) <- info$parnames
    
    colnames(betasamp) <- colnames(model.matrix(formula,data))[-1] #attr(Terms,"term.labels")
    retlist$betasamp <- betasamp
    retlist$omegasamp <- omegasamp

    retlist$Ysamp <- matrix(0,1000,nrow(X))
    
    #retlist$loglik <- loglik
    
    retlist$X <- X
    retlist$survivaldata <- survivaldata   
    retlist$gridded <- control$gridded    
    retlist$omegatrans <- control$omegatrans
    retlist$omegaitrans <- control$omegaitrans   
    retlist$control <- control
    retlist$censoringtype <- attr(survivaldata,"type")   
    retlist$time.taken <- Sys.time() - start
    
    cat("Time taken:",retlist$time.taken,"\n")
    
    class(retlist) <- c("list","mlspatsurv")

    return(retlist)
}


