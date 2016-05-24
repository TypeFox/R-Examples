##=============================================================================
## 'Prediction' and Inference on return levels using the so called 'delta
##  method'
##
## 'newdata' contains the return periods
##==============================================================================

predict.Renouv <- function(object,
                           newdata = c(10, 20, 50, 100, 200, 500, 1000),
                           cov.rate = TRUE,
                           level = c(0.95, 0.70),
                           prob = FALSE,
                           trace = 1,
                           eps = 1e-6,
                           ...) {
  
    if ( names(object$estimate)[1L] != "lambda") {
        stop("element 1 of 'estimate' must have name \"lambda\"")
    }
    
    pct.conf <- level * 100
    pred.period <- newdata
    
    nc <- 3L + 2L * length(pct.conf)
    
    ## Prepare a 'pred' matrix with dims and names
    
    cnames <- c("period", "prob", "quant",
                paste(rep(c("L", "U"), length(pct.conf)),
                      rep(pct.conf, each = 2L), sep = "."))
    
    pred.prob <- 1.0 - 1 / object$estimate[1L] / pred.period
    ind <- (pred.prob > 0.0) & (pred.prob < 1.0)
    pred.period <- pred.period[ind]
    pred.prob <- pred.prob[ind]

    p <- object$p
    lambda <-  object$estimate[1L]
    p.y <- object$p.y
    est.y <- object$estimate[-1L]
    fixed.y <- object$fixed.y
    cov.y <- object$cov[-1, -1, drop = FALSE]
    nb.OT <- object$nb.OT
    if (is.null(nb.OT)) nb.OT <- length(object$y.OT)
    
    pred <- matrix(NA, nrow = length(pred.period), ncol = nc)
    rownames(pred) <- pred.period
    colnames(pred) <- cnames
    
    ## compute quantiles
    y <- object$funs$q.y(parm = object$estimate[-1], p = pred.prob)
    pred[ , "period"] <- pred.period
    pred[ , "prob"]  <- pred.prob
    pred[ , "quant"] <- object$threshold + y
    
    if (object$transFlag) {
        threshold.trans <- object$funs$transfun(object$threshold)
        dth <- threshold.trans - object$threshold
        pred[ , "quant"] <- object$funs$invtransfun(dth +  pred[ , "quant"])
    }
    
    ##-------------------------------------------------------------------------
    ## Prepare if needed a matrix with rows matching the probabilities
    ## (or return periods) and columns mathing the estimated parameters
    ##-------------------------------------------------------------------------
    if (p > 0L) {
        Delta <- array(0, dim = c(length(pred.prob), p),
                       dimnames = list(round(pred.period, digits = 2L),
                           names(object$estimate)[!object$fixed]))
    }

    ##-------------------------------------------------------------------------
    ## Note that 'logf.y' MUST accept a vectorized 'x' formal 
    ##-------------------------------------------------------------------------
    if (!object$fixed[1L] && cov.rate) {
        f.y <- exp(object$funs$logf.y(parm = est.y, x = y))
        Delta[ , 1L] <- 1.0 / (lambda^2 * pred.period * f.y)
    }
    
    if ( (object$distname.y == "exponential") && !object$history.MAX$flag
        && !object$history.OTS$flag ) {
        
        if (trace) {
            cat("Special inference for the exponential case without history\n")
        }

        if (cov.rate) {
            warning("uncertainty on the rate not taken into account yet",
                    "  in the exponential with no history case")
        }
        
        ##---------------------------------------------------------------------
        ## Use the sampling distribution to derive confidence
        ## limits on quantiles. If 'lambda' was known, these limits
        ## would be exact.
        ##---------------------------------------------------------------------
        
        for (ipct in 1L:length(pct.conf)) {
            alpha.conf <- (100 - pct.conf[ipct]) / 100
           
            theta.L <- 2 * nb.OT / est.y /
                qchisq(1 - alpha.conf / 2, df = 2 * nb.OT)
            theta.U <- 2 * nb.OT / est.y /
                qchisq(alpha.conf / 2, df = 2 * nb.OT)
            
            pred[ , 2L * ipct + 2L] <-
                object$threshold + qexp(p = pred.prob, rate = 1.0 / theta.L)
            pred[ , 2L * ipct + 3L] <-
                object$threshold + qexp(p = pred.prob, rate = 1.0 / theta.U)  
        }
        
        if (object$transFlag) {
            ## dth was defined above
            for (ipct in 1L:length(pct.conf)) { 
                pred[ , 2L * ipct + 2L] <-
                    object$funs$invtransfun(dth + pred[ , 2L * ipct + 2L])
                pred[ , 2L * ipct + 3L] <-
                    object$funs$invtransfun(dth + pred[ , 2L * ipct + 3L]) 
            }
        }
        
        infer.method <-
            "chi-square for exponential distribution (no historical data)"
        
    } else {
        
        ##=====================================================================
        ## DELTA METHOD
        ## matrices for numerical derivation
        ## The column 'ip' of 'Parmat' contains the parameter value with
        ## a tiny modification of its ip component.
        ##=====================================================================
        
        if (p.y > 0L) { 
            
            parmMat <- matrix(est.y[!fixed.y], nrow = p.y, ncol = p.y)
            rownames(parmMat) <- object$parnames.y[!fixed.y]
            colnames(parmMat) <- object$parnames.y[!fixed.y]
            
            dparms <- abs(est.y[!fixed.y]) * eps
            dparms[dparms < eps] <- eps
            
            for (ip in 1L:p.y) parmMat[ip, ip] <- parmMat[ip, ip] + dparms[ip]
            
            ##=================================================================
            ## Compute Delta's elts : derivatives w.r.t. unknown params
            ##
            ## CAUTION
            ##
            ## The quantile function 'q.y' takes a vector as first arg which
            ## should morally be of length 'p.y', but is here of length
            ## parnb.y. This works because 'q.y' uses the elements in named
            ## form, irrespective of their position.
            ##
            ##=================================================================
            shift <- !object$fixed[1L]
            
            for (i in 1L:length(pred.prob)) {
                for (ip in 1L:p.y) {
                    est.prov <- est.y
                    est.prov[!fixed.y] <-  parmMat[ , ip]
                    Delta[i, ip + shift] <- (object$funs$q.y(est.prov, p = pred.prob[i]) -
                                                 y[i]) / dparms[ip]
                }
            }
            
        }

        if (trace >= 2L) {
            cat("Derivative of return levels w.r.t. parameters\n")
            print(Delta)
        }
        
        pred.sig <- rep(0, length(pred.prob))
        if (p > 0L) {
            for (i in 1L:length(pred.prob)) {
                dd <- Delta[i, ]
                pred.sig[i] <- sqrt(t(dd) %*% object$cov[!object$fixed, !object$fixed] %*% dd)
            }
        }       
                                  
        ##=====================================================================
        ## Apply "delta method" for the quantiles
        ##
        ## Note that 'q.y' MUST accept a vectorized prob 
        ## Use a loop on i to remove this constraint ???
        ##=====================================================================
        
        for (ipct in 1L:length(pct.conf)) {
            alpha.conf <- (100 - pct.conf[ipct]) / 100
            z.conf <- qnorm(1 - alpha.conf / 2.0)
            pred[ , 2L * ipct + 2L] <-
                object$threshold + object$funs$q.y(est.y, p = pred.prob) -
                    z.conf * pred.sig
            pred[ , 2L * ipct + 3L] <-
                object$threshold + object$funs$q.y(est.y, p = pred.prob) +
                    z.conf * pred.sig
        }
        
        if (object$transFlag) {
            for (ipct in 1L:length(pct.conf)) {
                pred[ , 2L * ipct + 2L] <-
                    object$funs$invtransfun(dth  + pred[ , 2L * ipct + 2L])
                pred[ , 2L * ipct + 3L] <-
                    object$funs$invtransfun(dth  + pred[ , 2L * ipct + 3L])
            }
        }
        
        infer.method <- "Delta method"
    }
    
    pred <- as.data.frame(pred)
    if (!prob) {
        pred$prob <- NULL
    }
    attr(pred, "infer.method") <- infer.method
    pred
    
}
