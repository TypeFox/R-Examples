
##=============================================================================
## default probabilities at return levels evaluation
##=============================================================================

.rl.prob <- c(0.0001,
              seq(from = 0.01, to = 0.09, by = 0.01),
              seq(from = 0.10, to = 0.80, by = 0.10),
              seq(from = 0.85, to = 0.99, by = 0.01),
              0.995, 0.996, 0.997, 0.998, 0.999,
              0.9995, 0.9996, 0.9997, 0.9998, 0.9999,
              0.99999, 0.999999)

`RenouvNoEst` <- function(threshold,
                          estimate = NULL,
                          distname.y = "exponential",
                          fixed.par.y = NULL,
                          trans.y = NULL,
                          pct.conf = c(95, 70),
                          rl.prob = NULL,
                          prob.max = 1.0 - 1e-4,
                          pred.period = NULL,
                          cov = NULL,
                          nb.OT = NULL,
                          infer.method = NULL) {
    
    mc <- match.call()
    if (!is.null(nb.OT)) {
        if ( is.na(nb.OT) || (nb.OT != as.integer(nb.OT)) || nb.OT < 1L )
            stop("when given, 'nb.OT' must be a positive integer")
    }
    
    if (is.null(fixed.par.y)) fixed.par.y <- list()

    ## drop all attributes except names
    ne <- names(estimate)
    estimate <- as.vector(estimate)
    names(estimate) <- ne
    
    ## make sure that 'lambda' is given!
    pos <- match("lambda", names(estimate)) 
    if (length(pos) != 1L)  stop("'estimate' must have exactly one element named \"lambda\"")
    parnames.cand <- names(estimate)[-pos]
    pnames <- c("lambda", parnames.cand)
    
    ## check the distribution name
    myDist <- checkDist(distname.y = distname.y)
    
    if (!myDist$special.y) myDist$parnames.y <- parnames.cand
    
    ##=========================================================================
    ## build probability functions and find the characteristics
    ## of the parameters
    ##=========================================================================
    myFuns <- makeFuns(funname.y = myDist$funname.y,
                       parnames.y = myDist$parnames.y,
                       fixed.par.y = myDist$fixed.par.y,
                       trace = 0) 
    
    ## build transformation functions, if necessary
    myTransFuns <- transFuns(trans.y = trans.y,
                             distname.y = distname.y) 
    
    ## clean and complete myFuns, if necessary
    funs <- list(transfun = myTransFuns$transfun,
                 invtransfun = myTransFuns$invtransfun,
                 dfun.y = myFuns$dfun.y,
                 pfun.y = myFuns$pfun.y,
                 qfun.y = myFuns$qfun.y,
                 logf.y = myFuns$logf.y,
                 q.y = myFuns$q.y, 
                 F.y = myFuns$F.y)
    
    ## Check that the provided covariance is good
    p <- myFuns$p.y + 1L
    covAll <- array(0, dim = c(p, p), dimnames = list(pnames, pnames))
    fixed <-  c("lambda" = FALSE, myFuns$fixed.y)
    
    if (!is.null(cov)) {
        if (!is.matrix(cov) || (nrow(cov) != ncol(cov))) {
            stop("'cov' must be a square matrix")
        }
        if (is.null(rownames(cov))) {
            if (is.null(colnames(cov))) {
                stop("'cov' must have rownames or colnames")
            }
            covNames <- colnames(cov)
        } else  covNames <- rownames(cov)
        m <- pnames[pnames %in% covNames]
        covAll[m, m] <- cov[m, m]
    }
    cov <- covAll
    sigma <- sqrt(diag(cov))
    
    ##=========================================================================
    ## Prepare results list
    ##=========================================================================
    
    res <- list(call = mc,
                nb.OT = nb.OT,
                threshold = threshold,
                distname.y = myDist$distname.y,
                p.y = myFuns$p.y,
                p = myFuns$p.y + 1L,
                df = myFuns$p.y + 1L,
                nobs = ifelse(is.null(nb.OT), NA, nb.OT),
                parnames.y = myDist$parnames.y,
                trans.y = NULL,
                sigma = sigma,
                cov = cov,
                estimate = estimate,
                fixed.y = myFuns$fixed.y,
                fixed = fixed, 
                history.MAX = list(flag = FALSE),
                history.OTS = list(flag = FALSE),
                transFlag = myTransFuns$transFlag,
                funs = funs,
                pct.conf = pct.conf)
    
    ## Compute a return level table
    if (is.null(rl.prob)) {
        rl.prob <- .rl.prob
        rl.prob <- rl.prob[rl.prob <= prob.max]
    } else {
        if (any(is.na(rl.prob))) stop("'rl.prob' values can not be NA") 
        if ( any(rl.prob <= 0.0) || any(rl.prob >= 1.0) ) stop("'rl.prob' values must be >0 and <1") 
        rl.prob <- sort(rl.prob)
    }
    
    if (is.null(pred.period)) {
        rr <- 2
        pred.period <- (10^rr) * c(0.1, 0.2, 0.5, 1:10)
    } else {
        if (any(is.na(pred.period))) stop("'pred.period' values can not be NA") 
        pred.period <- sort(pred.period)
    }
    
    rl.period <- 1 / estimate[1] / (1 - rl.prob)
    rl.sort <- sort(c(rl.period, pred.period), index.return = TRUE)
    rl.period <- rl.sort$x
    ind <- !duplicated(rl.period)
    rl.period <- rl.period[ind]
    ## rl.prob <- rl.prob[ind]
    
    ind.pred <- rl.period %in% pred.period
    ## cat("rl.period\n"); print(rl.period)
      
    ## Use pct.conf 
    ret.lev <- predict.Renouv(object = res,
                              newdata = rl.period,
                              cov.rate = FALSE,
                              level = round(pct.conf) / 100,
                              trace = 1,
                              eps = 1e-6)
    
    rownames(ret.lev) <- NULL
    res[["ret.lev"]] <- ret.lev
    res[["pred"]] <- as.data.frame(ret.lev[ind.pred, , drop = FALSE])
    res[["infer.method"]] <- attr(res[["pred"]], "infer.method")
    
    class(res) <- "Renouv"
    res 
    
}
