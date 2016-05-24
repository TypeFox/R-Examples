# Copyright 2014 by Dmitry Pavlyuk <Dmitry.V.Pavlyuk@gmail.com>

#
# Routines for maximum likelihood estimator
#

#
# Wrapper for 'optim' procedure from 'stats' package
#
optimEstimator <- function(formula, data, func.lf, ini,gr=NULL){
    ret <- prepareEstimates(status = 1000)
    if (is.null(ini)){
        logging("Ini is not defined",caller = match.call())
        return(prepareEstimates(status = 1001))
    }
    tryCatch({
        control <- envirGet("control")
        constr <- envirGet("constr")
        pscale <- envirGet("parscale")
        st <- 0
        stMax <- control$optim.control$repeatNM
        if (is.null(stMax)) stMax <- 1
        repeat{
            st <- st + 1
            p<-maxLik(func.lf, start=ini,grad=gr,constraints=constr,
                 tol=control$optim.control$tol,
                 #parscale=pscale,
                 iterlim=control$optim.control$iterlim)
            
            if (!is.null(gr) || st==stMax){
                break;
            }else{
                logging(p, level = "info")
                logging(paste("Restarting Nelder-Mead: ", st), level = "info")
                ini <- p$estimate
            }
        }
        #p <- psoptim(par=NA, fn=func.lf, lower=c(-Inf, -Inf, -Inf, 0,0,-1), upper=c(Inf, Inf, Inf, Inf,Inf,1))
        logging(summary(p), level = "debug")
        logging(paste("Convergence is",
                      ifelse(succcessCode(p)==0, "", "NOT"),"achieved",
                      "[",p$code,"]"),
                level = "info",
                caller = match.call())    
        logging(paste("Log-Likelihood value = ",p$maximum),
                level = "info",
                caller = match.call())  
        ret <- prepareEstimates(estimates = p)
    }, error = function(e){ 
        logging(e$message, level="warn")
    })
    return(ret)
}

succcessCode <- function(p){
    res <- p$code
    if (p$type=="Newton-Raphson maximisation"){
        if (res<3) res<-0
    }
    return(res)
}
#
# Creates the ModelEstimates class on the base of raw (optim) estimation results
#
prepareEstimates <- function(estimates = NULL, status = 0){
    ret <- new("ModelEstimates", status = status)
    if (!is.null(estimates)) {
        if (estimates$maximum == Infin){
            return(new("ModelEstimates", status = 1002))
        }
        ret <- new("ModelEstimates",
                   resultParams = estimates$estimate,
                   status = succcessCode(estimates),
                   logL = estimates$maximum,
                   logLcalls = estimates$iterations,
                   hessian = estimates$hessian
        )
        #Setting hessian and implicitly calculating standard errors
        #if (!is.null(estimates$hessian)){
        #    hessian(ret) <- estimates$hessian
        #}
    }
    return(ret)
}