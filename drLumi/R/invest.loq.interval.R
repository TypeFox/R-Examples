## Inverse of the response variable.
##
## @param model nls class model
## @param yvalue the MFI value to be applied the inverse
## @param inter.method interval method to be used for the inverse estimation 
## of the confidence intervals
## @param level level value of the interval method
## @param maxiter maximum number of iterations of the uniroot function
## @param tol tolerance value of the uniroot function
## 
## @details This function estimates the inverse of a value given a nls model
## based on the uniroot function. The SE is an aproximation based on the
## differences between the two concentration confidence intervals values 
## estimated.
invest.loq.interval <- function(x, analyte, yvalue,  inter.method="prediction", 
            level = 0.95, maxiter=1000,  tol=.Machine$double.eps^0.25 ){

    if(!inherits(x,"scluminex")) stop("model object must be 'scluminex' class")
    tinterval <- charmatch(inter.method, c("confidence","prediction"))
    if(is.na(tinterval)){
        stop("inter.method argument must be 'confidence' or 'prediction'")  
    }
    if(level<0 | level>1){
        stop("'level' argument must be a value between 0 and 1")  
    }
    
    model <- x[[analyte]]$model
    env.mod <- model$m$getEnv()
    allobj <- ls(env.mod)
    lhs <- as.character(model$m$formula()[[2]])
    pars <- names(model$m$getPars())
    rhs <- allobj[-match(c(pars,lhs),allobj)]
    ndf <- data.frame(get(lhs, model$m$getEnv()), get(rhs, model$m$getEnv()))
    names(ndf) <- c(lhs, rhs)
    lower <- min(ndf[,rhs])
    upper <- max(ndf[,rhs])
    ff <- function(xx, fest.type, fyvalue) {
        ans <- conf_bands(x, analyte, xx, level, inter.method)[,fest.type] - fyvalue
        ans
    } 
    estimation.inv <- sapply(yvalue, function(i) 
        lapply(c(paste0(lhs,".lci"),paste0(lhs,".uci"), lhs), 
            function(est.type) try(uniroot(ff,interval = c(lower, upper),
                                    fyvalue=i,fest.type=est.type,
                                    tol = tol, 
                                    maxiter = maxiter)$root,
                                    silent=TRUE)))  
        qc.inv <-  lapply(estimation.inv, function(x){
        if(inherits(x,"try-error")) x <- NA
        x
    })  
    mat.inv <- matrix(unlist(qc.inv), ncol = 3, byrow=TRUE)  
    pred.lower <- conf_bands(x, analyte, lower)[,1]
    pred.upper <- conf_bands(x, analyte, upper)[,1]  
    if (pred.lower<=pred.upper) {
        aux <- mat.inv[,2]
        inv.uci <- mat.inv[,1]
        inv.lci <- aux
        fit <- mat.inv[,3]    
    } else {
        inv.lci <- mat.inv[,1]
        inv.uci <- mat.inv[,2]
        fit <- mat.inv[,3] 
    }   
    tvalue <-  qt(1-((1-level)/2), df = summary(model)$df[2])    
    ase <- try( abs(inv.uci - inv.lci)/2 /tvalue, silent=TRUE)
    if(inherits(ase,"try-error")) ase <- NA  
    ans <- data.frame(cbind(yvalue, fit, inv.lci, inv.uci,  ase))
    names(ans) <- c(lhs, paste0(rhs,".fit"),
                    paste0(rhs,".lci"), 
                    paste0(rhs,".uci"),
                    paste0(rhs,".aprox.se"))
    return(ans)
}
