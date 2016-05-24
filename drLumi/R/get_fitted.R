## Get fitted values
get_fitted<-function(model, bkg.method, bkg.mean, fct) {
    lhs <- as.character(model$m$formula()[[2]])
    parameters <- names(model$m$getPars())
    allobj <- ls(model$m$getEnv())
    rhs <- allobj[-match(c(parameters,lhs),allobj)]
    ndf <- data.frame(get(lhs, model$m$getEnv()), get(rhs, model$m$getEnv()))
    names(ndf) <- c(lhs, rhs)
    
    if(bkg.method=="constraint")  log10.bkgmean <- log10(bkg.mean)
 
    yvalue <- ndf[,lhs]    
    if(bkg.method!="constraint"){
      inv <- invest.fun(model,"noconstraint", fct, yvalue, parameters, NULL)
    } 
    if(bkg.method=="constraint"){
      inv <- invest.fun(model,"constraint", 
                        fct, yvalue, parameters, log10.bkgmean)
    }
    est <- unlist(lapply(1:length(yvalue), function(x) inv$inv[[x]]$est))

    form <- unlist(lapply(1:length(yvalue), function(x) inv$form[[x]]))
    
    se <- lapply(form,
                 function(x) msm::deltamethod(x, coef(model), vcov(model)))  
    se <- unlist(se)
  
    ans <- as.data.frame(cbind(est, se))
    names(ans) <- c("log10_concentration.fit", "log10_concentratrion.se")
    return(ans)
}



