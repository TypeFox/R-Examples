.getDistCode<- function(dist){
  if(dist == "exponential"){
    .getDistCode <- 1
  } else if(dist == "weibull"){ 
    .getDistCode <- 2
  } else if(dist == "burr"){ 
    .getDistCode <- 3
  } else if(dist == "gengamma"){ 
    .getDistCode <- 4
  } else if(dist == "genf"){ 
    .getDistCode <- 5
  } else if(dist == "qweibull"){ 
    .getDistCode <- 6
  } else if(dist == "mixqwe"){ 
    .getDistCode <- 7
  } else if(dist == "mixqww"){ 
    .getDistCode <- 8
  } else if(dist == "mixinvgauss"){ 
    .getDistCode <- 9
  } else if(dist == "birnbaum-saunders"){ 
    .getDistCode <- 10 
  } else stop("the provided distribution does not exist")
}

.setBP <- function(Nbp){
  .setBP <- switch(as.character(Nbp),
                  "1" = c(1),
                  "2" = c(.5, 1.5),
                  "3" = c(.3, .8, 1.5),
                  "4" = c(.3, .8, 1.5, 2),
                  "5" = c(.3, .5, 1, 1.5, 2),
                  "6" = c(.2, .4, .6, 1, 1.5, 2))
}


.checkOrderAndPara <- function(order, para, distCode, model){
  #checks the order:
  if(model %in% c("ACD", "LACD1", "LACD2","ABACD")){
    if(length(order) != 2) stop("The order is not entered in the correct format, check the description")
  } else if(model %in% c("AMACD", "SNIACD", "LSNIACD")){
    if(length(order) != 3) stop("The order is not entered in the correct format, check the description")  
  } else stop("the model is wrongly entered or not supported!")
  
  if(any(order != round(order))) stop("The order must be integers")  
  if(any(order < 0)) stop("The order can't have negative entries")   
  
  #checks the number of parameters
  if(model %in% c("ACD", "LACD1", "LACD2")){
    Nmodelpara <- 1 +  order[1] + order[2]
  } else if(model %in% c("ABACD")){
    Nmodelpara <- 4 + 2*order[1] + order[2]
  }else if(model %in% c("AMACD")){
    Nmodelpara <- 1 +  order[1] + order[2] + order[3]
  } else if(model %in% c("BACD")){
    Nmodelpara <- 1 +  order[1] + order[2] + 2
  } else if(model %in% c("SNIACD", "LSNIACD")){
    Nmodelpara <- 1 +  order[1] + order[2] + order[3]
  }
    
  if(distCode == 1){
    Ndistpara <- 0
  } else if(distCode == 2){
    Ndistpara <- 1
  } else if(distCode == 3){
    Ndistpara <- 2
  } else if(distCode == 4){
    Ndistpara <- 2
  } else if(distCode == 5){
    Ndistpara <- 3
  } else if(distCode == 6){
    Ndistpara <- 2
  } else if(distCode == 7){
    Ndistpara <- 4
  } else if(distCode == 8){
    Ndistpara <- 5
  } else if(distCode == 9){
    Ndistpara <- 3
  } else if(distCode == 10){
    Ndistpara <- 1
  } else{
    stop("Wrong distCode!")
  } 
  
  if(length(para) != (Nmodelpara + Ndistpara)){
    errMsg <- paste("Wrong number of given parameters. The", model, "model should have", Nmodelpara, "and the distribution", Ndistpara, "parameters.")
    stop(errMsg)
  }   
}

.setStartPara <- function(model, distCode, mean, order){
    
  if(distCode == 1){
    distStartPara <- NULL
  } else if(distCode == 2){
    distStartPara <- .8
  } else if(distCode == 3){      
    distStartPara <- c(1.1, 0.3)  
  } else if(distCode == 4){      
    distStartPara <- c(2, .5)  
  } else if(distCode == 5){      
    distStartPara <- c(0.4, 0.7, 2.5)  
  } else if(distCode == 6){      
    distStartPara <- c(.8, 1.2)  
  } else if(distCode == 7){      
    distStartPara <- c(.7, 1.3, 1.2, 1.2)  
  } else if(distCode == 8){      
    distStartPara <- c(.7, 1.3, 1.2, 0.7, 1.2)  
  } else if(distCode == 9){      
    distStartPara <- c(.7, 0.2, 0.4)  
  } else if(distCode == 10){      
    distStartPara <- c(1)  
  } 
  
  if(model == "ACD"){
    startPara <- c(mean/10, rep(.15/order[1],order[1]), rep(.8/order[2],order[2]))
  } else if(model == "LACD1"){    
    startPara <- c(0.03, rep(.03/order[1],order[1]), rep(.98/order[2],order[2]))
  } else if(model == "LACD2"){    
    startPara <- c(0, rep(.03/order[1],order[1]), rep(.98/order[2],order[2]))
  } else if(model == "ABACD"){
    startPara <- c(mean/20, rep(.10/order[1],order[1]), rep(0, order[1]), rep(.8/order[2],order[2]), 0, 1, 1)
  } else if(model == "AMACD"){
    startPara <- c(mean/10, rep(.15/(order[1]+order[2]), (order[1]+order[2])), rep(.8/order[3],order[3]))
  } else if(model == "BACD"){
    startPara <- c(mean/20, rep(.10/order[1],order[1]), rep(.8/order[2],order[2]), 1, 1)
  } else if(model %in% c("SNIACD")){
    startPara <- c(mean/10, c(0.15, rep(.1,order[3])), rep(0,order[1] - 1), rep(.8/order[2],order[2]))
  } else if(model %in% c("LSNIACD")){
    startPara <- c(0, c(0.03, rep(0,order[3])), rep(0,order[1] - 1), rep(.8/order[2],order[2]))
  }
  
  return(list(startPara = c(startPara, distStartPara), modelStartPara = startPara, distStartPara = distStartPara))
}

.seperateStartPara <- function(startPara, model, distCode, order){  
   
  if(model == "ACD"){
    startMPara <- startPara[1:(1 + order[1] + order[2])]
  } else if(model == "LACD1"){    
    startMPara <- startPara[1:(1 + order[1] + order[2])]
  } else if(model == "LACD2"){    
    startMPara <- startPara[1:(1 + order[1] + order[2])]
  } else if(model == "ABACD"){
    startMPara <- startPara[1:(4 + 2 * order[1] + order[2])]
  } else if(model == "AMACD"){
    startMPara <- startPara[1:(1 + order[1] + order[2] + order[3])]
  } else if(model == "BACD"){
    startMPara <- startPara[1:(1 + order[1] + order[2] + 2)]
  } else if(model %in% c("SNIACD", "LSNIACD")){
    startMPara <- startPara[1:(1 + order[1] + order[2] + order[3])]
  }
  
  
  if(distCode == 1){
    distStartPara <- NULL
  } else{
    distStartPara <- startPara[(length(startMPara) + 1):length(startPara)]
  }
  
  return(list(startPara = startPara, modelStartPara = startMPara, distStartPara = distStartPara))
}

.checkOrder <- function(order, model){
  if(model %in% c("ACD", "LACD1", "LACD2","ABACD", "BACD")){
    if(length(order) != 2) stop("The order is not entered in the correct format, check the description")
  } else if(model %in% c("AMACD", "SNIACD", "LSNIACD")){
    if(length(order) != 3) stop("The order is not entered in the correct format, check the description")  
  } else stop("the model is wrongly entered or not supported!")
  
  if(any(order != round(order))) stop("The order must be integers")  
  if(any(order < 0)) stop("The order can't have negative entries")   
}

.setOrder <- function(model){
    if(model %in% c("ACD", "LACD1", "LACD2","ABACD", "BACD")){
      return(c(1, 1))
    } else if(model %in% c("AMACD")){
      return(c(1, 1, 1))
    }  else if(model %in% c("SNIACD", "LSNIACD")){
      return(c(1, 1, 2))
    }  
}

.getNewDay <- function(time){
  daysDiff <- as.Date(time[-1])-as.Date(time[1:(length(time)-1)])
  return(which(daysDiff != 0)+1)
}


#returns the full parameter vector from the shorter freePara and fixedParam
.returnfixedPara <- function(freePara, fixedParam, fixedParamPos){
  if(length(fixedParamPos) != length(freePara) + length(fixedParam)) stop(".returnfixedPara() error")
    
  returnPara <- numeric(length(fixedParamPos))
  fixedParamIndex <- 1
  fitParIndex <- 1  
  
  for(j in seq_along(fixedParamPos)){
    if(fixedParamPos[j]){
      returnPara[j] <-fixedParam[fixedParamIndex]
      names(returnPara)[j] <- names(fixedParam)[fixedParamIndex]
      fixedParamIndex <- fixedParamIndex + 1
    } else{
      returnPara[j] <- freePara[fitParIndex]
      names(returnPara)[j] <- names(freePara)[fitParIndex]
      fitParIndex <- fitParIndex + 1      
    }
  }
  
  return(returnPara)
}

#returns the full vector of SE from the shorter freeSE (sets the SE of fixed parameters to NA)
.returnfixedSE <- function(freeSE, fixedParamPos){
  if(sum(!fixedParamPos) != length(freeSE)) stop(".returnfixedSE() error")
  
  returnSE <- numeric(length(fixedParamPos))
  fixedParamIndex <- 1
  fitParIndex <- 1  
  
  for(j in seq_along(fixedParamPos)){
    if(fixedParamPos[j]){
      returnSE[j] <- NA
    } else{
      returnSE[j] <- freeSE[fitParIndex]
      fitParIndex <- fitParIndex + 1      
    }
  }
  
  return(returnSE)
}


.returnFixedMeanPara <- function(distCode, distPara){
  
  if(distCode == 1){ #Exponential    
    
    lambda <- 1
    names(lambda) <- "lambda"
    return(lambda)       
    
  } else if(distCode == 2){ #Weibull
    
    theta <- gamma(1+1/distPara)^distPara
    names(theta) <- "theta"
    return(theta)    
    
  } else if(distCode == 3){ #Burr
    
    kappa <- distPara[1]
    sig2 <- distPara[2]
    theta <- ((gamma(1+1/kappa)*gamma(1/sig2 - 1/kappa))/(sig2^(1+1/kappa)*gamma(1/sig2+1)))^(kappa)    
    names(theta) <- "theta"
    return(theta)    
    
  } else if(distCode == 4){ #generelized Gamma
    
    kappa <- distPara[1]
    gammaPara <- distPara[2]
    lambda <- exp(lgamma(kappa) - lgamma(kappa + 1 / gammaPara))
    names(lambda) <- "lambda"
    return(lambda)
    
  } else if(distCode == 5){ #generelized F
    
    kappa <- distPara[1]
    eta <- distPara[2]
    gammaPara <- distPara[3]
    
    lambda <- (lgamma(kappa) + lgamma(eta) 
               -(1/gammaPara)*log(eta)
               - lgamma(kappa + 1/gammaPara) - lgamma(eta - 1/gammaPara))
    lambda <- exp(lambda)
    names(lambda) <- "lambda"
    return(lambda)
    
  } else if(distCode == 6){ #q-Weibull
    
    a <- distPara[1]
    q <- distPara[2]
    
    b <- lgamma(1/(q-1)) - lgamma(1/a) - lgamma(1/(q-1) - 1/a - 1)
    b <- exp(b) * a * (q-1)^( (1+a) / a ) / (2 - q) 
    
    names(b) <- "b"
    return(b)
    
  } else if(distCode == 7){ #mixed q-Weibull and exponential
    
    p <- distPara[1]
    a <- distPara[2]
    q <- distPara[3]
    lambda <- distPara[4]
    
    b <- lgamma(1/(q-1)) - lgamma(1/a) - lgamma(1/(q-1) - 1/a - 1)
    b <- exp(b) * a * (q-1)^( (1+a) / a ) / (2 - q) 
    b <- b * (1 - (1 - p) * lambda) / p   
    
    names(b) <- "b"
    return(b)
    
  } else if(distCode == 8){ #mixed q-Weibull and Weibull
    
    p <- distPara[1]
    a <- distPara[2]
    q <- distPara[3]
    theta <- distPara[4]
    gamma <- distPara[5]
    
    b <- lgamma(1/(q-1)) - lgamma(1/a) - lgamma(1/(q-1) - 1/a - 1)
    b <- exp(b) * a * (q-1)^( (1+a) / a ) / (2 - q) 
    b <- b * (1 - (1 - p) * theta^(-1/gamma) * gamma(1/gamma + 1)) / p   
    
    names(b) <- "b"
    return(b)
    
  } else if(distCode == 9){ #"mixinvgauss" finite inverse Gaussian mixature
    
    return(NULL)
    
  } else if(distCode == 10){ #"birnbaum-saunders"  
    
    return(NULL)    
    
  } else stop("the provided distribution does not exist")    
}


.getCoef <- function(para, model = c("ACD","LACD1","LACD2","AMACD","SNIACD", "LSNIACD"), dist = c("exponential","weibull","burr"), 
                    hessian, order, bootError = NULL, bootCorr = NULL, bootMean = NULL, robustCorr = NULL, 
                    robustSE = NULL, fixedParam = NULL, fixedParamPos = NULL){
  
  #the standard error of the parameters, estimated from the numerical hessian of the log likelihood function:
  se <- sqrt(diag(solve(hessian)))
  
  #combines the para and fixedParam into the full para vector if there are fixed parameters:
  if(length(fixedParamPos) != 0){    
    para <- .returnfixedPara(freePara = para, fixedParam = fixedParam, fixedParamPos = fixedParamPos)
    se <- .returnfixedSE(freeSE = se, fixedParamPos = fixedParamPos)
  }
  
  
  comment <- NULL  
  if(model %in% c("ACD")){
    conDurPara <- para[1]
    paraNames <- "omega"
    NconDurPara <- 1
    if(order[1] != 0){
      for(j in 1:order[1]){
        conDurPara <- c(conDurPara, para[j+1])
        paraNames <- c(paraNames, paste("alpha", j, sep = ""))
        NconDurPara <- NconDurPara + 1
      } 
    }
    if(order[2] != 0){
      for(j in 1:order[2]){
        conDurPara <- c(conDurPara, para[j+1+order[1]])
        paraNames <- c(paraNames, paste("beta", j, sep = ""))
        NconDurPara <- NconDurPara + 1
      }
    }
    names(conDurPara) <-  paraNames
    pval <- 2*(1-stats::pnorm(abs(para/se)))[1:NconDurPara]
  } else if(model %in% c("LACD1", "LACD2")){
    conDurPara <- para[1]
    paraNames <- "omega"
    NconDurPara <- 1
    for(j in 1:order[1]){
      conDurPara <- c(conDurPara, para[j+1])
      paraNames <- c(paraNames, paste("alpha", j, sep = ""))
      NconDurPara <- NconDurPara + 1
    } 
    for(j in 1:order[2]){
      conDurPara <- c(conDurPara, para[j+1+order[1]])
      paraNames <- c(paraNames, paste("beta", j, sep = ""))
      NconDurPara <- NconDurPara + 1
    }
    names(conDurPara) <-  paraNames
    pval = 2*(1-stats::pnorm(abs(para/se)))[1:NconDurPara]
  } else if(model %in% c("AMACD")){
    conDurPara <- para[1]
    paraNames <- "omega"
    NconDurPara <- 1
    for(j in 1:order[1]){
      conDurPara <- c(conDurPara, para[j+1])
      paraNames <- c(paraNames, paste("alpha", j, sep = ""))
      NconDurPara <- NconDurPara + 1
    } 
    for(j in 1:order[2]){
      conDurPara <- c(conDurPara, para[j+1+order[1]])
      paraNames <- c(paraNames, paste("nu", j, sep = ""))
      NconDurPara <- NconDurPara + 1
    } 
    for(j in 1:order[3]){
      conDurPara <- c(conDurPara, para[j+1+order[1]+order[2]])
      paraNames <- c(paraNames, paste("beta", j, sep = ""))
      NconDurPara <- NconDurPara + 1
    }
    names(conDurPara) <-  paraNames
    pval = 2*(1-(stats::pnorm(abs(para/se))))[1:NconDurPara]
  } else if(model %in% c("ABACD")){
    conDurPara <- para[1]
    paraNames <- "omega"
    NconDurPara <- 1
    for(j in 1:order[1]){
      conDurPara <- c(conDurPara, para[j+1])
      paraNames <- c(paraNames, paste("alpha", j, sep = ""))
      NconDurPara <- NconDurPara + 1
    } 
    for(j in 1:order[1]){
      conDurPara <- c(conDurPara, para[j+1+order[1]])
      paraNames <- c(paraNames, paste("c", j, sep = ""))
      NconDurPara <- NconDurPara + 1
    } 
    for(j in 1:order[2]){
      conDurPara <- c(conDurPara, para[j+1+2*order[1]])
      paraNames <- c(paraNames, paste("beta", j, sep = ""))
      NconDurPara <- NconDurPara + 1
    }
    conDurPara <- c(conDurPara, para[2+2*order[1]+order[2]], para[3+2*order[1]+order[2]], para[4+2*order[1]+order[2]])
    paraNames <- c(paraNames, "nu", "delta1", "delta2")
    NconDurPara <- NconDurPara + 3
    
    names(conDurPara) <-  paraNames
    pval = 2*(1-(stats::pnorm(abs(para/se))))[1:NconDurPara]
  } else if(model %in% c("BACD")){
    conDurPara <- para[1]
    paraNames <- "omega"
    NconDurPara <- 1
    for(j in 1:order[1]){
      conDurPara <- c(conDurPara, para[j+1])
      paraNames <- c(paraNames, paste("alpha", j, sep = ""))
      NconDurPara <- NconDurPara + 1
    } 
    for(j in 1:order[2]){
      conDurPara <- c(conDurPara, para[j+1+order[1]])
      paraNames <- c(paraNames, paste("beta", j, sep = ""))
      NconDurPara <- NconDurPara + 1
    }
    conDurPara <- c(conDurPara, para[2+order[1]+order[2]], para[3+order[1]+order[2]])
    paraNames <- c(paraNames,"delta1", "delta2")
    NconDurPara <- NconDurPara + 2
    
    names(conDurPara) <-  paraNames
    pval = 2*(1-(stats::pnorm(abs(para/se))))[1:NconDurPara]
  } else if(model %in% c("SNIACD", "LSNIACD")){
    conDurPara <- para[1]
    paraNames <- "omega"
    NconDurPara <- 1
    for(j in 1:(order[3] + 1)){
      conDurPara <- c(conDurPara, para[j+1])
      paraNames <- c(paraNames, paste("c", j-1, sep = ""))
      NconDurPara <- NconDurPara + 1
    } 
    if(order[1] > 1)
      for(j in 1:(order[1] - 1)){
        conDurPara <- c(conDurPara, para[length(conDurPara) + 1])
        paraNames <- c(paraNames, paste("alpha", j, sep = ""))
        NconDurPara <- NconDurPara + 1
      }
    for(j in 1:order[2]){
      conDurPara <- c(conDurPara, para[length(conDurPara) + 1])
      paraNames <- c(paraNames, paste("beta", j, sep = ""))
      NconDurPara <- NconDurPara + 1
    }
    
    names(conDurPara) <-  paraNames
    pval = 2*(1-(stats::pnorm(abs(para/se))))[1:NconDurPara]
  } else stop("model not supported")
  
  if(dist == "exponential"){
    distPara <- NULL
    pval <- c(pval, NULL)
    comment <- c(comment, NULL)
  } else if(dist == "weibull") {
    distPara <- para[NconDurPara+1]
    paraNames <- c(paraNames, "gamma")
    names(distPara) <-  paraNames[NconDurPara+1]
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)]-1)/se[length(para)]))))
    comment <- c(comment, "The p-value for the distribution parameter gamma is from the 2-tailed test H0: gamma = 1.")
  } else if(dist == "burr") {
    distPara <- para[(NconDurPara+1):length(para)]
    paraNames <- c(paraNames, "kappa", "sigma2")
    names(distPara) <-  paraNames[(NconDurPara+1):length(para)]
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-1]-1)/se[length(para)-1]))))
    pval <- c(pval, (1-stats::pnorm((para[length(para)])/se[length(para)])))
    comment <- c(comment, "The p-value for the distribution parameter kappa is from the 2-tailed test H0: kappa = 1, and for sigma2 it is from the one sided test H0: sigma2 = 0 (or rather approching zero). If the two H0s are true, the Burr distribution reduces to the exponential distribution")
  } else if(dist == "gengamma") {
    distPara <- para[(NconDurPara+1):length(para)]
    paraNames <- c(paraNames, "kappa", "gamma")
    names(distPara) <-  paraNames[(NconDurPara+1):length(para)]
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-1]-1)/se[length(para)-1]))))
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)]-1)/se[length(para)]))))
    comment <- c(comment, "For the distribution parameters the null hypothesis is such that the parameter = 1 (2-sided). If the null is true, the generelized gamma distribution reduces to the exponential distribution")
  } else if(dist == "genf") {
    distPara <- para[(NconDurPara+1):length(para)]
    paraNames <- c(paraNames, "kappa", "eta", "gamma")
    names(distPara) <-  paraNames[(NconDurPara+1):length(para)]
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-2]-1)/se[length(para)-2]))))
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-1]-1)/se[length(para)-1]))))
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)]-1)/se[length(para)]))))
    comment <- c(comment, "The p-value for the distribution parameters are from the 2-tailed tests H0: distributionParameter = 1")
  } else if(dist == "qweibull") {
    distPara <- para[(NconDurPara+1):length(para)]
    paraNames <- c(paraNames, "a", "q")
    names(distPara) <-  paraNames[(NconDurPara+1):length(para)]
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-1]-1)/se[length(para)-1]))))
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)]-1)/se[length(para)]))))
    comment <- c(comment, "The p-value for the distribution parameters are from the 2-tailed tests H0: distributionParameter = 1")
  } else if(dist == "mixqwe") {
    distPara <- para[(NconDurPara+1):length(para)]
    paraNames <- c(paraNames, "p","a", "q", "lambda")
    names(distPara) <-  paraNames[(NconDurPara+1):length(para)]
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-3]-.5)/se[length(para)-3]))))
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-2]-1)/se[length(para)-2]))))
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-1]-1)/se[length(para)-1]))))
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-0]-1)/se[length(para)-0]))))
    comment <- c(comment, "The p-value for p is from the 2-tailed tests H0: p = .5, the rest of the distribution parameters are from H0: para = 1")
  } else if(dist == "mixqww") {
    distPara <- para[(NconDurPara+1):length(para)]
    paraNames <- c(paraNames, "p","a", "q", "theta", "gamma")
    names(distPara) <-  paraNames[(NconDurPara+1):length(para)]
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-4]-.5)/se[length(para)-4]))))
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-3]-1)/se[length(para)-3]))))
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-2]-1)/se[length(para)-2]))))
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-1]-1)/se[length(para)-1]))))
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-0]-1)/se[length(para)-0]))))
    comment <- c(comment, "The p-value for p is from the 2-tailed tests H0: p = .5, the rest of the distribution parameters are from H0: para = 1")
  } else if(dist == "mixinvgauss") {
    distPara <- para[(NconDurPara+1):length(para)]
    paraNames <- c(paraNames, "theta","lambda", "gamma")
    names(distPara) <-  paraNames[(NconDurPara+1):length(para)]
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-2]-0)/se[length(para)-2]))))
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-1]-0)/se[length(para)-1]))))
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-0]-0)/se[length(para)-0]))))
    comment <- c(comment, "The p-values for the distribution parameters are from the 2-tailed tests H0: para = 0")
  } else if(dist == "birnbaum-saunders") {
    distPara <- para[(NconDurPara+1):length(para)]
    paraNames <- c(paraNames, "kappa")
    pval <- c(pval, 2*(1-stats::pnorm(abs((para[length(para)-0]-1)/se[length(para)-0]))))
    comment <- c(comment, "The p-values for the distribution parameters are from the 2-tailed tests H0: kappa = 1")
  }
  
  if(length(fixedParamPos) != 0){
    paraNames <- ifelse(fixedParamPos, paste(paraNames, "(fixed)", sep = " "), paraNames)
  }
  
  parameterInference <- data.frame(Parameters = paraNames,
                                   Coef = para,
                                   SE = se,
                                   PV = round(pval, digits = 3),
                                   row.names = 1)
  
  #if bootstrapp where aviable: names the rows and columns of the correlation matrix and adds the mean and standard errors:
  if(length(bootError) != 0){
    parameterInference <- cbind(parameterInference, BootMean = bootMean, BootSE = bootError)
    dimnames(bootCorr) <- list(paraNames, paraNames)
  }  
  
  #names the rows and columns of the robust correlation matrix (if available):
  if(length(robustCorr) != 0){
    parameterInference <- cbind(parameterInference, robustSE = robustSE)
    dimnames(robustCorr) <- list(paraNames, paraNames)
  }   
  
  return(list(MPar = conDurPara, DPar = distPara, Inference = parameterInference, comment = comment, paraNames = paraNames, bootCorr = bootCorr, robustCorr = robustCorr))
}


.getdmudtheta_ACD <- function(param, x, order, mean = mean(x), newDay = c(0)){
  
  if(length(newDay) == 1 & newDay[1] == 0){
    NnewDays = 0;
  }else{
    NnewDays = length(newDay)
  }
  
  temp<-.C("getdmudtheta_ACD",
           as.double(x),
           as.integer(length(x)),
           as.double(param[1:(1+order[1]+order[2])]),
           as.integer(order),
           as.double(mean),
           as.double(numeric(length(x))),
           as.double(numeric(length(x))),
           as.integer(newDay),
           as.integer(NnewDays),
           as.double(numeric(length(x) * (1+order[1]+order[2]))), PACKAGE = "ACDm")
  
  .getdmudtheta_ACD <- matrix(temp[[10]], nrow = length(x), ncol = (1+order[1]+order[2]))    
}