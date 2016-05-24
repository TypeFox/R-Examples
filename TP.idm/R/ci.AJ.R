ci.AJ <-
function(s,t,level,ci.transformation="linear",dNs.id_tr,TP.AJs,cov.AJs,e.times.id_tr){
  
  if(s==0){
    p.transitions <- c("1 1", "1 2", "1 3")
    
    alpha<-qnorm(level+(1-level)/2)
    
    CI.tp<-array(0,dim=c(nrow(dNs.id_tr),4,length(p.transitions)),
                 dimnames=list(rows=e.times.id_tr,cols=c("probs","lower","upper","variance"),trans=p.transitions))
    
    # different transformations to built CI
    ci.transformation <- match.arg(ci.transformation, c("linear", "log", "cloglog", "log-log"))
    
    for (j in 1:length(p.transitions)) { ## loop through possible transitions
      
      idx <- unlist(strsplit(p.transitions[j], " "))
      CI.tp[ , 1, j] <- P <- TP.AJs[idx[1], idx[2] , ]
      CI.tp[ , 4, j] <- var <- cov.AJs[p.transitions[j], p.transitions[j], ]
      
      
      switch(ci.transformation[1],
             "linear" = {
               CI.tp[ , 2, j] <- P - alpha * sqrt(var)
               CI.tp[ , 3, j] <- P + alpha * sqrt(var)},
             "log" = {
               CI.tp[ , 2, j] <- exp(log(P) - alpha * sqrt(var) / P)
               CI.tp[ , 3, j] <- exp(log(P) + alpha * sqrt(var) / P)},
             "cloglog" = {
               CI.tp[ , 2, j] <- 1 - (1 - P)^(exp(alpha * (sqrt(var) /
                                                             ((1 - P) * log(1 - P)))))
               CI.tp[ , 3, j] <- 1 - (1 - P)^(exp(-alpha * (sqrt(var) /
                                                              ((1 - P) * log(1 - P)))))},
             "log-log" = {
               CI.tp[ , 2, j] <- P^(exp(-alpha * (sqrt(var) / (P * log(P)))))
               CI.tp[ , 3, j] <- P^(exp(alpha * (sqrt(var) / (P * log(P)))))})
      
      CI.tp[ , 2, j] <- pmax(CI.tp[ , 2, j], 0)
      CI.tp[ , 3, j] <- pmin(CI.tp[ , 3, j], 1)
      
    } ## end j loop
    
    CI.t<- matrix(0, nrow=length(p.transitions), ncol=4)
    colnames(CI.t) <- c("probs","lower","upper","variance")
    rownames(CI.t) <- p.transitions
    
    for(j in 1:length(p.transitions)){
      CI.t[j,]<-CI.tp[nrow(CI.tp[,,j]),,j] 
    }
    
    CI.tp <- round(CI.tp,7)
    CI.t <-round(CI.t,7)
  }else{
    p.transitions <- c("1 1", "1 2", "1 3", "2 2", "2 3")
    
    alpha<-qnorm(level+(1-level)/2)
    
    CI.tp<-array(0,dim=c(nrow(dNs.id_tr),4,length(p.transitions)),
                 dimnames=list(rows=e.times.id_tr,cols=c("probs","lower","upper","variance"),trans=p.transitions))
    
    # different transformations to built CI
    ci.transformation <- match.arg(ci.transformation, c("linear", "log", "cloglog", "log-log"))
    
    for (j in 1:length(p.transitions)) { ## loop through possible transitions
      
      idx <- unlist(strsplit(p.transitions[j], " "))
      CI.tp[ , 1, j] <- P <- TP.AJs[idx[1], idx[2] , ]
      CI.tp[ , 4, j] <- var <- cov.AJs[p.transitions[j], p.transitions[j], ]
      
      
      switch(ci.transformation[1],
             "linear" = {
               CI.tp[ , 2, j] <- P - alpha * sqrt(var)
               CI.tp[ , 3, j] <- P + alpha * sqrt(var)},
             "log" = {
               CI.tp[ , 2, j] <- exp(log(P) - alpha * sqrt(var) / P)
               CI.tp[ , 3, j] <- exp(log(P) + alpha * sqrt(var) / P)},
             "cloglog" = {
               CI.tp[ , 2, j] <- 1 - (1 - P)^(exp(alpha * (sqrt(var) /
                                                             ((1 - P) * log(1 - P)))))
               CI.tp[ , 3, j] <- 1 - (1 - P)^(exp(-alpha * (sqrt(var) /
                                                              ((1 - P) * log(1 - P)))))},
             "log-log" = {
               CI.tp[ , 2, j] <- P^(exp(-alpha * (sqrt(var) / (P * log(P)))))
               CI.tp[ , 3, j] <- P^(exp(alpha * (sqrt(var) / (P * log(P)))))})
      
      CI.tp[ , 2, j] <- pmax(CI.tp[ , 2, j], 0)
      CI.tp[ , 3, j] <- pmin(CI.tp[ , 3, j], 1)
      
    } ## end j loop
    
    CI.t<- matrix(0, nrow=length(p.transitions), ncol=4)
    colnames(CI.t) <- c("probs","lower","upper","variance")
    rownames(CI.t) <- p.transitions
    
    for(j in 1:length(p.transitions)){
      CI.t[j,]<-CI.tp[nrow(CI.tp[,,j]),,j] 
    }
    
    CI.tp <- round(CI.tp,7)
    CI.t <-round(CI.t,7)
  }
  
  return(list(CI=CI.t,all.CI=CI.tp))
  
}
