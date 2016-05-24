ControlCharts <-
function(data,xlabel,ylabel,CMP=TRUE,P=TRUE, CMPProb=TRUE){

  
  data.cmp <- ComputeLambdaAndNuHat.shiftCMPest(data)
  
  data.lambda <- data.cmp$par[1]
  data.nu <- data.cmp$par[2]
  data.cmpmean <- com.mean(data.lambda, data.nu) + min(data)
  data.cmpvar <- com.var(data.lambda, data.nu)
  data.cmpstdev <- sqrt(data.cmpvar)
  
  uppercontrolCMP <- data.cmpmean + 3*data.cmpstdev
  lowercontrolCMP <- data.cmpmean - 3*data.cmpstdev
  lowercontrolCMPPlot <- data.cmpmean - 3*data.cmpstdev
  
  if (lowercontrolCMP < 0){lowercontrolCMPPlot = 0}
  
  
  data.mean <- mean(data)
  data.stdev <- sqrt(data.mean)
  
  uppercontrolPoisson <- data.mean + 3*data.stdev
  lowercontrolPoisson <- data.mean - 3*data.stdev
  lowercontrolPoissonPlot <- data.mean - 3*data.stdev
  
  if (lowercontrolPoisson < 0){lowercontrolPoissonPlot = 0}
  
  TrueMax <- max(data)
  
  CumulativeProbMax <- sum(exp(com.log.density(0:max(data),data.lambda,data.nu)))
  
  if(CumulativeProbMax < 0.99865){
    
    j <- 1
    
    CumulativeProbMax_Cur <- CumulativeProbMax
  
  while((CumulativeProbMax_Cur) < 0.99865){
    
    CumulativeProbMax_Cur <- sum(exp(com.log.density((0:(max(data)+j)),data.lambda,data.nu)))
    
    TrueMax <- TrueMax + 1
    
    j <- j + 1
    
  }
  
  TrueMax <- TrueMax + 1
  
  }
  
  
  if(CumulativeProbMax >= 0.99865){
    
    k <- 1
    
    CumulativeProbMax_Cur <- CumulativeProbMax
    
    while(CumulativeProbMax_Cur > 0.99865){
      
      CumulativeProbMax_Cur <- sum(exp(com.log.density(0:(max(data-k)),data.lambda,data.nu)))
      
      TrueMax <- TrueMax - 1
      
      k <- k + 1
      
    }
  
    TrueMax <- TrueMax + 2
    
  }
  
  TrueMin <- min(data)
  
  CumulativeProbMin <- sum(exp(com.log.density(0:min(data),data.lambda,data.nu)))
  
  if(CumulativeProbMin > 0.00135 & min(data) > 0){
  
    i <- 1
    
    CumulativeProbMin_Cur <- CumulativeProbMin
  
      while (CumulativeProbMin_Cur > 0.00135){
    
      CumulativeProbMin_Cur <- sum(exp(com.log.density((0:(min(data)-i)), data.lambda, data.nu)))
    
      TrueMin <- TrueMin - 1
    
      i <- i + 1
    
        if (TrueMin==0){break}
    
  }
  
  TrueMin <- TrueMin
  
  }
  
  if(CumulativeProbMin <= 0.00135){
    
    l <- 1
    
    CumulativeProbMin_Cur <- CumulativeProbMin
    
      while(CumulativeProbMin_Cur < 0.00135){
        
      CumulativeProbMin_Cur <- sum(exp(com.log.density(0:(min(data)+l), data.lambda,data.nu)))
      
      TrueMin <- TrueMin + 1
      
      l <- l+1
        
      }
  
      TrueMin <- TrueMin - 1
  
  }
  
  plot(data, type="l", ylim=c(min(lowercontrolCMPPlot, min(data),lowercontrolPoissonPlot, TrueMin), max(uppercontrolCMP, max(data), uppercontrolPoisson, TrueMax)),xlab=xlabel, ylab=ylabel, sub="Red = CMP Shewhart ; Blue = Poisson Shewhart ; Green = CMP Probability Limits")
  
  if (CMP==TRUE){
  abline(h=uppercontrolCMP, col=2)
  abline(h=lowercontrolCMPPlot, col=2)
  abline(h=data.cmpmean, col=2)}
  
  if (P==TRUE){
  abline(h=uppercontrolPoisson, col=4)
  abline(h=lowercontrolPoissonPlot, col=4)
  abline(h=data.mean, col=4)}
  
  if (CMPProb==TRUE){
  abline(h=TrueMax, col=3)
  abline(h=TrueMin, col=3)}

  UpperOOCPoints <- c()
  LowerOOCPoints <- c()
  
  for (i in 1:length(data)){
    if (data[i] > uppercontrolPoisson || data[i] > uppercontrolCMP){UpperOOCPoints <- c(UpperOOCPoints,which(data==data[i]))}
    
    if (data[i] < lowercontrolPoisson || data[i] < lowercontrolCMP){LowerOOCPoints <- c(LowerOOCPoints,which(data==data[i]))}
  }
  
  if (length(UpperOOCPoints)==0){UpperOOCPoints <- c("No Upper Out of Control Points")}
  
  if (length(LowerOOCPoints)==0){LowerOOCPoints <- c("No Lower Out of Control Points")}
  
 
  list.cmp <- list(data.cmp[[1]], c(data.cmpmean, data.cmpstdev) ,c(uppercontrolCMP, lowercontrolCMP), c(data.mean, data.stdev), c(uppercontrolPoisson, lowercontrolPoisson), UpperOOCPoints, LowerOOCPoints, c(TrueMax,TrueMin))
  names(list.cmp) <- c("CMP Lambda Hat and Nu Hat","CMP Mean and Standard Deviation", "CMP Shewhart Upper and Lower Bounds", "Poisson Mean and Standard Deviation", "Poisson Shewhart Upper and Lower Bounds", "Upper Out of Control Observations", "Lower Out of Control Observations", "CMP Probability Limits")
  
  return(list.cmp)

}
