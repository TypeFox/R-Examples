prepareDataSP <- function(time, event, marker){
  
  #probably put some checks here
  
  outData <- as.data.frame(cbind(time, event, marker, 1)) # add vector of 1s for sample weights used later

  names(outData) <- c("xi", "di", "Y", "wi")
  completeCases <- complete.cases(outData)
  if(sum(completeCases) < nrow(outData)) warning(paste("NA's present, only complete cases will be used. New sample size is:", sum(completeCases)) )
  outData <- outData[complete.cases(outData), ]
  
  outData
  
}


prepareDataNP <- function(time, event, marker){
  N = length(time)
  
  #calculate censoring weights
  psi = rep(0, N)
  for( eventID in unique(event)){
     psi[is.element(event, eventID)] <- sum(is.element(event, eventID))/N
  }
  
  outData <- data.frame(cbind(time, event, marker, 1, 0, event+1, psi, 1))  
  names(outData) = c("xi","di","yi","vi","zi","si","psi", "wi")
  
  completeCases <- complete.cases(outData)
  if(sum(completeCases) < nrow(outData)) warning(paste("NA's present, only complete cases will be used. New sample size is:", sum(completeCases)) )
  outData <- outData[complete.cases(outData), ]

  outData
  
}

processRawOutput <- function(myests, CImethod, alpha){


#calculate confidence intervals
if(substr(CImethod, 1, 4)=="stan"){
  
  myests$CIbounds = data.frame(rbind(upper = myests$estimates - qnorm(alpha/2)*myests$se, 
                                     lower = myests$estimates - qnorm(1-alpha/2)*myests$se))
  names(myests$CIbounds) = names(myests$estimates)
  
}else{
  #logit transform everything but the coef

  mynames = names(myests$estimates)
  myests$CIbounds = data.frame(rbind(upper = expit(logit(myests$estimates[,mynames!="coef"]) - 
                                                     qnorm(alpha/2)*myests$se[,mynames!="coef"]/(myests$estimates[,mynames!="coef"]*(1-myests$estimates[,mynames!="coef"]))), 
                                     lower = expit(logit(myests$estimates[,mynames!="coef"]) - 
                                                     qnorm(1-alpha/2)*myests$se[,mynames!="coef"]/(myests$estimates[,mynames!="coef"]*(1-myests$estimates[,mynames!="coef"])))))
  
  if(is.element("coef", mynames)){
  myests$CIbounds = cbind(data.frame(rbind(upper = myests$estimates[,mynames=="coef"] - qnorm(alpha/2)*myests$se[,mynames=="coef"], 
                                           lower = myests$estimates[,mynames=="coef"] - qnorm(1-alpha/2)*myests$se[,mynames=="coef"])), myests$CIbounds)
  }
  names(myests$CIbounds) = names(myests$estimates)
  
}

myests$model.fit <- myests$fit; 

myests
}