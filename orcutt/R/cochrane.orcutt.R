cochrane.orcutt <-
function(reg){
  X <- model.matrix(reg)
  Y <- model.response(model.frame(reg))
  n<-length(Y)
  e<-reg$residuals  
  e2<-e[-1]
  e3<-e[-n]
  regP<-lm(e2~e3-1)
  rho<-summary(regP)$coeff[1]
  rho2<-c(rho)
  XB<-X[-1,]-rho*X[-n,]
  YB<-Y[-1]-rho*Y[-n]
  regCO<-lm(YB~XB-1)
  ypCO<-regCO$coeff[1]+as.matrix(X[,-1])%*%regCO$coeff[-1]  
  e1<-ypCO-Y
  e2<-e1[-1]
  e3<-e1[-n]
  regP<-lm(e2~e3-1) 
  rho<-summary(regP)$coeff[1]
  rho2[2]<-rho
  i<-2
  while (round(rho2[i-1],8)!=round(rho2[i],8)){
    XB<-X[-1,]-rho*X[-n,]
    YB<-Y[-1]-rho*Y[-n]
    regCO<-lm(YB~XB-1)
    ypCO<-regCO$coeff[1]+as.matrix(X[,-1])%*%regCO$coeff[-1]  
    e1<-ypCO-Y
    e2<-e1[-1]
    e3<-e1[-n]
    regP<-lm(e2~e3-1)
    rho<-summary(regP)$coeff[1]  
    i<-i+1
    rho2[i]<-rho
    }
  result<-list()
  result$Cochrane.Orcutt<-summary(regCO)
  result$rho<-rho
  result$number.interaction<-i-1
  return(result)
}
