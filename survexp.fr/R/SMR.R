SMR <-
function(futime,status,age,sex,entry_date,ratetable=survexp.fr,alpha=0.05){
  data=na.omit(data.frame(futime,status,age,sex,entry_date))
  data$risk=-log(survexp(futime~1,data=data,rmap=list(year=entry_date,age=age,sex=sex),cohort=F,ratetable=ratetable,conditional=TRUE))

  O=sum(data$status)
  E=sum(data$risk)
  
  # calcul classique
  SMR=O/E
SMR.lo=O/E*(1-1/9/O-qnorm(1-alpha/2)/3/sqrt(O))^3
SMR.up=(O+1)/E*(1-1/9/(O+1)+qnorm(1-alpha/2)/3/sqrt(O+1))^3
  if (E>=10){
  chisq=((abs(O-E)-0.5)^2)/E
  } else{
  Mstar=ifelse(O>=E,O,O+1)
  chisq=9*Mstar*(1-(1/(9*Mstar))-((E/Mstar)^(1/3)))^2
}
p.value=1-pchisq(chisq,1)
  SMR.classic=list(SMR=SMR,SMR.lo=SMR.lo,SMR.up=SMR.up,p.value=p.value)
    
  # modèle de Poisson
  fit=glm(O~1+offset(log(E)),family=poisson)    
  coefs=summary(fit)$coefficients
  SMR=exp(coefs[1,1])
  SMR.lo=exp(coefs[1,1]-qnorm(1-alpha/2)*coefs[1,2])
  SMR.up=exp(coefs[1,1]+qnorm(1-alpha/2)*coefs[1,2])  
  p.value=coefs[1,4]
  SMR.poisson=list(SMR=SMR,SMR.lo=SMR.lo,SMR.up=SMR.up,p.value=p.value)
  
  return(list(O=O,E=E,SMR.classic=SMR.classic,SMR.poisson=SMR.poisson))
}
