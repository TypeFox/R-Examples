AER <-
function(futime,status,age,sex,entry_date,PY.stand=10000,ratetable=survexp.fr,alpha=0.05){
  data=na.omit(data.frame(futime,status,age,sex,entry_date,person_year=futime/365.241))
  data$risk=-log(survexp(futime~1,data=data,rmap=list(year=entry_date,age=age,sex=sex),cohort=F,ratetable=ratetable,conditional=TRUE))
  
  mypoisson=poisson()
  mypoisson$link="glm with Poisson error and non standard link function"
  mypoisson$linkfun=function(mu) log(mu-E)
  mypoisson$linkinv=function(eta) exp(eta)+E
  mypoisson$initialize=expression({
    if (any(y < 0)) stop("Negative values not allowed for the Poisson family")
    n <- rep.int(1, nobs)
    mustart <- pmax(y, linkinv(-1000)) + 0.1                                    # linkinv(-1000) returns expected
  })  

  E=sum(data$risk)
  O=sum(data$status)
  PY.observed=sum(data$person_year)
  
  fit=glm(O~1+offset(log(PY.observed)),family=mypoisson)
  coefs=summary(fit)$coefficients
  AER=exp(coefs[1,1]+log(PY.observed))*PY.stand/PY.observed                     # le +log(PY) et /PY s'annulent
  AER.lo=exp(coefs[1,1]-qnorm(1-alpha/2)*coefs[1,2]+log(PY.observed))*PY.stand/PY.observed
  AER.up=exp(coefs[1,1]+qnorm(1-alpha/2)*coefs[1,2]+log(PY.observed))*PY.stand/PY.observed
  return(list(AER=AER,PY.stand=PY.stand,AER.lo=AER.lo,AER.up=AER.up,p.value=coefs[1,4],O=O,E=E,PY.observed=PY.observed))
}
