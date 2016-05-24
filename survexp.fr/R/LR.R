LR <-
function(futime,status,age,sex,entry_date,ratetable=survexp.fr){
  data=na.omit(data.frame(futime,status,age,sex,entry_date))
  data$risk=-log(survexp(futime~1,data=data,rmap=list(year=entry_date,age=age,sex=sex),cohort=F,ratetable=ratetable,conditional=TRUE))
  O=sum(data$status)
  E=sum(data$risk)
  LR=((O-E)^2)/E
  p.value=1-pchisq(LR,1)
  return(list(O=O,E=E,LR=LR,p.value=p.value))
}
