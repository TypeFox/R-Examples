

get.est <-
function(x  ,  data  ,  formula ){

  s1<-Surv(time = as.numeric(x) , event = data$right != Inf )
  formula1<-formula(paste('s1' , deparse(formula)))
  cph1<-coxph(formula1 , data = data)
  vc1<-vcov(cph1)
  coef1<-cph1$coefficients
  return(list( sigma = vc1 , beta = coef1 ))

}
