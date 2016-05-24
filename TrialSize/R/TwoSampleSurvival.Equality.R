TwoSampleSurvival.Equality <-
function(alpha,beta,lam1,lam2,k,ttotal,taccrual,gamma){

variance1<-lam1^2*(1+gamma*exp(-lam1*ttotal)*(1-exp((lam1-gamma)*taccrual))/((lam1-gamma)*(1-exp(-gamma*taccrual))))^-1
variance2<-lam2^2*(1+gamma*exp(-lam2*ttotal)*(1-exp((lam2-gamma)*taccrual))/((lam2-gamma)*(1-exp(-gamma*taccrual))))^-1

print(variance1)
print(variance2)

n2=(qnorm(1-alpha/2)+qnorm(1-beta))^2*(variance1/k+variance2)/(lam1-lam2)^2
n2
}
