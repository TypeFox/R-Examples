TwoSampleSurvival.Conditional <-
function(alpha,beta,lam1,lam2,eta1,eta2,k,ttotal,taccrual,g1,g2){

gamma1=g1
gamma2=g2

variance1<-lam1^2*(lam1/(lam1+eta1)+lam1*gamma1*exp(-(lam1+eta1)*ttotal)*(1-exp((lam1+eta1-gamma1)*taccrual))/((lam1+eta1-gamma1)*(1-exp(-gamma1*taccrual))*(lam1+eta1)))^-1
variance2<-lam2^2*(lam2/(lam2+eta1)+lam2*gamma2*exp(-(lam2+eta2)*ttotal)*(1-exp((lam2+eta2-gamma2)*taccrual))/((lam2+eta2-gamma2)*(1-exp(-gamma2*taccrual))*(lam2+eta2)))^-1

print(variance1)
print(variance2)

lam<-(k*lam1+lam2)/(k+1)
eta<-(k*eta1+eta2)/(k+1)
gamma<-(k*gamma1+gamma2)/(k+1)

variance<-lam^2*(lam/(lam+eta)+lam*gamma*exp(-(lam+eta)*ttotal)*(1-exp((lam+eta-gamma)*taccrual))/((lam+eta-gamma)*(1-exp(-gamma*taccrual))*(lam+eta)))^-1

n2=1/(lam2-lam1)^2*(qnorm(1-alpha/2)*variance*(1/k+1)+qnorm(1-beta)*(variance1/k+variance2)^0.5)^2

n2
}
