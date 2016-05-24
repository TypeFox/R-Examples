minxent.single<-function (q, G, eta, lambda) 
{
fk2<-function (lambda, q, G, eta) 
{
G1<-G[-1,]
lambda0<-log(sum(q*exp(-lambda*G1)))    
(q * exp(-lambda0)*exp(-lambda * G1))%*%t(G) - eta
}
fkss2<-function (p, q, G, eta) 
{
    sum(fk2(p, q, G, eta)^2)
}

out<-optimize(fkss2,c(-2,2),q=q,G=G,eta=eta)
lambda<-out$minimum
lambda0<-log(sum(q*exp(-lambda*G[-1,])))    
pi_solve<-(q * exp(-lambda0)*exp(-lambda * G[-1,]))
list(Langrangians= c(lambda0,lambda) , Estimates=pi_solve)

}
