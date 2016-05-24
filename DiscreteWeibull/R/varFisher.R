varFisher <-
function(x,zero=FALSE)
{
k<-x
par<-estdweibull(x,"ML",zero=zero)
q<-par[1]
beta<-par[2]
L<-ifelse(zero,expression(log(q^(k)^beta-q^(k+1)^beta)),expression(log(q^(k-1)^beta-q^k^beta)))
DLq<-D(L,"q")
Dq<-function(q,beta){if (k!=as.numeric(!zero)) eval(DLq) else -1/(1-q)}
DLbeta<-D(L,"beta")
Dbeta<-function(q,beta){if (k!=as.numeric(!zero)) eval(DLbeta) else 0}
DLqbeta<-D(DLq,"beta")
Dqbeta<-function(q,beta){if (k!=as.numeric(!zero)) eval(DLqbeta) else 0}
Dbetaq<-Dqbeta
DLqq<-D(DLq,"q")
Dqq<-function(q,beta){if (k!=as.numeric(!zero)) eval(DLqq) else -1/(1-q)^2}
DLbetabeta<-D(DLbeta,"beta")
Dbetabeta<-function(q,beta){if (k!=as.numeric(!zero)) eval(DLbetabeta) else 0}
I11<-0
I22<-0
I12<-0
for(t in 1:length(x))
{
k<-x[t]
I22<-I22-Dbetabeta(q,beta)
I11<-I11-Dqq(q,beta)
I12<-I12-Dqbeta(q,beta)
}
# Fisher Information matrix
IF<-matrix(c(I11,I12,I12,I22),2,2)/length(x)
# inverse of Fisher Information matrix
IFI<-qr.solve(IF, tol = 1e-12)
# asymptotic variance/covariance matrix
list(FisherInfMatrix=IF,InvFisherInfMatrix=IFI)
}

