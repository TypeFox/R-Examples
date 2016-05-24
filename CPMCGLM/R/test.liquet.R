test.liquet<-function(formula,data,codage,Z,Y,family,link1,family1)
{
 model.h0 <- glm(formula,family=family1(link=link1),data)
 tol <- 10E-40

#Pour la famille exponentielle

phi<-switch(family,
	"binomial"=1,
	"gaussian"=var(Y),
	"poisson"=1,
	)



 pi3 <- model.h0$fitted

V<-switch(family,
	"binomial"=diag(pi3*(1-pi3),length(Y)),
	"gaussian"=diag(var(Y),length(Y)),
	"poisson"=diag(model.h0$weight,length(Y)),
	)


 I<-diag(nrow(codage))
 ZZ<-cbind(rep(1,nrow(data)),Z)

 H<-V%*%ZZ%*%(solve(t(ZZ)%*%V%*%ZZ))%*%t(ZZ)
 #initialisation des param�tres
 num<-NULL
num1<-NULL
 J<-NULL
 t<-NULL
# t1<-NULL
 Var<-NULL
ZZ1<-NULL

for(i in 1:ncol(codage))
 {
 num[i]<- (1/phi)*t(codage[,i])%*%(Y-pi3)
 Var[i]<- (1/phi)^2*t(codage[,i])%*%(I-H)%*%V%*%codage[,i]
# t[i]<-num[i]/(sqrt(Var[i]))
		ZZ1[[i]] <- cbind(rep(1,nrow(data)),Z,codage[,i])
		num1[[i]] <- t(ZZ1[[i]])%*%(Y-pi3)
		temp <- t(ZZ1[[i]])%*%V%*%ZZ1[[i]]
		if((det(temp)==0)||(1/det(temp)<tol)){
			t[i] <- NA
		}else{
			J[[i]] <-  solve(temp,tol=tol)
			t[i] <- sqrt(t(num1[[i]])%*%J[[i]]%*%num1[[i]])
		}

 }

 Corr<-matrix(nrow=ncol(codage),ncol=ncol(codage))
 for (i in 1:ncol(codage))
 {
 for (j in 1:ncol(codage))
 {
 Corr[i,j]<-((1/phi)^2*codage[,i]%*%(I-H)%*%V%*%codage[,j])/(sqrt(Var[i])%*%sqrt(Var[j]))
 }
 }

pval.exact<-NULL
for(i in 1:ncol(codage))
 {
 pval.exact<-c(pval.exact,1-pmvnorm(lower=-rep(max(abs(t[1:i])),i),upper=rep(max(abs(t[1:i])),i),mean=rep(0,i),sigma=Corr[1:i,1:i]))
 }
 res <- pval.exact
 }
