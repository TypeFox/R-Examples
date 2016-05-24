test.score <- function(formula,data,codage,Z,Y,family1,family,link1,ind.ajust){
	tol=10E-40

	model.h0 <- glm(formula,family=family1(link=link1),data=data)

	pi2 <- model.h0$fitted
	
	V<-switch(family,
	"binomial"=diag(pi2*(1-pi2),length(Y)),
	"gaussian"=diag(var(Y),length(Y)),
	"poisson"=diag(model.h0$weight,length(Y)),
	)

	ZZ<-J<-t<-NULL
	num<-NULL
	for (i in 1:ncol(codage)){
		ZZ[[i]] <- cbind(rep(1,nrow(data)),Z,codage[,i])
		num[[i]] <- t(ZZ[[i]])%*%(Y-pi2)
		temp <- t(ZZ[[i]])%*%V%*%ZZ[[i]]
		if((det(temp)==0)||(1/det(temp)<tol)){
			t[i] <- NA
		}else{
			J[[i]] <-  solve(temp,tol=tol)
			t[i] <- t(num[[i]])%*%J[[i]]%*%num[[i]]
		}
	}
	res<-t
	
}

