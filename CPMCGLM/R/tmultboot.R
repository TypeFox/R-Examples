####test du score
tmultboot <- function(formula,data,Y1,Z,codage,family1,family=family,link1){

#Definition du seuil de tolérance
tol <- 10E-40

#Definition de la formule sous H0
form<-formula(paste("Y1",formula))

#Calcul des paramètres du modèle sous H0
model.h0 <- glm(formula=form,data=data,family=family1(link=link1))

pi1 <- model.h0$fitted
V<-switch(family,
	"binomial"=diag(pi1*(1-pi1),length(Y1)),
	"gaussian"=diag(var(Y1),length(Y1)),
	"poisson"=diag(model.h0$weight,length(Y1)),
	)

ZZ<-J<-t<-NULL

num<-NULL

#Calcul de la statistique de test pour chaque échantillon
for (j in 1:ncol(codage))
{
ZZ[[j]] <- cbind(rep(1,nrow(data)),Z,codage[,j])
num[[j]] <- t(ZZ[[j]])%*%(Y1-pi1)
temp <- t(ZZ[[j]])%*%V%*%ZZ[[j]]
		
		if((det(temp)==0)||(1/det(temp)<tol)){
		t[j] <- NA
		}else{
			J[[j]] <-  solve(temp,tol=tol)
			t[j] <- t(num[[j]])%*%J[[j]]%*%num[[j]]
		}
		
}
res<-t
}
