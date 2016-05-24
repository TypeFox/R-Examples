
####################################################################
#' Initial values for the variance components in Model 3
#'
#' This function is used in \code{\link[mme]{initial.values}} to calculate the initial values for the variance
#' components in the multinomial mixed model with two independent random effects
#' for each category of the response variable: one domain random effect and another correlated time and domain random effect (Model 3).
#'
#' @param beta.0 a list with the initial values for the fixed effects per category obtained from \code{\link[mme]{initial.values}}.
#' @param y matrix with the response variable obtained from \code{\link[mme]{data.mme}}. The rows are the domains and the columns are the categories of the response variable minus one.
#' @param Xk list of matrices with the auxiliary variables per category obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of domains.
#' @param M vector with the sample size of the areas.
#' @param u1 matrix with the values for the first random effect obtained from \code{\link[mme]{initial.values}}.
#' @param u2 matrix with the values for the second random effect obtained from \code{\link[mme]{initial.values}}.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#'\code{\link[mme]{wmatrix}}, \code{\link[mme]{prmu.time}},
#' \code{\link[mme]{Fbetaf.ct}}, \code{\link[mme]{phi.direct.ct}},
#' \code{\link[mme]{sPhikf.ct}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit3}}, \code{\link[mme]{msef.ct}},\code{\link[mme]{omega}},
#' \code{\link[mme]{mseb}}.
#' @return A list containing the following components.
#' \item{phi.0}{vector of the initial values for the variance components.}
#' \item{rho.0}{vector of the initial values for the correlation parameter.}
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Small area estimation of labour force indicators under a multinomial mixed model
#' with correlated time and area effects. Submitted for review.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' mod=3  #type of model
#' data(simdata3) #data
#' D=nrow(simdata3)
#' datar=data.mme(simdata3,k,pp,mod)
#' ###Fixed effects values
#' beta.new=list()
#' beta.new[[1]]=matrix(c( 1.3,-1),2,1)
#' beta.new[[2]]=matrix(c( -1.6,1),2,1)
#' ## Random effects values
#' u1.new=rep(0.01,((k-1)*datar$d))
#' dim(u1.new)=c(datar$d,k-1)
#' u2.new=rep(0.01,((k-1)*D))
#' dim(u2.new)=c(D,k-1)
#'
#' ## Initial variance components
#' phi=phi.mult.ct(beta.new,datar$y,datar$Xk,datar$n,u1.new,u2.new)

phi.mult.ct<-function(beta.0,y,Xk,M,u1,u2){
D=nrow(u2)
d=nrow(u1)
t=D/d
k=ncol(y)
theta.ol=matrix(0,D,(k-1))
prmul=prmu.time(M,Xk,beta.0,u1,u2)
theta=prmul[[3]]

for (i in 1:(k-1)){
	for (j in 1:D){
		if ((y[j,i]==0) |(y[j,k]==0)) {theta.ol[j,i]=log((y[j,i]+1)/(y[j,k]+1))} else {theta.ol[j,i]=log(y[j,i]/y[j,k])}}}
dif=(theta-theta.ol)^2
phi1.new=colSums(dif)/(2*(D-1))
dif=sqrt(dif)
a=rep(0,(k-1))
b=rep(0,(k-1))
for (kk in 1:(k-1)){
	i=1
	j=1
	while (i<(D-1)){
		a[kk]=a[[kk]]+dif[i,kk]*(dif[(i+1),kk]-dif[(i+2),kk])
		b[kk]=b[[kk]]+dif[i,kk]*(dif[(i),kk]-dif[(i+1),kk])
		if (j<(t-2)){
		i=i+1
		j=j+1}
		else
		{
			j=1
			i=i+3}
		}}

rho=a/b
if (max(abs(rho))>0.99){
	for ( i in 1:(k-1)){
		if (rho[i]< -1) {	rho[i]=-0.8}
		if (rho[i]> 1) {	rho[i]=0.8}}}
salida=list(phi.0=phi1.new,rho.0=rho)
return(salida)
}

####################################################################
#' Inverse of the Fisher information matrix of fixed and random effects in Model 3
#'
#' This function calculates the score vector S and the inverse of the Fisher information
#' matrix for the fixed (beta) and the random effects (u1, u2) in Model 3. This model has two independet sets of random effects. 
#' The first one contains independent random effects u1dk associated to each category and domain. The second set contains random effects
#' u2dkt associated to each category, domain and time period. Model 3 assumes that the u2dk are AR(1) correlated across time. 
#' \code{\link[mme]{modelfit3}} uses the output of this function to estimate the fixed and random effect by the PQL method.
#'
#' @param sigmap a list with the model variance-covariance matrices for each domain obtained from \code{\link[mme]{wmatrix}}.
#' @param X list of matrices with the auxiliary variables obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of categories of the response variable minus one.
#' @param Z design matrix  of random effects.
#' @param phi1 vector with the values of the variance components for the first random effects obtained from \code{\link[mme]{modelfit3}}.
#' @param phi2 vector with the values of the variance components for the second random effects obtained from \code{\link[mme]{modelfit3}}.
#' @param y matrix with the response variable, except the reference category. The rows are the domains and the columns are the categories of the response variable minus one.
#' @param mu matrix with the estimated mean of the response variable.
#' @param u1 matrix with the values of the first random effect obtained from \code{\link[mme]{modelfit3}}.
#' @param u2 matrix with the values of the second random effect obtained from \code{\link[mme]{modelfit3}}.
#' @param rho vector with the values of the correlation parameter obtained from \code{\link[mme]{modelfit3}}.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#' \code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult.ct}},
#' \code{\link[mme]{prmu.time}}, \code{\link[mme]{phi.direct.ct}},
#' \code{\link[mme]{sPhikf.ct}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit3}}, \code{\link[mme]{msef.ct}}, \code{\link[mme]{omega}},
#' \code{\link[mme]{mseb}}
#' @return A list containing the following components.
#' \item{F}{the inverse of the Fisher information matrix of (beta, u1, u2).}
#' \item{S}{(beta, u1, u2) score vectors}
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Small area estimation of labour force indicators under a multinomial mixed model
#' with correlated time and area effects. Submitted for review.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' mod=3 #type of model
#' data(simdata3)
#' datar=data.mme(simdata3,k,pp,mod)
#' initial=datar$initial
#' mean=prmu.time(datar$n,datar$Xk,initial$beta.0,initial$u1.0,initial$u2.0)
#' sigmap=wmatrix(datar$n,mean$estimated.probabilities) #variance-covariance
#'
#' ##The inverse of the Fisher information matrix and the score matrix
#' Fisher.beta=Fbetaf.ct(sigmap,datar$X,datar$Z,initial$phi1.0,initial$phi2.0,
#'          datar$y[,1:(k-1)],mean$mean,initial$u1.0,initial$u2.0,initial$rho.0)


Fbetaf.ct<- function(sigmap,X,Z,phi1,phi2,y,mu,u1,u2,rho){
D=nrow(u2)
d=nrow(u1)
t=D/d
k=ncol(u1)+1
A_d=list()
B_d=list()
B2_d=list()
D_d=list()
S_beta_d=list()
S_u2_d=rep(0,D*(k-1))
dim(S_u2_d)=c(D,(k-1))
S_u1=matrix(0,d,(k-1))
salida=list()
for(i in 1:D){
	S_beta_d[[i]]=crossprod((X[[i]]),(y[i,]-mu[i,]))
	S_u2_d[i,]=(y[i,]-mu[i,])
}
S_beta=add(S_beta_d)
ome=list()
for(i in 1:(k-1)){
	Omegad=omega(t,k,rho,phi2)
	ome[[i]]=solve(Omegad[[1]][[i]])}
jj=1
tt=1
S_u2_dd=matrix(0,D,(k-1))
for (i in 1:D){
	for (j in 1:(k-1)){
		S_u2_dd[i,j]=ome[[j]][tt,]%*%u2[jj:(jj+t-1),j]}
	tt=tt+1
	if (tt>t){
	tt=1
	jj=jj+t}}

S_u2=S_u2_d-(S_u2_dd/(matrix(rep(phi2,D),D,(k-1),byrow=TRUE)))
j=1
for (i in 1:d){
		S_u1[i,]=colSums(S_u2_d[j:(j+t-1),])
		j=j+t}
S_u1=S_u1-(u1/(matrix(rep(phi1,d),d,(k-1),byrow=TRUE)))

S_u11=matrix(t(S_u1),d*(k-1),1)
S_u22=matrix(t(S_u2),D*(k-1),1)
for(i in 1:D){
	A_d[[i]]=crossprod((X[[i]]),sigmap[[i]])%*%X[[i]]
	B_d[[i]]=crossprod((X[[i]]),sigmap[[i]])
	D_d[[i]]=crossprod((Z[[i]]),(sigmap[[i]]))
}

A=as(add(A_d[1:D]),"sparseMatrix")
#Hu2beta=B#
B=as(do.call(cbind,B_d),"sparseMatrix")
B2_d=addtolist(B_d,t,d)
#Hu1beta=B2#
B2=as(do.call(cbind,B2_d),"sparseMatrix")

#Hu2u2=DD#
DD=bdiag(sigmap)

D2_d=addtolist(D_d,t,d)
#Hu2u1=D2#
D2=as(do.call(cbind,D2_d),"sparseMatrix")

C2_d=addtolist(D_d,t,d)
C2=do.call(cbind,C2_d)
#Hu1u1=C22#
C22=addtomatrix(C2,d,t,k)

#H2u2u2#
for (i in 1:(k-1)){
	ome[[i]]=ome[[i]]/phi2[i]}

a=matrix(0,(t*(k-1)),(t*(k-1)))
tt2=1
kk2=1

v=list()
for (i in 1:(t*(k-1))){
	tt1=1
	kk1=1
	for (j in 1:(t*(k-1))){
		if (kk1 != kk2) {a[i,j]=0}
		if (kk1==kk2) {a[i,j]=ome[[kk1]][tt2,tt1]}
		kk1=kk1+1
		if (kk1>(k-1))
		{kk1=1
		tt1=tt1+1}
		if (tt1>(t)) {tt1=1}
		}
		kk2=kk2+1
		if (kk2>(k-1)) {kk2=1
		tt2=tt2+1}
		if (tt2>(t)) {tt2=1}}

for (i in 1:d){
	v[[i]]=a}

DDu2=bdiag(v)
C22=as(C22,"sparseMatrix")+as(diag(rep((1/phi1),d)),"sparseMatrix")
DD=DD+DDu2
rm(B2_d,A_d,C2,B_d,D_d,D2_d,DDu2,S_beta_d,S_u2_dd, S_u2_d,S_u1,S_u2,C2_d,v)
Fuu=as(rBind(cBind(C22,t(D2)),cBind(D2,DD)),"sparseMatrix")
Fuu=as(solve(Fuu),"sparseMatrix")
int=as(solve(DD),"sparseMatrix")
rm(DD)
Fuubb=rBind(t(B2),t(B))
A=solve(A-t(Fuubb)%*%Fuu%*%(Fuubb))
Fub=-1*tcrossprod(A,Fuubb)%*%Fuu
Fu1u1=solve(C22-t(D2)%*%int%*%D2)
rm(C22,Fuu,B,B2)
Fu1u2=-1*tcrossprod(Fu1u1,D2)%*%int
Fu2u2=int+int%*%D2%*%tcrossprod(Fu1u1,D2)%*%int
Fiuu=rBind(cBind(Fu1u1,Fu1u2),cBind(t(Fu1u2),Fu2u2))
rm(int,D2)
Fuu=Fiuu+Fiuu%*%Fuubb%*%tcrossprod(A,Fuubb)%*%Fiuu
F=rBind(cBind(A,(Fub)),cBind(t(Fub),Fuu))
S=rbind(S_beta,rbind(S_u11,S_u22))
rm(A,Fub,Fu1u1,Fu1u2,Fu2u2,Fiuu,Fuubb,S_beta,S_u11,S_u22,Fuu)

F=as(F,"matrix")
fisher=list(F=F,S=S)
return(fisher)
}

####################################################################
#' Variance components for Model 3
#'
#' This function calculates the variance components for the multinomial mixed model with two independent random effects
#' for each category of the response variable: one domain random effect and another correlated time and domain random effect (Model 3). This variance components
#' are used in the second part of the fitting algorithm
#' implemented in \code{\link[mme]{modelfit3}}. The algorithm adapts the ideas of Schall (1991) to a multivariate model. The variance components are
#' estimated by the REML method.
#'
#' @param p vector with the number of auxiliary variables per category.
#' @param sigmap a list with the model variance-covariance matrices for each domain obtained from \code{\link[mme]{wmatrix}}.
#' @param X list of matrices with the auxiliary variables obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of categories of the response variable minus one.
#' @param theta matrix with the estimated log-probabilites of each category in front of the reference category obtained from \code{\link[mme]{prmu.time}}.
#' @param phi1 vector with the initial values of the first variance component obtained from \code{\link[mme]{modelfit3}}.
#' @param phi2 vector with the initial values of the second variance component obtained from \code{\link[mme]{modelfit3}}.
#' @param u1 matrix with the values of the first random effect obtained from \code{\link[mme]{modelfit3}}.
#' @param u2 matrix with the values of the second random effect obtained from \code{\link[mme]{modelfit3}}.
#' @param rho vector with the initial values of the correlation parameter obtained from \code{\link[mme]{modelfit3}}.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#' \code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult.ct}},
#' \code{\link[mme]{prmu.time}}, \code{\link[mme]{Fbetaf.ct}}
#' \code{\link[mme]{sPhikf.ct}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit3}}, \code{\link[mme]{msef.ct}},
#' \code{\link[mme]{mseb}}, \code{\link[mme]{omega}}
#' @return a list containing the following components.
#' \item{phi1.new}{vector with the values of the variance component for the first random effect.}
#' \item{phi2.new}{vector with the values of the variance component for the second random effect.}
#' \item{rho.new}{vector with the correlation parameter.}
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Small area estimation of labour force indicator under a multinomial mixed model
#' with correlated time and area effects. Submitted for review.
#' @references Schall, R (1991). Estimation in generalized linear models with
#' random effects. Biometrika,
#' 78,719-727.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' mod=3 #type of model
#' data(simdata3) #data
#' datar=data.mme(simdata3,k,pp,mod)
#' initial=datar$initial
#' mean=prmu.time(datar$n,datar$Xk,initial$beta.0,initial$u1.0,initial$u2.0)
#' sigmap=wmatrix(datar$n,mean$estimated.probabilities)
#'
#' ##The variance components
#' phi.ct=phi.direct.ct(pp,sigmap,datar$X,mean$eta,initial$phi1.0,
#'        initial$phi2.0,initial$u1.0,initial$u2.0,initial$rho.0)


phi.direct.ct<- function(p,sigmap,X,theta,phi1,phi2,u1,u2,rho){
D=nrow(u2)
d=nrow(u1)
t=D/d
k=ncol(u1)+1
r=sum(p+1)
F=matrix(0,t,t)
u=rep(t,d)
mdcum=cumsum(u)
a=list()
a[[1]]=1:t
for(i in 2:d){
a[[i]] <- (mdcum[i-1]+1):mdcum[i]}

Xd=list()
for (i in 1:d){
	Xd[[i]]=X[a[[i]]]
	Xd[[i]]=do.call(rbind,Xd[[i]])
}

Xt=as(do.call(rbind,Xd),"sparseMatrix")
rm(Xd,a)
omegaa=omega(t,k,rho,phi2)

V22=list()
V1=list()
V2=list()
for (i in 1:(k-1)){
	V1[[i]]=(phi1[i])*diag(1,d)
	for (j in 1:d){
		V22[[j]]=phi2[i]*omegaa[[1]][[i]]}
	V2[[i]]=bdiag(V22)}


W=bdiag(sigmap)
AA=list()
for(i in 1:D){
	AA[[i]]=matrix(0,(k-1),(k-1)*D)
	for (j in 1:(k-1)){
		AA[[i]][j,((j-1)*D+i)]=1}
	}
Z2=as(do.call(rbind,AA),"sparseMatrix")

AA=list()
for(i in 1:d){
	AA[[i]]=matrix(0,t*(k-1),(k-1)*d)
	kk=1
	for (j in 1:(t*(k-1))){
		AA[[i]][j,((kk-1)*d+i)]=1
		kk=kk+1
		if (kk>(k-1)) {kk=k-1}}
	}
Z1=as(do.call(rbind,AA),"sparseMatrix")

rm(AA)

V=Z1%*%bdiag(V1)%*%t(Z1)+Z2%*%bdiag(V2)%*%t(Z2)+solve(W)
V=solve(V)

qq=t(Xt)%*%V%*%Xt
qq=solve(qq)

rm(V,V22)
T1=solve(t(Z1)%*%W%*%Z1+solve(bdiag(V1)))
T2=solve(t(Z2)%*%W%*%Z2+solve(bdiag(V2)))

T1=T1+T1%*%t(Z1)%*%W%*%Xt%*%qq%*%t(Xt)%*%W%*%Z1%*%T1
T2=T2+T2%*%t(Z2)%*%W%*%Xt%*%qq%*%t(Xt)%*%W%*%Z2%*%T2
rm(W,Xt,qq,Z1,Z2)
j=1
#T1=as(T1,"matrix")
#T2=as(T2,"matrix")
tau=rep(0,(k-1))
Trmll=list()
nr=nrow(T1)
for (i in 1:(k-1)){
	Trmll[[i]]=T1[j:(j+(nr/(k-1))-1),j:(j+(nr/(k-1))-1)]
	tau[i]=sum(diag(Trmll[[i]]))/phi1[i]
	j=j+(nr/(k-1))
}

tau=as.vector(tau)
phi1.new=diag((t(u1)%*%u1)/diag(d-tau))
nr=nrow(T2)
diagonal=c(0,rep(1,(t-2)),0)
E=diag(diagonal)
F=matrix(0,t,t)
F[lower.tri(F)]=(sequence((t-1):1))
F[F>1]=0
F=-F
F=F+t(F)
diag(F)=0
EE=list()
FF=list()
for( i in 1:d){
	EE[[i]]=E
	FF[[i]]=F}
E=bdiag(EE)
F=bdiag(FF)
num=rep(0,(k-1))
den=rep(0,(k-1))

j=1
for (i in 1:(k-1)){
	Trmll[[i]]=T2[j:(j+(nr/(k-1))-1),j:(j+(nr/(k-1))-1)]
	tau[i]=sum(diag(solve(V2[[i]]/phi2[i])%*%Trmll[[i]]))/phi2[i]
	num[i]=sum(diag(Trmll[[i]]%*%F))
	den[i]=sum(diag(Trmll[[i]]%*%E))
	j=j+(nr/(k-1))
}

den=as.vector(den)
num=as.vector(num)
phi2.new=rep(0,(k-1))
for (i in 1:(k-1)){
	phi2.new[i]=as.vector((t(u2[,i])%*%solve((V2[[i]])/phi2[i])%*%u2[,i])/((d*t)-tau[i]))
	div=as.vector((2*d)/(1-rho[i]^2)+2*(1/phi2[i])*(den[i]+(t(u2[,i])%*%E%*%u2[,i])))

	rho[i]=as.vector((-1/phi2[i])*(num[i]+(t(u2[,i])%*%F%*%u2[,i]))/div)}

rm(Trmll,tau, T1,T2,V1,V2,num,den,E,F,EE,FF)
resul=list(phi1.new=phi1.new,phi2.new=phi2.new,rho.new=rho)

}


####################################################################
#' Model correlation matrix for Model 3
#'
#' This function calculates the model correlation matrix and the first derivative of the model correlation matrix for Model 3. Model 3 is the multinomial mixed model with two independent random effects
#' for each category of the response variable: one domain random effect and another correlated time and domain random effect.
#'
#' @param t number of time periods.
#' @param k number of categories of the response variable.
#' @param rho vector with the correlation parameter obtained from \code{\link[mme]{modelfit3}}.
#' @param phi2 vector with the values of the second variance component obtained from \code{\link[mme]{modelfit3}}.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#'\code{\link[mme]{wmatrix}} ,\code{\link[mme]{phi.mult.ct}},
#' \code{\link[mme]{prmu.time}}, \code{\link[mme]{phi.direct.ct}},
#' \code{\link[mme]{Fbetaf.ct}}, \code{\link[mme]{sPhikf.ct}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit3}}, \code{\link[mme]{msef.ct}},
#' \code{\link[mme]{mseb}}
#' @return A list containing the following components.
#' \item{Omega.d}{correlation matrix.}
#' \item{First.derivative.Omegad}{Fisher derivative of the model correlation matrix.}
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Small area estimation of labour force indicator under a multinomial mixed model
#' with correlated time and area effects. Submitted for review.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' mod=3 #type of model
#' data(simdata3)   #data
#' datar=data.mme(simdata3,k,pp,mod)
#' initial=datar$initial
#' mean=prmu.time(datar$n,datar$Xk,initial$beta.0,initial$u1.0,initial$u2.0)
#'
#' ##The model correlation matrix
#' matrix.corr=omega(datar$t,k,initial$rho.0,initial$phi2.0)

omega<-function(t,k,rho,phi2){
ome=list()
Omegadd=list()
Omegadfirsti=list()
for(i in 1:(k-1)){
	Omegad=matrix(0,t,t)
	Omegad[lower.tri(Omegad)]=rho[i]^(sequence((t-1):1))
	Omegad=Omegad+t(Omegad)
	diag(Omegad)=1
	Omegad=(1/(1-rho[i]^2))*Omegad
	Omegadd[[i]]=Omegad
	ome[[i]]=phi2[i]*Omegad

	OmegadFirst<-matrix(0,t,t)
	OmegadFirst[lower.tri(OmegadFirst)]<-sequence((t-1):1)*rho[i]^(sequence((t-1):1)-1)
	OmegadFirst<-OmegadFirst+t(OmegadFirst)
	OmegadFirst<- (1/(1-rho[i]^2))*OmegadFirst
	OmegadFirst <- OmegadFirst + (2*rho[i]/(1-rho[i]^2))*Omegad
	Omegadfirsti[[i]] <- phi2[i]*OmegadFirst}
salida=list(Omega.d=Omegadd,First.derivative.Omegad=Omegadfirsti)
return(salida)
}


####################################################################
#' Fisher information matrix and score vectors of the variance components for Model 3
#'
#' This function computes the Fisher information matrix and the score vectors
#' of the variance components, for the multinomial mixed model with two independent random effects
#' for each category of the response variable: one domain random effect and another correlated time and domain random effect (Model 3).
#' These values are used in the fitting algorithm implemented in \code{\link[mme]{modelfit3}} to estimate the random effects. The algorithm adatps the
#' ideas of Schall (1991) to a multivariate
#' model. The variance components are estimated by the REML method.
#'
#' @param d number of areas.
#' @param t number of time periods.
#' @param pp vector with the number of the auxiliary variables per category.
#' @param sigmap a list with the model variance-covariance matrices for each domain obtained from \code{\link[mme]{wmatrix}}.
#' @param X list of matrices with the auxiliary variables obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of categories of the response variable minus one.
#' @param eta matrix with the estimated log-rates of probabilites of each category over the reference category obtained from \code{\link[mme]{prmu.time}}.
#' @param phi1 vector with the values of the first variance component obtained from \code{\link[mme]{modelfit3}}.
#' @param phi2 vector with the values of the second variance component obtained from \code{\link[mme]{modelfit3}}.
#' @param rho vector with the correlation parameter obtained from \code{\link[mme]{modelfit3}}.
#' @param pr matrix with the estimated probabilities of the response variable obtained from \code{\link[mme]{prmu.time}}.
#' @param M vector with the area sample sizes.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#' \code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult.ct}},
#' \code{\link[mme]{prmu.time}}, \code{\link[mme]{phi.direct.ct}},
#' \code{\link[mme]{Fbetaf.ct}}, \code{\link[mme]{omega}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit3}}, \code{\link[mme]{msef.ct}},
#' \code{\link[mme]{mseb}}
#' @return A list containing the following components.
#' \item{S}{(phi1, phi2, rho) score vector.}
#' \item{F}{Fisher information matrix of the variance components (phi1, phi2, rho).}
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Small area estimation of labour force indicators under a multinomial mixed model
#' with correlated time and area effects. Submitted for review.
#' @references Schall, R (1991). Estimation in generalized linear models with
#' random effects. Biometrika,
#' 78,719-727.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' mod=3 #type of model
#' data(simdata3) #data
#' datar=data.mme(simdata3,k,pp, mod)
#' initial=datar$initial
#' mean=prmu.time(datar$n,datar$Xk,initial$beta.0,initial$u1.0,initial$u2.0)
#' sigmap=wmatrix(datar$n,mean$estimated.probabilities)
#'
#' ## Fisher information matrix and the score vectors
#' Fisher.phi.ct=sPhikf.ct(datar$d,datar$t,pp,sigmap,datar$X,mean$eta,initial$phi1.0,
#'            initial$phi2.0,initial$rho.0,mean$estimated.probabilities,datar$n)


sPhikf.ct<- function(d,t,pp,sigmap,X,eta,phi1,phi2,rho,pr,M){
k=ncol(eta)+1
r=sum(pp+1)
Vd=list()
qq=matrix(0,r,r)

u=rep(t,d)
mdcum=cumsum(u)
a=list()
a[[1]]=1:t
for(i in 2:d){
a[[i]] <- (mdcum[i-1]+1):mdcum[i]}

Xd=list()
for (i in 1:d){
	Xd[[i]]=X[a[[i]]]
	Xd[[i]]=do.call(rbind,Xd[[i]])
}
thetaa=list()
j=1
for (i in 1:d){
	thetaa[[i]]=matrix((eta[(j:(j+t-1)),]),((k-1)*t),1)
		j=j+t}
j=1
Wd1_1=list()
Wd1=list()
for (i in 1:d){
	for (jj in 1:(k-1)){
		Wd1_1[[jj]]=diag(M[(j:(j+t-1))]*(pr[(j:(j+t-1)),jj]))}
	Wd1[[i]]=do.call("blockdiag",Wd1_1)
	j=j+t}

js=1
Wd2_2=list()
Wd2_22=list()
Wd2=list()
for (i in 1:d){
	jj=1
	for (j in 1:(k-1)){
		for (l in 1:(k-1)){
		Wd2_2[[l]]=-diag(M[(js:(js+t-1))]*(pr[(js:(js+t-1)),j])*(pr[(js:(js+t-1)),l]))
		jj=jj+1}
		Wd2_22[[j]]=do.call(cbind,Wd2_2)}
		Wd2[[i]]=do.call(rbind,Wd2_22)
		js=js+t}
for (i in 1:d){
	Wd1[[i]]=Wd2[[i]]+Wd1[[i]]}
ome=list()
Omegadd=list()
Omegaddd=list()
for(i in 1:(k-1)){
	a=omega(t,k,rho,phi2)
	ome[[i]]=phi2[i]*a[[1]][[i]]
	Omegadd[[i]]=a[[2]][[i]]
}
Vu2=list()
for (i in 1:d){
	Vu2[[i]]=do.call("blockdiag",ome)
	Omegaddd[[i]]=do.call("blockdiag",Omegadd)}

j=1
I1=list()
for (i in 1:(k-1)){
	I1[[i]]=matrix(rep(1,t),t,1)}
I1=do.call("blockdiag",I1)
for (i in 1:d){
	Vd[[i]]=I1%*%diag(phi1)%*%t(I1)+Vu2[[i]]+solve(Wd1[[i]])
	Vd[[i]]=solve(Vd[[i]])
	qq=qq+crossprod((Xd[[i]]),(Vd[[i]]))%*%Xd[[i]]}

qq=solve(qq)

S1=rep(0,3*(k-1))
S2=rep(0,3*(k-1))

delta11=list()
sum1=list()
sum2=list()
sum3=list()
sum4=list()
for (i in 1:(3*(k-1))){
	sum1[[i]]=matrix(0,1,r)
	sum2[[i]]=matrix(0,r,1)
	sum3[[i]]=matrix(0,1,r)
	sum4[[i]]=matrix(0,r,r)}


dia=list()
dia2=list()
o=1
for(i in 1:3){
	for (j in 1:(k-1)){

			delta=matrix(0,(k-1),(k-1))
			for (ll in 1:(k-1)){
				dia[[ll]]=matrix(0,t,t)
				dia2[[ll]]=matrix(0,t,t)
				if (ll==j) {
				delta[ll,ll]=1
				a=omega(t,k,rho,phi2)
				dia[[ll]]=a[[1]][[ll]]
				dia2[[ll]]=a[[2]][[ll]]}}

			if (i==1) {Vakd=I1%*%delta%*%t(I1)}
			if (i==2) {Vakd=do.call("blockdiag",dia)}
			if (i==3) {Vakd=do.call("blockdiag",dia2)}

			for (l in 1:d){
			S1[o]=S1[o]+sum(diag((Vd[[l]]-Vd[[l]]%*%Xd[[l]]%*%qq%*%t(Xd[[l]])%*%Vd[[l]])%*%Vakd))
			S2[o]=S2[o]+t(thetaa[[l]])%*%Vd[[l]]%*%Vakd%*%(Vd[[l]])%*%thetaa[[l]]
			sum1[[o]]=sum1[[o]]+t(thetaa[[l]])%*%Vd[[l]]%*%Vakd%*%(Vd[[l]])%*%Xd[[l]]
			sum2[[o]]=sum2[[o]]+t(Xd[[l]])%*%Vd[[l]]%*%thetaa[[l]]
			sum3[[o]]=sum3[[o]]+t(thetaa[[l]])%*%Vd[[l]]%*%Xd[[l]]
			sum4[[o]]=sum4[[o]]+t(Xd[[l]])%*%Vd[[l]]%*%Vakd%*%(Vd[[l]])%*%Xd[[l]]}
			o=o+1}}
tt=rep(0,(3*(k-1)))
for(i in 1:(3*(k-1))){
	tt[i]=-2*(sum1[[i]]%*%qq%*%sum2[[i]])+sum3[[i]]%*%qq%*%sum4[[i]]%*%qq%*%sum2[[i]]}

S2=S2+tt
S=-0.5*S1+0.5*S2

F1=rep(0,((k-1)*(k-1)*9))
F2=rep(0,((k-1)*(k-1)*9))
F3=rep(0,((k-1)*(k-1)*9))
F4=rep(0,((k-1)*(k-1)*9))
aa=1
o=1
for (i in 1:3){
	for(j in 1:3){
		for(ii in 1:(k-1)){
			for (jj in 1:(k-1)){
				delta=matrix(0,(k-1),(k-1))
				for (ll in 1:(k-1)){
					dia[[ll]]=matrix(0,t,t)
					dia2[[ll]]=matrix(0,t,t)
					if (ll==ii) {
						delta[ll,ll]=1
						dia[[ll]]=a[[1]][[ll]]
						dia2[[ll]]=a[[2]][[ll]]
				}	}

				if (i==1) {Vakd1=I1%*%delta%*%t(I1)}
				if (i==2) {Vakd1=do.call("blockdiag",dia)}
				if (i==3) {Vakd1=do.call("blockdiag",dia2)}

				delta=matrix(0,(k-1),(k-1))
				for (ll in 1:(k-1)){
						dia[[ll]]=matrix(0,t,t)
						dia2[[ll]]=matrix(0,t,t)
						if (ll==jj) {
						delta[ll,ll]=1
						dia[[ll]]=a[[1]][[ll]]
						dia2[[ll]]=a[[2]][[ll]]
						}}

				if (j==1) {Vakd2=I1%*%delta%*%t(I1)}
				if (j==2) {Vakd2=do.call("blockdiag",dia)}
				if (j==3) {Vakd2=do.call("blockdiag",dia2)}
				for (l in 1:d){
					F1[aa]=F1[aa]+0.5*sum(diag(Vd[[l]]%*%Vakd1%*%Vd[[l]]%*%Vakd2))
					F2[aa]=F2[aa]+0.5*sum(diag(Vd[[l]]%*%Vakd1%*%Vd[[l]]%*%Xd[[l]]%*%qq%*%t(Xd[[l]])%*%Vd[[l]]%*%Vakd2))
					F3[aa]=F3[aa]+0.5*sum(diag(Vd[[l]]%*%Xd[[l]]%*%qq%*%t(Xd[[l]])%*%Vd[[l]]%*%Vakd1%*%Vd[[l]]%*%Vakd2))
					F4[aa]=F3[aa]+0.5*sum(diag(Vd[[l]]%*%Xd[[l]]%*%qq%*%sum4[[o]]%*%qq%*%t(Xd[[l]])%*%Vd[[l]]%*%Vakd2))}
					aa=aa+1
				}}
					}
					o=o+1}

F=F1+F2+F3+F4

FF=list()
FFF=list()
o=1
for (i in 1:3){
	for (j in 1:3){
		FF[[j]]=matrix(F[o:(o+(k-1)*(k-1)-1)],(k-1),(k-1))
		o=o+(k-1)*(k-1)}
		FFF[[i]]=do.call(cbind,FF)}
F=do.call(rbind,FFF)
F=solve(F)
S=matrix(t(S),3*(k-1),1)
A=F%*%S
rm(sum1,sum2,sum3,sum4, Vakd1,Vakd2, Vakd, S1,S2,Xd,Vd,thetaa, I1)

score.F=list(S=S,F=F)
return(score.F)
}

####################################################################
####################################################################
#' Function used to fit Model 3
#'
#' This function fits the multinomial mixed model with two independent random effects
#' for each category of the response variable: one domain random effect and another correlated time and domain random effect (Model 3). The formulation is described in Lopez-Vizcaino et al. (2013).
#' The fitting algorithm combine the penalized quasi-likelihood method (PQL) for estimating
#' and predicting the fixed and random effects, respectively, with the residual maximun likelihood method (REML)
#' for estimating the variance components. This function uses as initial values the output of the function
#' \code{\link[mme]{initial.values}}.
#'
#' @param d number of areas.
#' @param t number of time periods.
#' @param pp vector with the number of the auxiliary variables per category.
#' @param Xk list of matrices with the auxiliary variables per category obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of domains.
#' @param X list of matrices with the auxiliary variables obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of categories of the response variable minus one.
#' @param Z design matrix  of random effects obtained from \code{\link[mme]{data.mme}}.
#' @param initial output of the function \code{\link[mme]{initial.values}}.
#' @param y matrix with the response variable obtained from \code{\link[mme]{data.mme}}, except the reference category. The rows are the domains and the columns are the categories of the response variable minus one.
#' @param M vector with the area sample sizes.
#' @param MM vector with the population sample sizes.
#' @param b parameter that indicates the bootstrap.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#' \code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult.ct}},
#' \code{\link[mme]{prmu.time}}, \code{\link[mme]{phi.direct.ct}},
#' \code{\link[mme]{sPhikf.ct}}, \code{\link[mme]{omega}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{Fbetaf.ct}}, \code{\link[mme]{msef.ct}},
#' \code{\link[mme]{mseb}}
#' @return A list containing the following components.
#' \item{Estimated.probabilities}{matrix with the estimated probabilities
#' for the categories of response variable.}
#' \item{Fisher.information.matrix.phi}{Fisher information matrix of phi.}
#' \item{Fisher.information.matrix.beta}{Fisher information matrix of beta.}
#' \item{u1}{matrix with the estimated first random effect.}
#' \item{u2}{matrix with the estimated second random effect.}
#' \item{mean}{matrix with the estimated mean of the response variable.}
#' \item{warning1}{0=OK,1=The model could not be fitted.}
#' \item{warning2}{0=OK,1=The value of the variance component is negative: the initial value
#' is taken.}
#' \item{beta.Stddev.p.value}{matrix with the estimated fixed effects, its standard
#' deviations and its p-values.}
#' \item{phi.Stddev.p.value}{matrix with the estimated variance components, its
#' standard deviations and its p-values.}
#' \item{rho}{estimated correlation parameter.}
#' \item{rho.Stddev.p.value}{matrix with the estimated correlation parameter, its
#' standard deviations and its p-values.}
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Small area estimation of labour force indicators under a multinomial mixed model
#' with correlated time and area effects. Submitted for review.
#' @examples
#' \dontrun{
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' mod=3 #type of model
#' data(simdata3)   #data
#' datar=data.mme(simdata3,k,pp,mod)
#'
#' ##Model fit
#' result=modelfit3(datar$d,datar$t,pp,datar$Xk,datar$X,datar$Z,datar$initial,datar$y[,1:(k-1)],
#'        datar$n,datar$N,0)}

modelfit3<-function(d,t,pp,Xk,X,Z,initial,y,M,MM,b){
maxiter1=50
maxiter2=50
iter1=0
eps=1e-3
D=d*t
beta.new=initial[[1]]
phi1.new=initial[[2]]
phi2.new=initial[[3]]
u1.new=initial[[4]]
u2.new=initial[[5]]
rho=initial[[6]]
k=ncol(u1.new)+1
phi1.prim=phi1.new
phi2.prim=phi2.new
rho.prim=rho
p=sum(pp+1)
cont=0
ii=1
ll=0
aviso=0
Fk=matrix(0,3*(k-1),3*(k-1))
mm=1

resulbeta=matrix(0,sum(pp+1),3)
resulphi=matrix(0,2*(k-1),3)
resulrho=matrix(0,(k-1),3)
while (iter1<maxiter1) {
	phi1.old=phi1.new
	phi2.old=phi2.new
	rho.old=rho
	iter2=0
	p1=matrix(1,(p+1),(k-1))
	p2=matrix(1,D,(k-1))
	while(iter2<maxiter2 & mm>eps){
		beta.old=beta.new
		u1.old=u1.new
		u2.old=u2.new
		prmul=prmu.time(M,Xk,beta.old,u1.old,u2.old)
		theta=prmul[[3]]
		pr=prmul[[1]]
		mu=prmul[[2]]
		sigmap=wmatrix(M,pr)
		for (i in 1:D){
			comp=det(sigmap[[i]])
			if (is.na(comp)){cont=1}
			if (is.na(comp)==FALSE & (abs(comp))<0.000001 ) {cont=1}}
		if (cont==0){
		inversaFbeta=Fbetaf.ct(sigmap,X,Z,phi1.old,phi2.old,y,mu,u1.old,u2.old,rho.old)
		F=inversaFbeta[[1]]
		S=inversaFbeta[[2]]
		beta_u_old=rbind(do.call(rbind,beta.old),rbind(matrix(t(u1.old),d*(k-1),1),matrix(t(u2.old),D*(k-1),1)))

		beta_u_new=beta_u_old+F%*%S
		p1=(beta_u_new-beta_u_old)/beta_u_old
		mm=max(abs(p1))
		if (mm<eps){iter2=maxiter2} else {iter2=iter2+1}
		ucum=cumsum(pp+1)
		beta.new=list()
		beta.new[[1]]=as.matrix(beta_u_new[(1:ucum[1]),1])
		for (i in 2:(k-1)){
			beta.new[[i]]=as.matrix(beta_u_new[((ucum[i-1]+1):ucum[i]),1])
		}
		u1.new=matrix(beta_u_new[(ucum[k-1]+1):(ucum[k-1]+d*(k-1)),],d,k-1,byrow=TRUE)
		u2.new=matrix(beta_u_new[(ucum[k-1]+1+d*(k-1)):nrow(beta_u_new),],D,k-1,byrow=TRUE)
		if (ii==1) {
		beta.prim=beta.new
		F.prim=F
		u1.prim=u1.new
		u2.prim=u2.new}}

		if (cont==1){iter2=maxiter2}

	}
	if (cont==0){

	prmul=prmu.time(M,Xk,beta.new,u1.new,u2.new)
	theta=prmul[[3]]
	pr=prmul[[1]]
	mu=prmul[[2]]
	sigmap=wmatrix(M,pr)

	#sk_F=sPhikf.ct(d,t,pp,sigmap,X,theta,phi1.old,phi2.old,rho.old,pr,M)
	#sk=sk_F[[1]]
	#Fk=sk_F[[2]]
	if (ii==1){Fk.prim=Fk}
	ii=ii+1
	phi.old=rbind(data.matrix(phi1.old),data.matrix(phi2.old),data.matrix(rho.old))
	resul=phi.direct.ct(pp,sigmap,X,theta,phi1.old,phi2.old,u1.new,u2.new,rho.old)
	phi1.new=resul[[1]]
	phi2.new=resul[[2]]
	rho=resul[[3]]
	phi.new=matrix(cbind(phi1.new,phi2.new,rho),(3*(k-1)),1)
  #phi.new=phi.old+F%*%sk
	phi.new=data.matrix(phi.new)
	p3=(phi.new-phi.old)/phi.old
	phi1.new=phi.new[1:(k-1),]
	phi2.new=phi.new[k:(2*(k-1)),]
	rho=phi.new[(2*k-1):(3*(k-1)),]
	mmm=max(abs(c(p1,p3)))
	if (min(phi1.new,phi2.new)<0.0001 | max(abs(rho))>0.99 ){
	aviso=1
	ll=ll+1
	phi1.new=phi1.prim
	phi2.new=phi2.prim
	#beta.new=beta.prim
	u1.new=u1.prim
	u2.new=u2.prim
	rho=rho.prim
	#F=F.prim
	if (ll>1){
	#Fk=Fk.prim
	phi1.new=phi1.prim
	phi2.new=phi2.prim
	#beta.new=beta.prim
	rho=rho.prim
	iter1=maxiter1}
	iter1=iter1+1
		}
	if (aviso==0){
	mmm=max(abs(c(p1,p3)))
	mm=1}
	else
	{
	mmm=max(abs(c(p3)))
	mm=eps/10}
	if (mmm>eps){iter1=iter1+1} else {iter1=maxiter1}}
	if (cont==1) {
	iter1=maxiter1}

	}
if (cont==0){
prmul=prmu.time(MM,Xk,beta.new,u1.new,u2.new)
mu=prmul[[2]]
pr=prmul[[1]]
colnames(mu)=paste("Yest",1:(k-1),sep="")
b2=do.call(rbind,beta.new)
cii=ci(b2,F[1:(sum(pp+1)),1:(sum(pp+1))])
resulbeta=cbind(Beta=b2,Std.dev=cii[[1]],p.value=cii[[2]])
colnames(resulbeta)=c("Estimate","Std.Error","p.value")
u=list()
for (i in 1:(k-1)){
		u[[i]]=as.matrix(colnames(Xk[[i]]))
	u[[i]][1]="Intercept"}
u=do.call(rbind,u)
rownames(resulbeta)=u

if (b==0){
sk_F=sPhikf.ct(d,t,pp,sigmap,X,theta,phi1.new,phi2.new,rho,pr,M)
Fk=sk_F[[2]]
phi.new=(rbind(as.matrix(phi1.new),as.matrix(phi2.new)))
cii=ci(phi.new,Fk[(1:(2*(k-1))),(1:(2*(k-1)))])
resulphi=cbind(phi.est=phi.new,Std.dev=cii[[1]],p.value=cii[[2]])
colnames(resulphi)=c("Estimate","Std.dev","p.value")
cii=ci(rho,Fk[(2*(k-1)+1):(3*(k-1)),(2*(k-1)+1):(3*(k-1))])
resulrho=cbind(rho.est=rho,Std.dev=cii[[1]],p.value=cii[[2]])
colnames(resulrho)=c("Estimate","Std.Error","p.value")}
if (b==1) {
phi.new=(rbind(as.matrix(phi1.new),as.matrix(phi2.new)))
resulphi=cbind(phi.est=phi.new)
colnames(resulphi)=c("Estimate")
resulrho=cbind(rho.est=rho)
colnames(resulrho)=c("Estimate")
}
}
result=list()
result=list(Estimated.probabilities=pr,u1=u1.new,u2=u2.new,mean=mu,warning1=cont,Fisher.information.matrix.beta=F,Fisher.information.matrix.phi=Fk,beta.Stddev.p.value=resulbeta,phi.Stddev.p.value=resulphi,warning2=aviso,rho=rho,rho.Stddev.p.value=resulrho)
class(result)="mme"
return(result)
}

####################################################################
#' Analytic MSE for Model 3
#'
#' This function calculates the analytic MSE for the multinomial mixed model with two independent random effects
#' for each category of the response variable: one  random effect associated with the domain and another correlated random effect associated with time and domain (Model 3). See details of the model and the expresion of mse in Lopez-Vizcaino et al. (2013). The formulas
#' of Prasad and Rao (1990) are adapted to Model 3. This function uses the output of \code{\link[mme]{modelfit3}}.
#'
#' @param p vector with the number of the auxiliary variables per category.
#' @param X list of matrices with the auxiliary variables obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of categories of the response variable minus one.
#' @param result the output of the function \code{\link[mme]{modelfit3}}.
#' @param M vector with the area sample sizes.
#' @param MM vector with the population sample sizes.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#'\code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult.ct}},
#' \code{\link[mme]{prmu.time}}, \code{\link[mme]{phi.direct.ct}},
#' \code{\link[mme]{sPhikf.ct}}, \code{\link[mme]{modelfit3}},
#' \code{\link[mme]{Fbetaf.ct}}, \code{\link[mme]{ci}}, \code{\link[mme]{omega}},
#' \code{\link[mme]{mseb}}.
#' @return mse.analitic is a matrix with the MSE estimator calculated by adapting the explicit
#' formulas of Prasad and  Rao (1990).
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Small area estimation of labour force indicators under a multinomial mixed model
#' with correlated time and area effects. Submitted for review.
#' @references Prasad, NGN, Rao, JNK (1990).The estimation of the mean squared error of small
#' area estimators. Journal of the American Statistical Association, 85, 163-171.
#' @examples
#' \dontrun{
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' mod=3 #type of model
#' data(simdata3) #data
#' datar=data.mme(simdata3,k,pp,mod)
#' ##Model fit
#' result=modelfit3(d,t,pp,datar$Xk,datar$X,datar$Z,datar$initial,datar$y[,1:(k-1)],
#'        datar$n,datar$N,0)
#'
#' ##Analytic MSE
#' msef=msef.ct(pp,datar$X,result,datar$n,datar$N)}




msef.ct<-function(p,X,result,M,MM){
pr=result$Estimated.probabilities
k=ncol(pr)
phi1.old=result[[9]][1:(k-1),1]
phi2.old=result[[9]][k:(2*k-2),1]
rho=result$rho
F=result$Fisher.information.matrix.phi
k=ncol(pr)
D=nrow(MM)
d=nrow(result$u1)
t=D/d
r=sum(p+1)
Vd=list()
qq=as(matrix(0,(r),(r)),"sparseMatrix")
Ld=list()

j=1
Wd1_1=list()
Wd1=list()
Hd1=Hdd1=list()
for (i in 1:d){
	for (jj in 1:(k-1)){
		Wd1_1[[jj]]=diag(M[(j:(j+t-1))]*(pr[(j:(j+t-1)),jj]))
		Hdd1[[jj]]=diag(MM[(j:(j+t-1))]*(pr[(j:(j+t-1)),jj]))}
	Wd1[[i]]=do.call("blockdiag",Wd1_1)
	Hd1[[i]]=do.call("blockdiag",Hdd1)
	j=j+t}

js=1
Wd2_2=Hd2_2=list()
Wd2_22=Hd2_22=list()
Wd2=Hd2=list()
for (i in 1:d){
	jj=1
	for (j in 1:(k-1)){
		for (l in 1:(k-1)){
		Wd2_2[[l]]=-diag(M[(js:(js+t-1))]*(pr[(js:(js+t-1)),j])*(pr[(js:(js+t-1)),l]))
		Hd2_2[[l]]=-diag(MM[(js:(js+t-1))]*(pr[(js:(js+t-1)),j])*(pr[(js:(js+t-1)),l]))
		jj=jj+1}
		Wd2_22[[j]]=do.call(cbind,Wd2_2)
		Hd2_22[[j]]=do.call(cbind,Hd2_2)}
		Wd2[[i]]=do.call(rbind,Wd2_22)
		Hd2[[i]]=do.call(rbind,Hd2_22)
		js=js+t}
for (i in 1:d){
	Wd1[[i]]=as(Wd2[[i]]+Wd1[[i]],"sparseMatrix")
	Hd1[[i]]=as(Hd2[[i]]+Hd1[[i]],"sparseMatrix")
}


u=rep(t,d)
mdcum=cumsum(u)
a=list()
a[[1]]=1:t
for(i in 2:d){
a[[i]] <- (mdcum[i-1]+1):mdcum[i]}
Xd=list()
Xaa=list()
for (i in 1:d){
	Xaa[[i]]=X[a[[i]]]
	Xd[[i]]=as(do.call(rbind,Xaa[[i]]),"sparseMatrix")
}

ome=list()
Omegadd=list()
Omegaddd=list()
for(i in 1:(k-1)){
	a=omega(t,k,rho,phi2.old)
	ome[[i]]=phi2.old[i]*a[[1]][[i]]
	Omegadd[[i]]=a[[2]][[i]]
}
Vu2=list()
for (i in 1:d){
	Vu2[[i]]=bdiag(ome)
	Omegaddd[[i]]=bdiag(Omegadd)}

j=1
I1=list()
for (i in 1:(k-1)){
	I1[[i]]=matrix(rep(1,t),t,1)}
I1=bdiag(I1)
for (i in 1:d){
	Vd[[i]]=I1%*%as(diag(phi1.old),"sparseMatrix")%*%t(as(I1,"matrix"))+Vu2[[i]]+solve(Wd1[[i]])
	Vd[[i]]=solve(as(Vd[[i]],"sparseMatrix"))
	qq=qq+crossprod(as(Xd[[i]],"matrix"),as(Vd[[i]],"matrix"))%*%Xd[[i]]}
qq=solve(qq)


H=bdiag(Hd1)

z2=as(diag(rep(1,(t*(k-1)))),"sparseMatrix")

T11=list()
T22=list()
T12=list()
Z1i=list()
TT=list()
for (i in 1:d){
	T11[[i]]=as(diag(phi1.old)-diag(phi1.old)%*%t(I1)%*%Vd[[i]]%*%I1%*%diag(phi1.old),"sparseMatrix")
	T12[[i]]=as(-diag(phi1.old)%*%t(I1)%*%Vd[[i]]%*%z2%*%Vu2[[i]],"sparseMatrix")
	T22[[i]]=as(Vu2[[i]]-Vu2[[i]]%*%t(z2)%*%Vd[[i]]%*%z2%*%Vu2[[i]],"sparseMatrix")
	#TT[[i]]=as(I1%*%T11[[i]]%*%t(I1)+I1%*%T12[[i]]%*%t(z2)+z2%*%t(T12[[i]])%*%t(I1)+z2%*%t(T22[[i]])%*%t(z2),"sparseMatrix")
	Z1i[[i]]=as(diag(phi1.old),"sparseMatrix")
}
t1=matrix(1,t,1)
M11=list()
M22=list()
M12=list()
a=b=c=list()
for (j in 1:(k-1)){
	a[[j]]=phi1.old[j]*t1%*%t(t1)
	b[[j]]=phi1.old[j]*t1
	c[[j]]=phi1.old[j]%*%t(t1)}
for (i in 1:d){

	M11[[i]]=bdiag(a)-bdiag(b)%*%t(I1)%*%Vd[[i]]%*%I1%*%bdiag(c)
	M12[[i]]=-bdiag(b)%*%t(I1)%*%Vd[[i]]%*%z2%*%Vu2[[i]]
	M22[[i]]=Vu2[[i]]-Vu2[[i]]%*%t(z2)%*%Vd[[i]]%*%z2%*%Vu2[[i]]
	#TT[[i]]=I1%*%T11[[i]]%*%t(I1)+I1%*%T12[[i]]%*%t(z2)+z2%*%t(T12[[i]])%*%t(I1)+z2%*%t(T22[[i]])%*%t(z2)
	#Z1i[[i]]=Iq1
}
G1=list()
G11=list()
G12=list()
G22=list()
A21=list()
A22=list()
G2=list()
del=delt=list()
jj=1
for (i in 1:d){
	for(j in 1:t){
		delta=matrix(0,t,1)
		delta[j,1]=1
		for(ij in 1:(k-1)){
			del[[ij]]=delta
			delt[[ij]]=t(delta)}
		G11[[jj]]=bdiag(delt)%*%Hd1[[i]]%*%M11[[i]]%*%t(Hd1[[i]])%*%bdiag(del)
		G12[[jj]]=bdiag(delt)%*%Hd1[[i]]%*%M12[[i]]%*%t(Hd1[[i]])%*%bdiag(del)
		G22[[jj]]=bdiag(delt)%*%Hd1[[i]]%*%M22[[i]]%*%t(Hd1[[i]])%*%bdiag(del)
		G1[[jj]]=G11[[jj]]+G12[[jj]]+t(G12[[jj]])+G22[[jj]]
		A21[[jj]]=bdiag(delt)%*%Hd1[[i]]%*%Xd[[i]]
		A22[[jj]]=bdiag(delt)%*%Hd1[[i]]%*%(M11[[i]]+M12[[i]]+t(M12[[i]])+M22[[i]])%*%Wd1[[i]]%*%Xd[[i]]
		G2[[jj]]=(A21[[jj]]-A22[[jj]])%*%qq%*%(t(A21[[jj]])-t(A22[[jj]]))
		jj=jj+1}}

Vu1=bdiag(Z1i)
V=bdiag(Vd)
Z1=matrix(0,1,((k-1)*d))
Z2=matrix(0,1,((k-1)*D))
jj=1
for (i in 1:d){
	a=c(rep(0,(i-1)),1,rep(0,(d-i)))
	b=cbind(rbind(a,rep(0,d)),rbind(rep(0,d),a))
	for (j in 1:t){
		dd=c(rep(0,((jj)-1)),1,rep(0,(D-(jj))))
		e=cbind(rbind(dd,rep(0,D)),rbind(rep(0,D),dd))
		Z1=rbind(Z1,b)
		Z2=rbind(Z2,e)
		jj=jj+1}}
Z1=as(Z1[2:((k-1)*D+1),],"sparseMatrix")
Z2=as(Z2[2:((k-1)*D+1),],"sparseMatrix")

R2=bdiag(Vu2)%*%V
R1=Z1%*%Vu1%*%t(Z1)%*%V

#Esto hay que rehacerlo#
aa=nrow(R1)
bb=nrow(R2)
delta11=list()
Vakd=list()


dia=list()
dia2=list()
Lk=va=list()
jj=1
for(i in 1:3){
	for (j in 1:(k-1)){
			delta=matrix(0,(k-1),(k-1))
			for (ll in 1:(k-1)){
				dia[[ll]]=matrix(0,t,t)
				dia2[[ll]]=matrix(0,t,t)
				if (ll==j) {
				delta[ll,ll]=1
				a=omega(t,k,rho,phi2.old)
				dia[[ll]]=a[[1]][[ll]]
				dia2[[ll]]=a[[2]][[ll]]}}

			if (i==1) {Vakd=I1%*%as(delta,"sparseMatrix")%*%t(I1)}
			if (i==2) {Vakd=bdiag(dia)}
			if (i==3) {Vakd=bdiag(dia2)}
			for (l in 1:d){
				va[[l]]=Vakd}
			if (i==1){Lk[[jj]]=as((diag(rep(1,aa))-R1),"sparseMatrix")%*%bdiag(va)%*%V}
			if (i>1){Lk[[jj]]=as((diag(rep(1,bb))-R2),"sparseMatrix")%*%bdiag(va)%*%V}
			jj=jj+1
				}}


a=list()
A=list()
jj=1
for (i in 1:d){
	for(j in 1:t){
		delta=matrix(0,t,1)
		delta[j,1]=1
		for (l in 1:(k-1)){
			a[[l]]=delta}
		A[[jj]]=1/1000*cbind(matrix(0,nrow=k-1,ncol=t*(i-1)*(k-1)),t(as(bdiag(a),"matrix")),matrix(0,nrow=k-1,ncol=t*(k-1)*(d-i)))
		jj=jj+1}}



jj=1
aa=3*(k-1)
gg3=list()
g3=matrix(0,(d*t*(k-1)),(d*t*(k-1)))
for (i in 1:aa){
	for (j in 1:aa){
		g33=H%*%Lk[[i]]%*%V%*%t(Lk[[j]])%*%t(H)
		g3=g3+F[i,j]*g33
		jj=jj+1
		}}
for (kk in 1:D){
	gg3[[kk]]=A[[kk]]%*%g3%*%t(A[[kk]])}

g=list()
gg=list()
for (i in 1:D){
g[[i]]=G1[[i]]+G2[[i]]+2*gg3[[i]]
gg[[i]]=G1[[i]]+G2[[i]]}

mse.analitic=matrix(0,D,(k-1))
for (i in 1:(k-1)){
	for (j in 1:D){
	mse.analitic[j,i]=g[[j]][i,i]}}
colnames(mse.analitic)=paste("mse.",1:(k-1),sep="")
mse=list(mse.analitic=mse.analitic)
rm(A,G1,G2,gg3,R1,R2,Z1,Z2,A21,A22,M11,M22,M12,a,b,c,G11,G22,G12,T11,T12,T22,Wd1,Wd2,Hd1,Hd2)
return(mse)
}
