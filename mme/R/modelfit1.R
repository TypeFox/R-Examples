
####################################################################
#' Initial values for the variance components for Model 1
#'
#' This function is used in \code{\link[mme]{initial.values}} to calculate the initial values for the variance
#' components in the multinomial mixed model with one independent random effect in each category
#' of the response variable (Model 1).
#'
#' @param beta.0 initial values for the fixed effects obtained in \code{\link[mme]{initial.values}}.
#' @param y matrix with the response variable obtained from \code{\link[mme]{data.mme}}. The rows are the domains and the columns are the categories of the response variable.
#' @param Xk list of matrices with the auxiliary variables per category obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of domains.
#' @param M vector with the sample size of the areas.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#'\code{\link[mme]{wmatrix}}, \code{\link[mme]{prmu}},
#' \code{\link[mme]{Fbetaf}}, \code{\link[mme]{phi.direct}},
#' \code{\link[mme]{sPhikf}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit1}}, \code{\link[mme]{msef}},
#' \code{\link[mme]{mseb}}.
#' @return phi.0 vector of inicial values for the variance components
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Multinomial-based small area estimation of labour force indicators.
#' Statistical Modelling, 13, 153-178.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' data(simdata) #data
#' mod=1 #type of model
#' datar=data.mme(simdata,k,pp,mod)
#' ###beta values
#' beta.new=list()
#' beta.new[[1]]=matrix(c( 1.3,-1),2,1)
#' beta.new[[2]]=matrix(c( -1.6,1),2,1)
#'
#' ##Initial variance components
#' phi=phi.mult(beta.new,datar$y,datar$Xk,datar$n)

phi.mult<-function(beta.0,y,Xk,M){
D=nrow(y)
k=ncol(y)
u.old=rep(0,D*(k-1))
dim(u.old)=c(D,(k-1))
theta_ol=matrix(0,D,(k-1))
prmul=prmu(M,Xk,beta.0,u.old)
theta=prmul[[3]]
for (i in 1:(k-1)){
	for (j in 1:D){
		if ((y[j,i]==0) |(y[j,k]==0)) {theta_ol[j,i]=log((y[j,i]+1)/(y[j,k]+1))} else {theta_ol[j,i]=log(y[j,i]/y[j,k])}}}

dif=(theta-theta_ol)^2
phi.new=colSums(dif)/(D-1)
return(phi.0=phi.new)
}



####################################################################
#' Estimated mean and probabilities for Model 1
#'
#' This function calculates the estimated probabilities and the estimated mean
#' of the response variable, in the multinomial mixed model with one independent random effect in each category
#' of the response variable (Model 1).
#'
#' @param M vector with the area sample sizes.
#' @param Xk list of matrices with the auxiliary variables per category obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of domains.
#' @param beta fixed effects obtained from \code{\link[mme]{modelfit1}}.
#' @param u values of random effects obtained from \code{\link[mme]{modelfit1}}.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#'\code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult}},
#' \code{\link[mme]{Fbetaf}}, \code{\link[mme]{phi.direct}},
#' \code{\link[mme]{sPhikf}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit1}}, \code{\link[mme]{msef}},
#' \code{\link[mme]{mseb}}.
#####################
#' @return A list containing the following components:
#' \item{Estimated.probabilities}{matrix with the estimated probabilities
#' for the categories of response variable.}
#' \item{mean}{ matrix with the estimated mean of the response variable.}
#' \item{eta}{matrix with the estimated log-rates of the probabilities of each category over the reference category.}
###########################
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Multinomial-based small area estimation of labour force indicators.
#' Statistical Modelling, 13, 153-178.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' data(simdata) #data
#' mod=1 #type of model
#' D=nrow(simdata)
#' datar=data.mme(simdata,k,pp,mod)
#' initial=datar$initial
#'
#' ##Estimated mean and probabilities
#' mean=prmu(datar$n,datar$Xk,initial$beta.0,initial$u.0)


prmu<-function(M,Xk,beta,u){
		D=nrow(u)
		k=ncol(u)+1
		pr=matrix(1,D,k)
		theta=matrix(0,D,(k-1))
		mu=matrix(0,D,k)
		salida=list()
		for (i in 1:(k-1)){
			theta[,i]=Xk[[i]]%*%beta[[i]]+u[,i]}
		suma=rowSums(exp(theta))
		pr[,k]=(1+suma)^(-1)
		for (i in 1:(k-1)){
			pr[,i]=pr[,k]*exp(theta[,i])}
		for (j in 1:D){
			mu[j,]=M[j]*pr[j,]}
		mu=mu[,1:k-1]
		est=list(estimated.probabilities=pr,mean=mu,eta=theta)
		return(est)}



####################################################################
#' Inverse of the Fisher information matrix of the fixed and random effects in Model 1
#'
#' This function calculates the inverse of the Fisher information
#' matrix of the fixed effects (beta) and the random effects (u) and the score vectors S.beta and S.u, for the model with
#' one independent random effect in each category
#' of the response variable (Model 1). \code{\link[mme]{modelfit1}} uses the output of this function
#' to estimate the fixed and random effects by the PQL method.
#'
#' @param sigmap a list with the model variance-covariance matrices for each domain.
#' @param X list of matrices with the auxiliary variables obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of categories of the response variable.
#' @param Z design matrix  of random effects.
#' @param phi vector with the values of the variance components obtained from \code{\link[mme]{modelfit1}}.
#' @param y matrix with the response variable except the reference category. The rows are the domains and the columns are the categories of the  
#' response variable minus one.
#' @param mu matrix with the estimated mean of the response variable obtained from \code{\link[mme]{prmu}}.
#' @param u matrix with the values of random effects obtained from \code{\link[mme]{modelfit1}}.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#' \code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult}},
#' \code{\link[mme]{prmu}}, \code{\link[mme]{phi.direct}},
#' \code{\link[mme]{sPhikf}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit1}}, \code{\link[mme]{msef}},
#' \code{\link[mme]{mseb}}.
#' @return A list containing the following components.
#' \item{F.beta.beta}{the first diagonal element of the inverse of the
#' Fisher information matrix.}
#' \item{F.beta.u}{the element out of the diagonal of the inverse of the
#' Fisher information matrix.}
#' \item{F.u.u}{the second diagonal element of the inverse of the Fisher
#' information matrix.}
#' \item{S.beta}{beta scores.}
#' \item{S.u}{u scores.}
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Multinomial-based small area estimation of labour force indicators.
#' Statistical Modelling, 13 ,153-178.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' data(simdata) #data
#' mod=1 #type of model
#' datar=data.mme(simdata,k,pp,mod)
#' initial=datar$initial
#' mean=prmu(datar$n,datar$Xk,initial$beta.0,initial$u.0)
#' sigmap=wmatrix(datar$n,mean$estimated.probabilities)
#'
#' #Inverse of the Fisher information matrix
#' Fisher=Fbetaf(sigmap,datar$X,datar$Z,initial$phi.0,datar$y[,1:(k-1)],
#'        mean$mean,initial$u.0)

Fbetaf<- function(sigmap,X,Z,phi,y,mu,u){
D=nrow(mu)
A_d=list()
B_d=list()
C_d=list()
S_beta_d=list()
k=ncol(mu)+1
S_u_d=rep(0,D*(k-1))
dim(S_u_d)=c(D,(k-1))
salida=list()
for(i in 1:D){
	aa=t(X[[i]])
	bb=data.matrix((y[i,]-mu[i,]))
	S_beta_d[[i]]=aa%*%(bb)
	S_u_d[i,]=as.matrix((y[i,]-mu[i,]))
}
S.beta=add(S_beta_d)
S.u=S_u_d-(u/(matrix(rep(phi,D),D,(k-1),byrow=TRUE)))


for(i in 1:D){
	A_d[[i]]=crossprod((X[[i]]),sigmap[[i]])%*%X[[i]]
	B_d[[i]]=crossprod((X[[i]]),sigmap[[i]])
	C_d[[i]]=crossprod((Z[[i]]),solve(sigmap[[i]]+diag(1/phi)))
	}
A=add(A_d[1:D])
B=do.call(cbind,B_d)
C=do.call(cbind,C_d)

Fbb=solve(A-B%*%C%*%t(B))
Fbu=-Fbb%*%B%*%C
Fuu=C+C%*%crossprod((B),Fbb)%*%(B)%*%C
rm(A_d,B_d,C_d,A,B,C,S_beta_d,S_u_d)
fisher=list(F.beta.beta=Fbb,F.beta.u=Fbu,F.u.u=Fuu,S.beta=S.beta,S.u=S.u)
return(fisher)}

#######################
#####################


####################################################################
#' Fisher information matrix and score vectors of the variance components for Model 1
#'
#' This function computes the Fisher information matrix and the score vectors
#' of the variance components, for the multinomial mixed model with
#' one independent random effect in each category
#' of the response variable (Model 1). These values are used in the fitting algorithm
#' implemented in \code{\link[mme]{modelfit1}} to estimate the random effects. The algorithm adatps the ideas of Schall (1991) to a multivariate
#' model. The variance components are estimated by the REML method.
#'
#' @param pp vector with the number of the auxiliary variables per category.
#' @param sigmap a list with the model variance-covariance matrices for each domain obtained from \code{\link[mme]{wmatrix}}.
#' @param X list of matrices with the auxiliary variables obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of categories of the response variable minus one.
#' @param eta matrix with the estimated log-rates of probabilities of each category over the reference category obtained from  
#' \code{\link[mme]{prmu}}.
#' @param phi vector with the values of the variance components obtained from \code{\link[mme]{modelfit1}}.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#' \code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult}},
#' \code{\link[mme]{prmu}}, \code{\link[mme]{phi.direct}},
#' \code{\link[mme]{Fbetaf}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit1}}, \code{\link[mme]{msef}},
#' \code{\link[mme]{mseb}}.
#' @return A list containing the following components.
#' \item{S.k}{phi score vector.}
#' \item{F}{Fisher information matrix of the variance component phi.}
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Multinomial-based small area estimation of labour force indicators.
#' Statistical Modelling, 13 ,153-178.
#' @references Schall, R (1991). Estimation in generalized linear models with
#' random effects. Biometrika,
#' 78,719-727.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' data(simdata) #data
#' mod=1 #type of model
#' datar=data.mme(simdata,k,pp, mod)
#' initial=datar$initial
#' mean=prmu(datar$n,datar$Xk,initial$beta.0,initial$u.0)
#' sigmap=wmatrix(datar$n,mean$estimated.probabilities)
#'
#' ##Fisher information matrix and score vectors
#' Fisher.phi=sPhikf(pp,sigmap,datar$X,mean$eta,initial$phi.0)

sPhikf<- function(pp,sigmap,X,eta,phi){
Vd=list()
D=nrow(eta)
k=ncol(eta)+1
sk1=matrix(0,(k-1),1)
sk2_1=matrix(0,(k-1),1)
sk2_2=list()
sk2_3=list()
sk2_4=list()
sk2_5=list()
sk2=matrix(0,(k-1),1)
F=matrix(0,(k-1),(k-1))
s=rep(0,((k-1)*(k-1)))
dim(s)=c((k-1),(k-1))
salida=list()
rr=sum(pp+1)
qq=matrix(0,rr,rr)

for (i in 1:D){
	Vd[[i]]=solve(diag(phi)+solve(sigmap[[i]]))
	qq=qq+crossprod((X[[i]]),(Vd[[i]]))%*%X[[i]]}
qq=solve(qq)
for (j in 1:(k-1)){
	d=matrix(0,1,rr)
	a=matrix(0,rr,1)
	b=matrix(0,1,rr)
	c=matrix(0,rr,rr)
	sk1[j,1]=0
	sk2_1[j,1]=0
	for (i in 1:D){
		delta=matrix(0,(k-1),(k-1))
		for (l in 1:(k-1)){
			if (l==j) {delta[l,l]=1}}
		s=(Vd[[i]]-Vd[[i]]%*%X[[i]]%*%qq%*%crossprod((X[[i]]),Vd[[i]]))%*%delta
		sk1[j,1]=sk1[j,1]+sum(diag(s))
		sk2_1[j,1]=sk2_1[j,1]+t(eta[i,])%*%Vd[[i]]%*%delta%*%Vd[[i]]%*%eta[i,]
		d=d+t(eta[i,])%*%Vd[[i]]%*%delta%*%Vd[[i]]%*%X[[i]]
		a=a+crossprod((X[[i]]),Vd[[i]])%*%eta[i,]
		b=b+t(eta[i,])%*%Vd[[i]]%*%X[[i]]
		c=c+crossprod((X[[i]]),Vd[[i]])%*%delta%*%Vd[[i]]%*%X[[i]]
			}
		sk2_2[[j]]=d
		sk2_3[[j]]=a
		sk2_4[[j]]=b
		sk2_5[[j]]=c
		}
for (j in 1:(k-1)){
	sk2[j,1]=sk2_1[j,1]-2*sk2_2[[j]]%*%qq%*%sk2_3[[j]]+sk2_4[[j]]%*%qq%*%sk2_5[[j]]%*%qq%*%sk2_3[[j]]}

sk=(-0.5)*sk1+(0.5)*sk2
for (j in 1:(k-1)){
			for (l in 1:(k-1)){
				a=0
				b=0
				c=0
				for (i in 1:D){
					delta=matrix(0,(k-1),(k-1))
					for (ll in 1:(k-1)){
						if (ll==j) {delta[ll,ll]=1}}
					delta1=matrix(0,(k-1),(k-1))
					for (ll in 1:(k-1)){
						if (ll==l) {delta1[ll,ll]=1}}
					a=a+(0.5)*sum(diag(Vd[[i]]%*%delta%*%Vd[[i]]%*%delta1))
					b=b-sum(diag(Vd[[i]]%*%X[[i]]%*%qq%*%crossprod((X[[i]]),Vd[[i]])%*%delta%*%Vd[[i]]%*%delta1))
					c=c+(0.5)*sum(diag(Vd[[i]]%*%X[[i]]%*%qq%*%sk2_5[[j]]%*%qq%*%crossprod((X[[i]]),Vd[[i]])%*%delta1))
				}

				F[j,l]=a+b+c}}

score.F=list(S.k=sk,F=F)
rm(sk1,sk2_1,sk2_2,sk2_3,sk2_4,sk2_5,qq,sk2,a,b,c,Vd)
return(score.F)
}



####################################################################
#' Variance components for Model 1
#'
#' This function calculates the variance components for the multinomial mixed model with
#' one independent random effect in each category
#' of the response variable (Model 1). These values are used in the second part of the fitting algorithm
#' implemented in \code{\link[mme]{modelfit1}}. The algorithm adapts the ideas of Schall (1991) to a
#' multivariate model and the variance components are estimated by the REML method.
#'
#' @param sigmap a list with the model variance-covariance matrices for each domain obtained from \code{\link[mme]{wmatrix}}.
#' @param X list of matrices with the auxiliary variables obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of categories of the response variable minus one.
#' @param phi vector with the initial values of the variance components obtained from \code{\link[mme]{modelfit1}}.
#' @param u matrix with the values of the random effects obtained from \code{\link[mme]{modelfit1}}.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#' \code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult}},
#' \code{\link[mme]{prmu}}, \code{\link[mme]{Fbetaf}},
#' \code{\link[mme]{sPhikf}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit1}}, \code{\link[mme]{msef}},
#' \code{\link[mme]{mseb}}.
#' @return a list containing the following components.
#' \item{phi.new}{vector with the variance components.}
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Multinomial-based small area estimation of labour force indicators.
#' Statistical Modelling, 13 ,153-178.
#' @references Schall, R (1991). Estimation in generalized linear models with
#' random effects. Biometrika,
#' 78,719-727.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' data(simdata)   #data
#' mod=1 #type of model
#' datar=data.mme(simdata,k,pp,mod)
#' initial=datar$initial
#' mean=prmu(datar$n,datar$Xk,initial$beta.0,initial$u.0)
#' #model variance-covariance matrix
#' sigmap=wmatrix(datar$n,mean$estimated.probabilities)
#'
#' ##Variance components
#' phi=phi.direct(sigmap,initial$phi.0,datar$X,initial$u.0)


phi.direct<-function(sigmap,phi,X,u){
k=ncol(u)+1
D=nrow(u)
Trmll=list()
tau=rep(0,(k-1))
Vu_d=list()
V=list()
Vd=list()
AA=list()
for(i in 1:D){
	AA[[i]]=matrix(0,(k-1),(k-1)*D)
	for (j in 1:(k-1)){
		AA[[i]][j,((j-1)*D+i)]=1}
	}
ZZ=do.call(rbind,AA)
XX=do.call(rbind,X)

C=as(bdiag(sigmap),"matrix")

for (i in 1:(k-1)){
	Vu_d[[i]]=phi[i]*diag(D)}
Vu=as(bdiag(Vu_d),"matrix")

V=ZZ%*%Vu%*%t(ZZ)+solve(C)
T=solve(ZZ%*%C%*%t(ZZ)+solve(Vu))
Trml=T+T%*%t(ZZ)%*%C%*%XX%*%solve(t(XX)%*%solve(V)%*%XX)%*%t(XX)%*%C%*%ZZ%*%T
j=1
nr=nrow(Trml)
for (i in 1:(k-1)){
	Trmll[[i]]=T[j:(j+(nr/(k-1))-1),j:(j+(nr/(k-1))-1)]
	tau[i]=sum(diag(Trmll[[i]]))
	j=j+(nr/(k-1))
}

phi.new=diag((t(u)%*%u)/(D-tau))
rm(V,T,Trml,Vu,Trmll,tau,ZZ,AA,XX,Vu_d)
return(phi.new=phi.new)
}
#######################################################################

####################################################################
#' Function used to fit Model 1
#'
#' This function fits the multinomial mixed model with one independent random effect per category
#' of the response variable (Model 1), like in the formulation described in Lopez-Vizcaino et al. (2013).
#' The fitting algorithm combines the penalized quasi-likelihood method (PQL) for estimating
#' and predicting the fixed and random effects with the residual maximum likelihood method (REML)
#' for estimating the variance components. This function uses as initial values the output of the function
#' \code{\link[mme]{initial.values}}
#'
#' @param pp vector with the number of the auxiliary variables per category.
#' @param Xk list of matrices with the auxiliary variables per category obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of domains.
#' @param X list of matrices with the auxiliary variables obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of categories of the response variable minus one.
#' @param Z design matrix of random effects obtained from \code{\link[mme]{data.mme}}.
#' @param initial output of the function \code{\link[mme]{initial.values}}.
#' @param y matrix with the response variable except the reference category obtained from \code{\link[mme]{data.mme}}. The rows are the domains and the columns are the categories of the
#' response variable minus 1.
#' @param M vector with the area sample sizes.
#' @param MM vector with the population sample sizes.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#'\code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult}},
#' \code{\link[mme]{prmu}}, \code{\link[mme]{phi.direct}},
#' \code{\link[mme]{sPhikf}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{Fbetaf}}, \code{\link[mme]{msef}},
#' \code{\link[mme]{mseb}}.
#' @return A list containing the following components.
#' \item{Estimated.probabilities}{matrix with the estimated probabilities
#' for the categories of response variable.}
#' \item{Fisher.information.matrix.phi}{Fisher information matrix of the random effect.}
#' \item{Fisher.information.matrix.beta}{Fisher information matrix of the fixed effect.}
#' \item{u}{matrix with the estimated random effects.}
#' \item{mean}{matrix with the estimated mean of the response variable.}
#' \item{warning1}{0=OK,1=The model could not be fitted.}
#' \item{warning2}{0=OK,1=The value of the variance component is negative: the initial value
#' is taken.}
#' \item{beta.Stddev.p.value}{matrix with the estimated fixed effects, its standard
#' deviations and its p-values.}
#' \item{phi.Stddev.p.value}{matrix with the estimated variance components, its
#' standard deviations and its p-values.}
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Multinomial-based small area estimation of labour force indicators.
#' Statistical Modelling, 13, 153-178.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' data(simdata)  #data
#' mod=1 #type of model
#' datar=data.mme(simdata,k,pp,mod)
#'
#' #Model fit
#' result=modelfit1(pp,datar$Xk,datar$X,datar$Z,datar$initial,datar$y[,1:(k-1)],
#'       datar$n,datar$N)

modelfit1<-function(pp,Xk,X,Z,initial,y,M,MM){
maxiter1=100
maxiter2=100
D=nrow(y)
k=ncol(y)+1
eps=1e-4
aviso=0
cont=0
beta.new=initial[[1]]
phi.new=initial[[2]]
u.new=initial[[3]]
phi.prim=phi.new
ii=1
pasos2=0
mm=1
mmm=1
iter1=0
ll=0
resulbeta=list()
resulphi=list()
while (iter1<maxiter1) {
	phi.old=phi.new
	iter2=0
	pasos2=pasos2+1
	while(iter2<maxiter2 & mm>eps){
		beta.old=beta.new
		u.old=u.new
		prmul=prmu(M,Xk,beta.old,u.old)
		theta=prmul[[3]]
		pr=prmul[[1]]
		mu=prmul[[2]]
		sigmap=wmatrix(M,pr)
		for (i in 1:D){
			comp=det(sigmap[[i]])
			if (is.na(comp)){cont=1}
			if (is.na(comp)==FALSE & (abs(comp))<0.000001 ) {cont=1}}
		if (cont==0){
		y=data.matrix(y)
		inversaFbeta=Fbetaf(sigmap,X,Z,phi.old,y,mu,u.old)
		S.beta=inversaFbeta[[4]]
		S.u=inversaFbeta[[5]]
		Fbb=inversaFbeta[[1]]
		Fbu=inversaFbeta[[2]]
		Fuu=inversaFbeta[[3]]
		beta.new2=do.call(rbind,beta.new)
		beta.old2=do.call(rbind,beta.old)
		S.u2=matrix(t(S.u),((k-1)*D),1)
		u.new2=matrix(t(u.new),((k-1)*D),1)
		u.old2=matrix(t(u.old),((k-1)*D),1)
		beta.new2=beta.old2+Fbb%*%S.beta+Fbu%*%S.u2
		u.new2=u.old2+t(Fbu)%*%S.beta+Fuu%*%S.u2
		ucum=cumsum(pp+1)
		beta_new=list()
		beta.new[[1]]=as.matrix(beta.new2[(1:ucum[1]),1])
		for (i in 2:(k-1)){
			beta.new[[i]]=as.matrix(beta.new2[((ucum[i-1]+1):ucum[i]),1])
		}
		u.new=t(matrix((u.new2),(k-1),D))
		p1=(beta.new2-beta.old2)/beta.old2
		p2=(u.new-u.old)/u.old
		mm=max(abs(c(p1,p2)))
		if (ii==1) {
		beta.prim2=beta.new2
		beta.prim=beta.new
		Fbb.prim=Fbb
		u.prim=u.new
		}
		if (mm<eps){iter2=maxiter2} else {iter2=iter2+1}}
		if (cont==1){
		iter2=maxiter2
		valido=0}

	}

	if (cont==0){
	ii=ii+1
	prmul=prmu(M,Xk,beta.new,u.new)
	theta=prmul[[3]]
	pr=prmul[[1]]
	mu=prmul[[2]]
	sigmap=wmatrix(M,pr)
	phi.new=phi.direct(sigmap,phi.old,X,u.new)
	phi.new=as.vector(phi.new)
	p3=(phi.new-phi.old)/phi.old
	if (min(phi.new)<0.0001){
	aviso=1
	ll=ll+1
	phi.new=phi.prim
	#beta.new=beta.prim
	#beta.new2=beta.prim2
	#Fbb=Fbb.prim
	u.new=u.prim
	if (ll>1)
	{
	phi.new=phi.prim
	iter1=maxiter1}
	iter1=iter1+1
	}
	if (aviso==0)
	{mmm=max(abs(c(p1,p2,p3)))
	mm=1}
	else
	{mmm=max(abs(p3))
	mm=eps/10}
	if (mmm>eps){iter1=iter1+1} else {iter1=maxiter1}}
	if (cont==1){iter1=maxiter1}
}

if (cont==0){
prmul=prmu(MM,Xk,beta.new,u.new)
mu=prmul[[2]]
pr=prmul[[1]]
colnames(mu)=paste("Yest.",1:(k-1),sep="")
resul=list()
cii=ci(beta.new2,Fbb)
resulbeta=cbind(Beta=beta.new2,Std.dev=cii[[1]],p.value=cii[[2]])
colnames(resulbeta)=c("Estimate","Std.Error","p.value")
u=list()
for (i in 1:(k-1)){
#	u[[i]]=as.matrix((1:(pp[i]+1)))
	u[[i]]=as.matrix(colnames(Xk[[i]]))
	u[[i]][1]="Intercept"}
u=do.call(rbind,u)
rownames(resulbeta)=u
sk.F=sPhikf(pp,sigmap,X,theta,phi.new)
F=solve(sk.F[[2]])
cii=ci(phi.new,F)
resulphi=cbind(Estimate=phi.new,Std.Error=cii[[1]],p.value=cii[[2]])
}
result=list(Estimated.probabilities=pr,Fisher.information.matrix.phi=F,Fisher.information.matrix.beta=Fbb,u=u.new,mean=mu,warning1=cont,warning2=aviso,beta.Stddev.p.value=resulbeta,phi.Stddev.p.value=resulphi)
class(result)="mme"
return(result)
}
 ####################################################################
#' Analytic MSE for Model 1
#'
#' This function calculates the analytic MSE for the multinomial mixed model with one independent random effect per category
#' of the response variable (Model 1). See Lopez-Vizcaino et al. (2013), section 4, for details. The formulas
#' of Prasad and Rao (1990) are adapted to Model 1. This function uses the output of \code{\link[mme]{modelfit1}}.
#'
#' @param resul the output of the function \code{\link[mme]{modelfit1}}.
#' @param X list of matrices with the auxiliary variables obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of categories of the response variable minus one.
#' @param Z design matrix of random effects obtained from \code{\link[mme]{data.mme}}.
#' @param pp vector with the number of the auxiliary variables per category.
#' @param M vector with the area sample sizes.
#' @param MM vector with the population sample sizes.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#'\code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult}},
#' \code{\link[mme]{prmu}}, \code{\link[mme]{phi.direct}},
#' \code{\link[mme]{sPhikf}}, \code{\link[mme]{modelfit1}},
#' \code{\link[mme]{Fbetaf}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{mseb}}.
#' @return mse is a matrix with the MSE estimator calculated by adapting the explicit
#' formulas of Prasad and  Rao (1990). 
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Multinomial-based small area estimation of labour force indicators.
#' Statistical Modelling, 13, 153-178.
#' @references Prasad, NGN, Rao, JNK (1990).The estimation of the mean squared error of small
#' area estimators. Journal of the American Statistical Association, 85, 163-171.
#' @examples
#' require(Matrix)
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' data(simdata)  #data
#' mod=1 # type of model
#' datar=data.mme(simdata,k,pp,mod)
#' # Model fit
#' result=modelfit1(pp,datar$Xk,datar$X,datar$Z,datar$initial,datar$y[,1:(k-1)],
#'        datar$n,datar$N)
#'
#' #Analytic MSE
#' mse=msef(pp,datar$X,datar$Z,result,datar$N,datar$n)
#'


msef<-function(pp,X,Z, resul,MM,M){
phi.new=resul[[9]][,1]
F=resul[[2]]
pr=resul[[1]]
k=ncol(pr)
D=nrow(pr)
L=list()
F=solve(F)
Vu_d=list()
C_d=list()
C_dMM=list()
V=list()
Vk=list()
Vd=list()
W=list()
rr=sum(pp+1)
qq=matrix(0,rr,rr)
g2=list()
g1=list()
g3=list()
mse.analitic=matrix(0,D,(k-1))
ZZ=do.call(rbind,Z)
XX=do.call(rbind,X)
sigmap=wmatrix(M,pr)
sigmapMM=wmatrix(MM,pr)
for(i in 1:D){
        Vu_d[[i]]=t(Z[[i]])%*%diag(phi.new)
        C_d[[i]]=t(Z[[i]])%*%sigmap[[i]]
        C_dMM[[i]]=t(Z[[i]])%*%sigmapMM[[i]]}

Vu=do.call(cbind,Vu_d)
C=do.call(cbind,C_d)
CMM=do.call(cbind,C_dMM)
V=ZZ%*%Vu%*%t(ZZ)+solve(C)
for (i in 1:D){
        Vd[[i]]=solve(diag(phi.new)+solve(sigmap[[i]]))
        qq=qq+t(X[[i]])%*%solve(Vd[[i]])%*%(X[[i]])}
qq=solve(qq)
T=Vu-Vu%*%solve(V)%*%Vu
j=1
for (i in 1:D){
        g1[[i]]=sigmapMM[[i]]%*%T[j:(j+k-2),j:(j+k-2)]%*%sigmapMM[[i]]
        j=j+k-1}
R=Vu%*%solve(V)
j=1
for (i in 1:D){

        A21=sigmapMM[[i]]%*%X[[i]]
        A22=sigmapMM[[i]]%*%T[j:(j+k-2),j:(j+k-2)]%*%sigmap[[i]]%*%X[[i]]
        g2[[i]]=(A21-A22)%*%qq%*%(t(A21)-t(A22))
        j=j+k-1}

ww=matrix(0,(k-1)*D,(k-1)*D)
for (i in 1:(k-1)){
        delta=matrix(0,(k-1),(k-1))
        for (l in 1:(k-1)){
                if (l==i) {delta[l,l]=1}}
                tt=1
                for (j in 1:D){
                        ww[(tt:(tt+k-2)),(tt:(tt+k-2))]=delta
                        tt=tt+k-1}
                W[[i]]=ww
                }
for (i in 1:(k-1)){
                L[[i]]=(diag(1,(k-1)*D,(k-1)*D)-R)%*%W[[i]]%*%solve(V)}
for (kk in 1:D){
g3[[kk]]=0
for (i in 1:(k-1)){
        for (j in 1:(k-1)){
                g33=Z[[kk]]%*%CMM%*%L[[i]]%*%V%*%t(L[[j]])%*%t(CMM)%*%t(Z[[kk]])
                g3[[kk]]=g3[[kk]]+F[i,j]*g33
                }}}
g=list()
for (i in 1:D){
        g[[i]]=(g1[[i]]+2*g3[[i]]+g2[[i]])}
for (i in 1:D){
        mse.analitic[i,]=diag(g[[i]])
        }
colnames(mse.analitic)=paste("mse.",1:(k-1),sep="")
return(list(mse=mse.analitic))
rm(XX,ZZ,g1,g2,g3,A21,A22,L,W,ww,T,Vu,V,qq,Vu_d,C_d,C_dMM,R,C,CMM)
}
