
####################################################################
#' Initial values for the variance components in Model 2
#'
#' This function is used in \code{\link[mme]{initial.values}} to calculate the initial values for the variance
#' components in the multinomial mixed model with two independent random effects
#' for each category of the response variable: one domain random effect (u1) and another independent time and domain random effect (u2) (Model 2).
#'
#'
#' @param beta.0 initial values for the fixed effects obtained from \code{\link[mme]{initial.values}}.
#' @param y matrix with the response variable obtained from \code{\link[mme]{data.mme}}. The rows are the domains and the columns are the categories of the response variable.
#' @param Xk list of matrices with the auxiliary variables per category obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of domains.
#' @param M vector with the sample size of the areas.
#' @param u1 vector with the initial values for the first random effect obtained from \code{\link[mme]{initial.values}}.
#' @param u2 vector with the initial values for the second random effect obtained from \code{\link[mme]{initial.values}}.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#'\code{\link[mme]{wmatrix}}, \code{\link[mme]{prmu.time}},
#' \code{\link[mme]{Fbetaf.it}}, \code{\link[mme]{phi.direct.it}},
#' \code{\link[mme]{sPhikf.it}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit2}}, \code{\link[mme]{msef.it}},
#' \code{\link[mme]{mseb}}.
#' @return phi.0 vector of the initial values for the variance components.
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Multinomial-based small area estimation of labour force indicators.
#' Statistical Modelling, 13, 153-178.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' data(simdata2)  #data
#' mod=2 #Type of model
#' datar=data.mme(simdata2,k,pp,mod)
#' D=nrow(simdata2)
#' ###fixed effects values
#' beta.new=list()
#' beta.new[[1]]=matrix(c( 1.3,-1),2,1)
#' beta.new[[2]]=matrix(c( -1.6,1),2,1)
#' ## random effects values
#' u1.new=rep(0.01,((k-1)*datar$d))
#' dim(u1.new)=c(datar$d,k-1)
#' u2.new=rep(0.01,((k-1)*D))
#' dim(u2.new)=c(D,k-1)
#'
#' ##Initial variance components
#' phi=phi.mult.it(beta.new,datar$y,datar$Xk,datar$n,u1.new,u2.new)

phi.mult.it<-function(beta.0,y,Xk,M,u1,u2){
D=nrow(y)
k=ncol(y)
d=nrow(u1)
t=D/d
theta.ol=matrix(0,D,(k-1))
prmul=prmu.time(M,Xk,beta.0,u1,u2)
theta=prmul[[3]]
for (i in 1:(k-1)){
	for (j in 1:D){
		if ((y[j,i]==0) |(y[j,k]==0)) {theta.ol[j,i]=log(y[j,i]+1/(y[j,k]+1))} else {theta.ol[j,i]=log(y[j,i]/y[j,k])}}}

dif=(theta-theta.ol)^2
phi1.new=colSums(dif)/(2*(D-1))
return(phi.0=phi1.new)
}


############################################################################
#' Initial values for fitting algorithm to estimate the fixed and random effects and the variance components
#'
#' This function sets the initial values. An iterative algorithm fits the multinomial mixed models
#' that requires initial values for the fixed effects, the random
#' effects and the variance components. This initial values are used in \code{\link[mme]{modelfit1}},
#' \code{\link[mme]{modelfit2}} and \code{\link[mme]{modelfit3}}.
#'
#' @param d number of areas.
#' @param pp  vector with the number of auxiliary variables per category.
#' @param datar output of function \code{\link[mme]{data.mme}}.
#' @param mod a number specifying the type of model: 1=multinomial mixed model with one independent random effect for each category of the response variable 
#' (Model 1), 2=multinomial mixed model with two independent random effects for each category of the response variable: one domain random effect and another independent time and domain random effect (Model 2) and
#'  3= multinomial mixed model with two independent random effects for each category of the response variable: one domain random effect and another correlated time and domain random effect (Model 3).
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{wmatrix}},
#'\code{\link[mme]{phi.mult.it}}, \code{\link[mme]{prmu.time}},
#' \code{\link[mme]{Fbetaf.it}}, \code{\link[mme]{phi.direct.it}},
#' \code{\link[mme]{sPhikf.it}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit2}}, \code{\link[mme]{msef.it}},
#' \code{\link[mme]{mseb}}
#' @return A list containing the following components, depending on the chosen model.
#' \item{beta.0}{a list with the initial values for the fixed effects beta per category.}
#' \item{phi.0}{vector with the initial values for the variance components phi of Model 1.}
#' \item{phi1.0}{vector with the initial values for the variance components phi1 of Model 2 or 3.}
#' \item{phi2.0}{vector with the initial values for the variance components phi2 of Model 2 or 3.}
#' \item{u}{matrix with the initial values for the random effect for Model 1.}
#' \item{u1.0}{matrix with the initial values for the first random effect for Model 2 or 3.}
#' \item{u2.0}{matrix with the initial values for the second random effect for Model 2 or 3.}
#' \item{rho.0}{vector with the initial values for the correlation parameter for Model 3.}
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Multinomial-based small area estimation of labour force indicators.
#' Statistical Modelling, 13, 153-178.
#' @references
#' Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013). Small area estimation of labour force indicators under a multinomial mixed model
#' with correlated time and area effects. Submitted for review.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' data(simdata)
#' D=nrow(simdata)
#' mod=1 #Type of model
#' datar=data.mme(simdata,k,pp,mod)
#'
#' ## Initial values for fixed, random effects and variance components
#' initial=datar$initial

initial.values<-function(d,pp,datar,mod){
Xk=datar$Xk
M=datar$n
resp=datar$y
k=ncol(resp)
D=nrow(resp)
beta=beta.new=list()
for (i in 1:(k-1)){
dep=cbind(round(resp[,i]),round(rowSums(resp)-resp[,i]))
mod1=glm(dep ~ (Xk[[i]])-1,family=binomial)
tt=coef(mod1)
beta[[i]]=as.matrix(tt)
beta.new[[i]]=(beta[[i]])
}
u=list()
for (i in 1:(k-1)){
	u[[i]]=as.matrix((1:(pp[i]+1)))
	rownames(beta.new[[i]])=u[[i]]}
if (mod==1){
phi.new=phi.mult(beta.new,datar$y,datar$Xk,datar$n)
phi.new=as.vector(phi.new)
u.new=rep(0.5,((k-1)*D))
dim(u.new)=c(D,k-1)
result=list(beta.0=beta.new,phi.0=phi.new,u.0=u.new)}
if (mod==2){
u1.new=rep(0.01,((k-1)*d))
dim(u1.new)=c(d,k-1)
u2.new=rep(0.01,((k-1)*D))
dim(u2.new)=c(D,k-1)
phi1.new=phi.mult.it(beta.new,resp,Xk,M,u1.new,u2.new)
phi2.new=phi1.new
result=list(beta.0=beta.new,phi1.0=phi1.new,phi2.0=phi2.new,u1.0=u1.new,u2.0=u2.new)}
if (mod==3) {
u1.new=rep(0.01,((k-1)*d))
dim(u1.new)=c(d,k-1)
u2.new=rep(0.01,((k-1)*D))
dim(u2.new)=c(D,k-1)
resul=phi.mult.ct(beta.new,resp,Xk,M,u1.new,u2.new)
phi1.new=resul[[1]]
phi2.new=phi1.new
rho=resul[[2]]
result=list(beta.0=beta.new,phi1.0=phi1.new,phi2.0=phi2.new,u1.0=u1.new,u2.0=u2.new,rho.0=rho) }
return(result)
}


####################################################################
#' Estimated mean and probabilities for Model 2 and 3
#'
#' This function calculates the estimated probabilities and the estimated mean
#' of the response variable, in the multinomial mixed models with two independent random effects, one random effect
#' associated with the area and the other associated with the time,
#' for each category of the response variable. The first model assumes independent time and domain random effect (Model 2)
#' and the second model assumes correlated time and domain random effect (Model 3).
#'
#' @param M vector with the area sample sizes.
#' @param Xk list of matrices with the auxiliary variables per category obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of domains.
#' @param beta a list with the values for the fixed effects beta per category obtained from \code{\link[mme]{modelfit2}}.
#' @param u1 a vector with the values of the first random effect obtained from \code{\link[mme]{modelfit2}} or \code{\link[mme]{modelfit3}}.
#' @param u2 a vector with the values of the second random effect obtained from \code{\link[mme]{modelfit2}} or \code{\link[mme]{modelfit3}}.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#'\code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult.it}},
#' \code{\link[mme]{Fbetaf.it}}, \code{\link[mme]{phi.direct.it}},
#' \code{\link[mme]{sPhikf.it}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit2}}, \code{\link[mme]{msef.it}},
#' \code{\link[mme]{mseb}}
#####################
#' @return A list containing the following components:
#' \item{Estimated.probabilities}{matrix with the estimated probabilities
#' for the categories of response variable.}
#' \item{mean}{matrix with the estimated mean of the response variable.}
#' \item{eta}{matrix with the estimated log-rates of the probabilities of each category over the reference category.}
###########################
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Small area estimation of labour force indicator under a multinomial mixed model
#' with correlated time and area effects. Submited for review.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' mod=2 #Type of model
#' data(simdata2) # data
#' datar=data.mme(simdata2,k,pp,mod)
#' initial=datar$initial
#'
#' ## Estimated mean and estimated probabilities
#' mean=prmu.time(datar$n,datar$Xk,initial$beta.0,initial$u1.0,initial$u2.0)

prmu.time<-function(M,Xk,beta,u1,u2){
		k=ncol(u1)+1
		d=nrow(u1)
		D=nrow(u2)
		t=D/d
		pr=matrix(1,D,k)
		theta=matrix(0,D,(k-1))
		mu=matrix(0,D,k)
		salida=list()
		u11=matrix(0,D,(k-1))
		jj=1
		for (i in 1:d){
			for(j in 1:t){
				u11[jj,]=u1[i,]
				jj=jj+1}}
		for (i in 1:(k-1)){
			theta[,i]=Xk[[i]]%*%beta[[i]]+u2[,i]+u11[,i]}
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
#' The inverse of the Fisher information matrix of the fixed and random effects for Model 2
#'
#' This function calculates the score vector S and the inverse of the Fisher information
#' matrix for the fixed (beta) and the random effects (u1, u2) in Model 2. This model has two independet sets of random effects. 
#' The first one contains independent random effects u1dk associated to each category and domain. The second set contains random effects
#' u2dkt associated to each category, domain and time period. Model 2 assumes that the u2dk are independent across time. 
#' \code{\link[mme]{modelfit2}} uses the output of this function to estimate the fixed and random effect by the PQL method.
#'
#' @param sigmap a list with the model variance-covariance matrices for each domain.
#' @param X list of matrices with the auxiliary variables obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of categories of the response variable minus one.
#' @param Z design matrix  of random effects obtained from \code{\link[mme]{data.mme}}.
#' @param phi1 vector with the first variance component obtained from \code{\link[mme]{modelfit2}}.
#' @param phi2 vector with the second variance component obtained from \code{\link[mme]{modelfit2}}.
#' @param y matrix with the response variable, except the reference category obtained from \code{\link[mme]{data.mme}}. The rows are the domains and the columns are the categories of the response variable minus one.
#' @param mu matrix with the estimated mean of the response variable obtained from \code{\link[mme]{prmu.time}}.
#' @param u1 matrix with the values of the first random effect obtained from \code{\link[mme]{modelfit2}}.
#' @param u2 matrix with the values of the second random effect obtained from \code{\link[mme]{modelfit2}}.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#'\code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult.it}},
#' \code{\link[mme]{prmu.time}}, \code{\link[mme]{phi.direct.it}},
#' \code{\link[mme]{sPhikf.it}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit2}}, \code{\link[mme]{msef.it}},
#' \code{\link[mme]{mseb}}.
#' @return A list containing the following components.
#' \item{F}{the inverse of the Fisher information matrix.}
#' \item{S}{(beta, u1, u2) scores}
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Small area estimation of labour force indicator under a multinomial mixed model
#' with correlated time and area effects. Submitted for review.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' mod=2  #Type of model
#' data(simdata2) #data
#' datar=data.mme(simdata2,k,pp,mod)
#' initial=datar$initial
#' mean=prmu.time(datar$n,datar$Xk,initial$beta.0,initial$u1.0,initial$u2.0)
#' sigmap=wmatrix(datar$n,mean$estimated.probabilities)
#'
#' ##The inverse of the Fisher information matrix of the fixed effects
#' Fisher=Fbetaf.it(sigmap,datar$X,datar$Z,initial$phi1.0,initial$phi2.0,
#'        datar$y[,1:(k-1)],mean$mean,initial$u1.0,initial$u2.0)

Fbetaf.it<- function(sigmap,X,Z,phi1,phi2,y,mu,u1,u2){
d=nrow(u1)
D=nrow(u2)
k=ncol(u1)+1
t=D/d
A_d=list()
B_d=list()
B2_d=list()
C_d=list()
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
S_u2=S_u2_d-(u2/(matrix(rep(phi2,D),D,(k-1),byrow=TRUE)))
j=1
for (i in 1:d){
		S_u1[i,]=colSums(S_u2_d[j:(j+t-1),])
		j=j+t}
S_u1=S_u1-(u1/(matrix(rep(phi1,d),d,(k-1),byrow=TRUE)))

for(i in 1:D){
	A_d[[i]]=crossprod((X[[i]]),sigmap[[i]])%*%X[[i]]
	B_d[[i]]=crossprod((X[[i]]),sigmap[[i]])
	D_d[[i]]=crossprod((Z[[i]]),(sigmap[[i]]))
}

A=as(add(A_d[1:D]),"sparseMatrix")
#Hu2beta=B#
B=as(do.call(cbind,B_d), "sparseMatrix")

B2_d=addtolist(B_d,t,d)
#Hu1beta=B2#
B2=as(do.call(cbind,B2_d),"sparseMatrix")
#Hu2u2=DD#
DD=as(do.call(cbind,D_d),"sparseMatrix")
D2_d=addtolist(D_d,t,d)
#Hu2u1=D2#
D2=as(do.call(cbind,D2_d),"sparseMatrix")
C2_d=addtolist(D_d,t,d)
C2=do.call(cbind,C2_d)
#Hu1u1=C22#
C22=addtomatrix(C2,d,t,k)
F11=as(C22,"sparseMatrix")+as(diag(rep((1/phi1),d)),"sparseMatrix")
F22=DD+as(diag(rep((1/phi2),D)),"sparseMatrix")
Fuu=as(rBind(cBind(F11,as(t(as(D2,"matrix")),"sparseMatrix")),cBind(D2,F22)),"sparseMatrix")
Fuui=as(solve(Fuu),"sparseMatrix")
int=as(solve(F22),"sparseMatrix")
Fuubb=rBind(as(t(as(B2,"matrix")),"sparseMatrix"),as(t(as(B,"matrix")),"sparseMatrix"))
Fbb=as(solve(A-t(as(Fuubb,"matrix"))%*%Fuui%*%(Fuubb)),"sparseMatrix")
Fub=-1*Fbb%*%t(as(Fuubb,"matrix"))%*%Fuui
Fu1u1=solve(F11-t(as(D2,"matrix"))%*%int%*%D2)
Fu1u2=-1*Fu1u1%*%t(as(D2,"matrix"))%*%int
Fu2u2=int+int%*%D2%*%Fu1u1%*%t(as(D2,"matrix"))%*%int
Fiuu=rBind(cBind(Fu1u1,Fu1u2),cBind(t(as(Fu1u2,"matrix")),Fu2u2))
Fuu=Fiuu+Fiuu%*%Fuubb%*%Fbb%*%t(as(Fuubb,"matrix"))%*%Fiuu
F=rBind(cBind(Fbb,(Fub)),cBind(t(as(Fub,"matrix")),Fuu))
S_u11=matrix(t(S_u1),d*(k-1),1)
S_u22=matrix(t(S_u2),D*(k-1),1)
S=rbind(S_beta,rbind(S_u11,S_u22))
rm(A,B,B2_d,B2,DD,D2_d,D2,C2_d,C2,C22, A_d,B_d,F11,F22,Fuu,Fuui,int,Fbb,Fub,Fu1u1,Fu1u2,Fu2u2,Fiuu)
F=as(F,"matrix")
fisher=list(F=F,S=S)
return(fisher)
}

####################################################################
#' Variance components for Model 2
#'
#' This function calculates the variance components for the multinomial mixed model with two independent random effects
#' for each category of the response variable: one domain random effect and another independent time and domain random effect (Model 2). This variance components
#' are used in the second part of the fitting algorithm
#' implemented in \code{\link[mme]{modelfit2}}. The algorithm adapts the ideas of Schall (1991) to a multivariate model. The variance components are
#' estimated by the REML method.
#'
#' @param pp vector with the number of auxiliary variables per category.
#' @param sigmap a list with the model variance-covariance matrices for each domain obtained from \code{\link[mme]{wmatrix}}.
#' @param X list of matrices with the auxiliary variables obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of categories of the response variable.
#' @param phi1 vector with the initial values of the first variance component obtained from \code{\link[mme]{modelfit2}}.
#' @param phi2 vector with the initial values of the second variance component obtained from \code{\link[mme]{modelfit2}}.
#' @param u1 matrix with the values of the first random effect obtained from \code{\link[mme]{modelfit2}}.
#' @param u2 matrix with the values of the second random effect obtained from \code{\link[mme]{modelfit2}}.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#' \code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult.it}},
#' \code{\link[mme]{prmu.time}}, \code{\link[mme]{Fbetaf.it}}
#' \code{\link[mme]{sPhikf.it}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit2}}, \code{\link[mme]{msef.it}},
#' \code{\link[mme]{mseb}}
#' @return a list containing the following components.
#' \item{phi1.new}{vector with the values of the variance component for the first random effect.}
#' \item{phi2.new}{vector with the values of the variance component for the second random effect.}
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
#' d=10 #number of areas
#' mod=2  #Type of model
#' data(simdata2)   #data
#' datar=data.mme(simdata2,k,pp,mod)
#' initial=datar$initial
#' mean=prmu.time(datar$n,datar$Xk,initial$beta.0,initial$u1.0,initial$u2.0)
#' sigmap=wmatrix(datar$n,mean$estimated.probabilities) #variance-covariance
#'
#' ## The variance components
#' phi.it=phi.direct.it(pp,sigmap,datar$X,initial$phi1.0,initial$phi2.0,
#'        initial$u1.0,initial$u2.0)


phi.direct.it<- function(pp,sigmap,X,phi1,phi2,u1,u2){
d=nrow(u1)
D=nrow(u2)
t=D/d
k=ncol(u1)+1
r=sum(pp+1)
Vd=list()
qq=matrix(0,r,r)
Ld=list()
sigmapt=list()

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
j=1
Iq1=matrix(rep(diag(rep(1,(k-1))),t),(t*(k-1)),(k-1),byrow=TRUE)
for (i in 1:d){
	phi22=diag(1/phi2)
	a=matrix(0,(k-1),(k-1))
	for (l in 1:t){
		sigmapt[[l]]=sigmap[[j]]-sigmap[[j]]%*%solve(sigmap[[j]]+phi22)%*%sigmap[[j]]
		a=a+sigmapt[[l]]
		j=j+1}
	Ld[[i]]=do.call("blockdiag",sigmapt)
	Vd[[i]]=Ld[[i]]-Ld[[i]]%*%Iq1%*%solve(diag(1/phi1)+a)%*%t(Iq1)%*%Ld[[i]]
	qq=qq+crossprod((Xd[[i]]),(Vd[[i]]))%*%Xd[[i]]
	}
qq=as(solve(qq),"sparseMatrix")

Xt=as(do.call(rbind,Xd),"sparseMatrix")
W=bdiag(sigmap)
V1=list()
V2=list()

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

for (i in 1:(k-1)){
	V1[[i]]=(1/phi1[i])*diag(1,d)
	V2[[i]]=(1/phi2[[i]])*diag(1,D)}

T1=as(solve(t(as(Z1,"matrix"))%*%W%*%Z1+bdiag(V1)),"sparseMatrix")
T2=as(solve(t(as(Z2,"matrix"))%*%W%*%Z2+bdiag(V2)),"sparseMatrix")

Treml1=T1+as(T1%*%t(as(Z1,"matrix"))%*%W%*%Xt%*%qq%*%t(as(Xt,"matrix"))%*%W%*%Z1%*%T1,"sparseMatrix")
Treml2=T2+as(T2%*%t(as(Z2,"matrix"))%*%W%*%Xt%*%qq%*%t(as(Xt,"matrix"))%*%W%*%Z2%*%T2,"sparseMatrix")
j=1
Treml1=as(Treml1,"matrix")
Treml2=as(Treml2,"matrix")
tau=list()
Trmll=list()
nr=nrow(Treml1)

for (i in 1:(k-1)){
	Trmll[[i]]=Treml1[j:(j+(nr/(k-1))-1),j:(j+(nr/(k-1))-1)]
	tau[i]=sum(diag(Trmll[[i]]))/phi1[i]
	j=j+(nr/(k-1))
}
tau=do.call(rbind,tau)
tau=as.vector(tau)
phi1.new=diag((t(u1)%*%u1)/diag(d-tau))

tau=list()
Trmll=list()
nr=nrow(Treml2)
j=1
for (i in 1:(k-1)){
	Trmll[[i]]=Treml2[j:(j+(nr/(k-1))-1),j:(j+(nr/(k-1))-1)]
	tau[i]=sum(diag(Trmll[[i]]))/phi2[i]
	j=j+(nr/(k-1))
}
tau=do.call(rbind,tau)
tau=as.vector(tau)
phi2.new=diag((t(u2)%*%u2)/diag((d*t)-tau))
rm(Trmll,tau,Treml2,Treml1, T1, T2,V1,V2, Z1, Z2,Xt,W,Vd,Ld,qq)
resul=list(phi1.new=phi1.new,phi2.new=phi2.new)
}



####################################################################
#' Fisher information matrix and score vectors of the variance components for Model 2
#'
#' This function computes the Fisher information matrix and the score vectors
#' of the variance components, for the multinomial mixed model with two independent random effects
#' for each category of the response variable: one domain random effect and another independent time and domain random effect (Model 2).
#' These values are used in the fitting algorithm implemented in \code{\link[mme]{modelfit2}} to estimate the random effects. The algorithm adatps the ideas of Schall (1991) to a multivariate
#' model. The variance components are estimated by the REML method.
#'
#' @param d number of areas.
#' @param t number of time periods.
#' @param pp vector with the number of the auxiliary variables per category.
#' @param sigmap a list with the model variance-covariance matrices for each domain obtained from \code{\link[mme]{wmatrix}}.
#' @param X list of matrices with the auxiliary variables obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of categories of the response variable minus one.
#' @param eta matrix with the estimated log-rates of probabilities of each category over the reference category obtained from \code{\link[mme]{prmu.time}}.
#' @param phi1 vector with the values of the first variance component obtained from \code{\link[mme]{modelfit2}}.
#' @param phi2 vector with the values of the second variance component obtained from \code{\link[mme]{modelfit2}}.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#' \code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult.it}},
#' \code{\link[mme]{prmu.time}}, \code{\link[mme]{phi.direct.it}},
#' \code{\link[mme]{Fbetaf.it}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{modelfit2}}, \code{\link[mme]{msef.it}},
#' \code{\link[mme]{mseb}}
#' @return A list containing the following components.
#' \item{S}{phi score vector.}
#' \item{F}{Fisher information matrix of the variance components phi1 and phi2.}
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
#' mod=2 #Type of model
#' data(simdata2) #data
#' datar=data.mme(simdata2,k,pp,mod)
#' initial=datar$initial
#' mean=prmu.time(datar$n,datar$Xk,initial$beta.0,initial$u1.0,initial$u2.0)
#' sigmap=wmatrix(datar$n,mean$estimated.probabilities)
#'
#' ##Fisher information matrix and score vectors
#' Fisher.phi=sPhikf.it(datar$d,datar$t,pp,sigmap,datar$X,mean$eta,initial$phi1.0,
#'            initial$phi2.0)

sPhikf.it<- function(d,t,pp,sigmap,X,eta,phi1,phi2){
k=ncol(eta)+1
r=sum(pp+1)
Vd=list()
F=matrix(0,(2*(k-1)),(2*(k-1)))
qq=matrix(0,r,r)
Ld=list()
sigmapt=list()

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
	thetaa[[i]]=matrix(t(eta[(j:(j+t-1)),]),((k-1)*t),1,byrow=TRUE)
		j=j+t}

j=1
Iq1=matrix(rep(diag(rep(1,(k-1))),t),(t*(k-1)),(k-1),byrow=TRUE)
for (i in 1:d){
	phi22=diag(1/phi2)
	a=matrix(0,(k-1),(k-1))
	for (l in 1:t){
		sigmapt[[l]]=sigmap[[j]]-sigmap[[j]]%*%solve(sigmap[[j]]+phi22)%*%sigmap[[j]]
		a=a+sigmapt[[l]]
		j=j+1}
	Ld[[i]]=do.call("blockdiag",sigmapt)
	Vd[[i]]=Ld[[i]]-Ld[[i]]%*%Iq1%*%solve(diag(1/phi1)+a)%*%t(Iq1)%*%Ld[[i]]
	qq=qq+crossprod((Xd[[i]]),(Vd[[i]]))%*%Xd[[i]]}

qq=solve(qq)

S1=rep(0,2*(k-1))
S2=rep(0,2*(k-1))
delta11=list()
sum1=list()
sum2=list()
sum3=list()
sum4=list()
for (i in 1:(2*(k-1))){
	sum1[[i]]=matrix(0,1,r)
	sum2[[i]]=matrix(0,r,1)
	sum3[[i]]=matrix(0,1,r)
	sum4[[i]]=matrix(0,r,r)}

o=1
for(i in 1:2){
	for (j in 1:(k-1)){
			delta=matrix(0,(k-1),(k-1))
			for (ll in 1:(k-1)){
				if (ll==j) {delta[ll,ll]=1}}
			for (ll in 1:t){
				delta11[[ll]]=delta
				}
			if (i==1) {Vakd=Iq1%*%delta%*%t(Iq1)}
			if (i==2) {Vakd=do.call("blockdiag",delta11)}
			for (l in 1:d){
			S1[o]=S1[o]+sum(diag((Vd[[l]]-Vd[[l]]%*%Xd[[l]]%*%qq%*%t(Xd[[l]])%*%Vd[[l]])%*%Vakd))
			S2[o]=S2[o]+t(thetaa[[l]])%*%Vd[[l]]%*%Vakd%*%(Vd[[l]])%*%thetaa[[l]]
			sum1[[o]]=sum1[[o]]+t(thetaa[[l]])%*%Vd[[l]]%*%Vakd%*%(Vd[[l]])%*%Xd[[l]]
			sum2[[o]]=sum2[[o]]+t(Xd[[l]])%*%Vd[[l]]%*%thetaa[[l]]
			sum3[[o]]=sum3[[o]]+t(thetaa[[l]])%*%Vd[[l]]%*%Xd[[l]]
			sum4[[o]]=sum4[[o]]+t(Xd[[l]])%*%Vd[[l]]%*%Vakd%*%(Vd[[l]])%*%Xd[[l]]}
			o=o+1}}
tt=rep(0,(2*(k-1)))
for(i in 1:(2*(k-1))){
	tt[i]=-2*(sum1[[i]]%*%qq%*%sum2[[i]])+sum3[[i]]%*%qq%*%sum4[[i]]%*%qq%*%sum2[[i]]}

S2=S2+tt
S=-0.5*S1+0.5*S2

F1=rep(0,((k-1)*(k-1)*4))
F2=rep(0,((k-1)*(k-1)*4))
F3=rep(0,((k-1)*(k-1)*4))
aa=1
o=1
for (i in 1:2){
		for(ii in 1:(k-1)){
			for(j in 1:2){
			for (jj in 1:(k-1)){
				delta=matrix(0,(k-1),(k-1))
				for (ll in 1:(k-1)){
						if (ll==ii) {delta[ll,ll]=1}}
				for (ll in 1:t){
						delta11[[ll]]=delta}

				if (i==1) {Vakd1=Iq1%*%delta%*%t(Iq1)}
				if (i==2) {Vakd1=do.call("blockdiag",delta11)}

				delta=matrix(0,(k-1),(k-1))
				for (ll in 1:(k-1)){
						if (ll==jj) {delta[ll,ll]=1}}
				for (ll in 1:t){
						delta11[[ll]]=delta
				}
				if (j==1) {Vakd2=Iq1%*%delta%*%t(Iq1)}
				if (j==2) {Vakd2=do.call("blockdiag",delta11)}

				for (l in 1:d){
					F1[aa]=F1[aa]+0.5*sum(diag(Vd[[l]]%*%Vakd1%*%Vd[[l]]%*%Vakd2))
					F2[aa]=F2[aa]+0.5*sum(diag(Vd[[l]]%*%Vakd1%*%Vd[[l]]%*%Xd[[l]]%*%qq%*%t(Xd[[l]])%*%Vd[[l]]%*%Vakd2))
					F3[aa]=F3[aa]+0.5*sum(diag(Vd[[l]]%*%Xd[[l]]%*%qq%*%sum4[[o]]%*%qq%*%t(Xd[[l]])%*%Vd[[l]]%*%Vakd2))}
					aa=aa+1
				}}
					o=o+1}}

F=F1+2*F2+F3

FF=list()
FFF=list()
o=1
for (i in 1:2){
	for (j in 1:2){
		FF[[j]]=matrix(F[o:(o+(k-1)*(k-1)-1)],(k-1),(k-1))
		o=o+(k-1)*(k-1)}
		FFF[[i]]=do.call(cbind,FF)}
F=do.call(rbind,FFF)

F=solve(F)

#F=solve(matrix(F,4,4,byrow=TRUE))
S=matrix(t(S),4,1)
A=F%*%S

rm(sum1,sum2,sum3,sum4, Vakd1,Vakd2, Vakd, S1,S2,Xd,Vd,Ld,sigmapt,thetaa, Iq1)

score.F=list(S=S,F=F)
return(score.F)
}

####################################################################
#' Function to fit Model 2
#'
#' This function fits the multinomial mixed model with two independent random effects
#' for each category of the response variable: one domain random effect and another independent time and domain random effect (Model 2). The formulation is described in Lopez-Vizcaino et al. (2013).
#' The fitting algorithm combines the penalized quasi-likelihood method (PQL) for estimating
#' and predicting the fixed and random effects, respectively, with the residual maximum likelihood method (REML)
#' for estimating the variance components. This function uses as initial values the output of the function
#' \code{\link[mme]{initial.values}}.
#'
#' @param d number of areas.
#' @param t number of time periods.
#' @param pp vector with the number of the auxiliary variables per category.
#' @param Xk list of matrices with the auxiliary variables per category obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of domains.
#' @param X list of matrices with the auxiliary variables obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of categories of the response variable minus one.
#' @param Z design matrix  of random effects \code{\link[mme]{data.mme}}.
#' @param initial output of the function \code{\link[mme]{initial.values}}.
#' @param y matrix with the response variable obtained from \code{\link[mme]{data.mme}}, except the reference category. The rows are the domains and the columns are the categories of the response variable minus one.
#' @param M vector with the area sample sizes.
#' @param MM vector with the population sample sizes.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#'\code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult.it}},
#' \code{\link[mme]{prmu.time}}, \code{\link[mme]{phi.direct.it}},
#' \code{\link[mme]{sPhikf.it}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{Fbetaf.it}}, \code{\link[mme]{msef.it}},
#' \code{\link[mme]{mseb}}
#' @return A list containing the following components.
#' \item{Estimated.probabilities}{matrix with the estimated probabilities
#' for the categories of response variable.}
#' \item{Fisher.information.matrix.phi}{Fisher information matrix of the variance components.}
#' \item{Fisher.information.matrix.beta}{Fisher information matrix of the fixed effects. }
#' \item{u1}{matrix with the estimated first random effect.}
#' \item{u2}{matrix with the estimated second random effect.}
#' \item{mean}{matrix with the estimated mean of response variable.}
#' \item{warning1}{0=OK,1=The model could not be fitted.}
#' \item{warning2}{0=OK,1=The value of the variance component is negative: the initial value
#' is taken.}
#' \item{beta.Stddev.p.value}{matrix with the estimated fixed effects, its standard
#' deviations and its p-values.}
#' \item{phi.Stddev.p.value}{matrix with the estimated variance components, its
#' standard deviations and its p-values.}
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Small area estimation of labour force indicators under a multinomial mixed model
#' with correlated time and area effects. Submitted for review.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' mod=2 #type of model
#' data(simdata2)  #data
#' datar=data.mme(simdata2,k,pp,mod)
#'
#' ##Model fit
#' result=modelfit2(datar$d,datar$t,pp,datar$Xk,datar$X,datar$Z,datar$initial,datar$y[,1:(k-1)],
#'        datar$n,datar$N)

modelfit2<-function(d,t,pp,Xk,X,Z,initial,y,M,MM){
D=d*t
beta.new=initial[[1]]
phi1.new=initial[[2]]
phi2.new=initial[[3]]
u1.new=initial[[4]]
u2.new=initial[[5]]
k=ncol(u1.new)+1
maxiter1=100
maxiter2=100
iter1=0
eps=1e-3
phi1.prim=phi1.new
phi2.prim=phi2.new
p=sum(pp+1)
cont=0
ii=1
ll=0
mm=1
aviso=0
Fk=matrix(0,(k-1),(k-1))
while (iter1<maxiter1) {
	phi1.old=phi1.new
	phi2.old=phi2.new
	iter2=0
	p1=matrix(1,p,(k-1))
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
		inversaFbeta=Fbetaf.it(sigmap,X,Z,phi1.old,phi2.old,y,mu,u1.old,u2.old)
		F=inversaFbeta[[1]]
		S=inversaFbeta[[2]]
		rm(inversaFbeta)
		beta.u.old=rbind(do.call(rbind,beta.old),rbind(matrix(t(u1.old),d*(k-1),1),matrix(t(u2.old),D*(k-1),1)))
		beta.u.new=beta.u.old+F%*%S

		p1=(beta.u.new-beta.u.old)/beta.u.old
		mm=max(abs(p1))
		if (mm<eps){iter2=maxiter2} else {iter2=iter2+1}
		ucum=cumsum(pp+1)
		beta_new=list()
		beta.new[[1]]=as.matrix(beta.u.new[(1:ucum[1]),1])
		for (i in 2:(k-1)){
			beta.new[[i]]=as.matrix(beta.u.new[((ucum[i-1]+1):ucum[i]),1])
		}

		u1.new=matrix(beta.u.new[(ucum[k-1]+1):(ucum[k-1]+d*(k-1)),],d,k-1,byrow=TRUE)
		u2.new=matrix(beta.u.new[(ucum[k-1]+1+d*(k-1)):nrow(beta.u.new),],D,k-1,byrow=TRUE)}
		if (ii==1) {
		beta.prim=beta.new
		F.prim=F
		u1.prim=u1.new
		u2.prim=u2.new
		}

		if (cont==1){iter2=maxiter2}

	}
	if (cont==0){

	prmul=prmu.time(M,Xk,beta.new,u1.new,u2.new)
	theta=prmul[[3]]
	pr=prmul[[1]]
	mu=prmul[[2]]
	rm(prmul)
	sigmap=wmatrix(M,pr)

	#sk_F=sPhikf.ti(d,t,pp,sigmap,X,Z,theta,phi1.old,phi2.old)

	#sk=sk_F[[1]]
	#Fk=sk_F[[2]]
	#if( ii==1) {Fk.prim=Fk}
	ii=ii+1
	phi.old=rbind(data.matrix(phi1.old),data.matrix(phi2.old))
	resul=phi.direct.it(pp,sigmap,X,phi1.old,phi2.old,u1.new,u2.new)
	phi1.new=resul[[1]]
	phi2.new=resul[[2]]
	phi.new=matrix(cbind(phi1.new,phi2.new),2*(k-1),1)
	#phi_new=phi.old+Fk%*%sk
	phi.new=data.matrix(phi.new)
	p3=(phi.new-phi.old)/phi.old
	phi1.new=phi.new[1:(k-1),]
	phi2.new=phi.new[k:(2*(k-1)),]
	mmm=max(abs(c(p1,p3)))
	if (min(phi.new)<0.001){
	aviso=1
	ll=ll+1
	phi1.new=phi1.prim
	phi2.new=phi2.prim
	#beta.new=beta.prim
	u1.new=u1.prim
	u2.new=u2.prim
	if (ll>1)
	{
	phi1.new=phi1.prim
	phi2.new=phi2.prim
	iter1=maxiter1}
	iter1=iter1+1

	}
	if (aviso==0)
	{mmm=max(abs(c(p1,p3)))
	mm=1}
	else
	{mmm=max(abs(p3))
	mm=eps/10}
	if (mmm>eps){iter1=iter1+1} else {iter1=maxiter1}}

	if (cont==1) {
	iter1=maxiter1}

	}
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
resulbeta
sk_F=sPhikf.it(d,t,pp,sigmap,X,theta,phi1.new,phi2.new)
Fk=sk_F[[2]]
phi.new=(rbind(as.matrix(phi1.new),as.matrix(phi2.new)))
cii=ci(phi.new,Fk)
resulphi=cbind(phi.est=phi.new,Std.dev=cii[[1]],p.value=cii[[2]])
colnames(resulphi)=c("Estimate","Std.Error","p.value")
result=list()
result=list(Estimated.probabilities=pr,u1=u1.new,u2=u2.new,mean=mu,warning1=cont,Fisher.information.matrix.beta=F,Fisher.information.matrix.phi=Fk,beta.Stddev.p.value=resulbeta,phi.Stddev.p.value=resulphi,warning2=aviso)
class(result)="mme"
return(result)
}


####################################################################
#' Analytic MSE for Model 2
#'
#' This function calculates the analytic MSE for the multinomial mixed model with two independent random effects
#' for each category of the response variable: one  random effect associated with the domain and another independent random effect associated with time and domain (Model 2).
#' See details of the model and the expresion of mse in Lopez-Vizcaino et al. (2013). The formulas
#' of Prasad and Rao (1990) are adapted to Model 2. This function uses the output of \code{\link[mme]{modelfit2}}.
#'
#' @param p vector with the number of the auxiliary variables per category.
#' @param X list of matrices with the auxiliary variables obtained from \code{\link[mme]{data.mme}}. The dimension of the list is the number of categories of the response variable minus one.
#' @param result the output of the function \code{\link[mme]{modelfit2}}.
#' @param M vector with the area sample sizes.
#' @param MM vector with the population sample sizes.
#' @keywords models
#' @export
#' @seealso \code{\link[mme]{data.mme}}, \code{\link[mme]{initial.values}},
#'\code{\link[mme]{wmatrix}}, \code{\link[mme]{phi.mult.it}},
#' \code{\link[mme]{prmu.time}}, \code{\link[mme]{phi.direct.it}},
#' \code{\link[mme]{sPhikf.it}}, \code{\link[mme]{modelfit2}},
#' \code{\link[mme]{Fbetaf.it}}, \code{\link[mme]{ci}},
#' \code{\link[mme]{mseb}}
#' @return mse.analitic is a matrix with the MSE estimator calculated by adapting the explicit
#' formulas of Prasad and  Rao (1990). The matrix dimension is the number of domains multiplied by the number of categories minus 1.
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Small area estimation of labour force indicator under a multinomial mixed model
#' with correlated time and area effects. Submitted for review.
#' @references Prasad, NGN, Rao, JNK (1990).The estimation of the mean squared error of small
#' area estimators. Journal of the American Statistical Association, 85, 163-171.
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' mod=2 #type of model
#' data(simdata2)
#' datar=data.mme(simdata2,k,pp,mod)
#' ##Model fit
#' result=modelfit2(datar$d,datar$t,pp,datar$Xk,datar$X,datar$Z,datar$initial,datar$y[,1:(k-1)],
#'        datar$n,datar$N)
#'
#' ##Analytic MSE
#' msef=msef.it(pp,datar$X,result,datar$n,datar$N)

msef.it<-function(p,X,result,M,MM){
pr=result$Estimated.probabilities
k=ncol(pr)
phi1.old=(result[[9]][1:(k-1),1])
phi2.old=(result[[9]][k:(2*k-2),1])
F=result$Fisher.information.matrix.phi
D=nrow(MM)
d=nrow(result$u1)
t=D/d
r=sum(p+1)
Vd=list()
qq=matrix(0,(r),(r))
Ld=list()
sigmap=wmatrix(MM,pr)
sigmapMM=wmatrix(M,pr)
sigmapt=list()
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
	Xd[[i]]=do.call(rbind,Xaa[[i]])
}

j=1
Iq1=matrix(rep(diag(rep(1,(k-1))),t),(t*(k-1)),(k-1),byrow=TRUE)
for (i in 1:d){
	phi22=diag(1/phi2.old)
	a=matrix(0,(k-1),(k-1))
	for (l in 1:t){
		sigmapt[[l]]=sigmapMM[[j]]-sigmapMM[[j]]%*%solve(sigmapMM[[j]]+phi22)%*%sigmapMM[[j]]
		a=a+sigmapt[[l]]
		j=j+1}
	Ld[[i]]=do.call("blockdiag",sigmapt)
	Vd[[i]]=Ld[[i]]-Ld[[i]]%*%Iq1%*%solve(diag(1/phi1.old)+a)%*%t(Iq1)%*%Ld[[i]]
	qq=qq+crossprod((Xd[[i]]),(Vd[[i]]))%*%Xd[[i]]
	}

qq=solve(qq)

Xt=do.call(rbind,Xd)
W=bdiag(sigmapMM)
H=bdiag(sigmap)
phi22=list()
for (i in 1:t){
	phi22[[i]]=diag(phi2.old)}
phi22=bdiag(phi22)

z2=diag(rep(1,(t*(k-1))))
T11=list()
T22=list()
T12=list()
Z1i=list()
TT=list()
for (i in 1:d){
	T11[[i]]=as(diag(phi1.old)-diag(phi1.old)%*%t(Iq1)%*%Vd[[i]]%*%Iq1%*%diag(phi1.old),"sparseMatrix")
	T12[[i]]=as(-diag(phi1.old)%*%t(Iq1)%*%Vd[[i]]%*%z2%*%phi22,"sparseMatrix")
	T22[[i]]=as(phi22-phi22%*%t(z2)%*%Vd[[i]]%*%z2%*%phi22,"sparseMatrix")
	TT[[i]]=as(Iq1%*%T11[[i]]%*%t(Iq1)+Iq1%*%T12[[i]]%*%t(as(z2,"matrix"))+z2%*%t(as(T12[[i]],"matrix"))%*%t(Iq1)+z2%*%t(as(T22[[i]],"matrix"))%*%t(as(z2,"matrix")),"sparseMatrix")
	Z1i[[i]]=as(Iq1,"sparseMatrix")}

pro=list()
prod=list()
pro2=list()
pro22=list()
prod2=list()

jj=1
for (i in 1:d){
	for (j in 1:t){
		pro[[jj]]=as(t(cbind(matrix(0,(k-1),(j-1)*(k-1)),sigmap[[jj]],matrix(0,(k-1),(t-j)*(k-1)))),"sparseMatrix")
		pro22[[j]]=as(sigmapMM[[jj]]%*%X[[jj]],"sparseMatrix")
		jj=jj+1}
		pro2[[i]]=do.call(rBind,pro22)}

G1=list()
G11=list()
G12=list()
G22=list()
A21=list()
A22=list()
G2=list()
jj=1
for (i in 1:d){
	for (l in 1:t){
		G11[[jj]]=as(sigmap[[jj]]%*%(diag(phi1.old)-diag(phi1.old)%*%t(Iq1)%*%Vd[[i]]%*%Iq1%*%diag(phi1.old))%*%t(as(sigmap[[jj]],"matrix")),"sparseMatrix")
		G12[[jj]]=as(-sigmap[[jj]]%*%diag(phi1.old)%*%t(Iq1)%*%Vd[[i]]%*%z2%*%phi22%*%(pro[[jj]]),"sparseMatrix")
		G22[[jj]]=as(t(as(pro[[jj]],"matrix"))%*%(phi22-phi22%*%t(as(z2,"matrix"))%*%Vd[[i]]%*%z2%*%phi22)%*%pro[[jj]],"sparseMatrix")
		A21[[jj]]=as(sigmap[[jj]]%*%X[[jj]],"sparseMatrix")
		A22[[jj]]=as(t(as(pro[[jj]],"matrix"))%*%TT[[i]]%*%(pro2[[i]]),"sparseMatrix")

		G2[[jj]]=as((A21[[jj]]-A22[[jj]])%*%qq%*%(t(as(A21[[jj]],"matrix"))-t(as(A22[[jj]],"matrix"))),"sparseMatrix")
		G1[[jj]]=G11[[jj]]+G12[[jj]]+t(as(G12[[jj]],"matrix"))+G22[[jj]]

		jj=jj+1
		}}


Vu1=diag(rep(phi1.old,d))
Vu2=diag(rep(phi2.old,(d*t)))
V=bdiag(Vd)
Z1=bdiag(Z1i)
R2=Vu2%*%V
R1=Z1%*%Vu1%*%t(as(Z1,"matrix"))%*%V

#Esto hay que rehacerlo#
aa=nrow(R1)
bb=nrow(R2)
delta11=list()
Vakd=list()
Va=list()
Lk=list()
o=1
for(i in 1:2){
	for (j in 1:(k-1)){

			delta=matrix(0,(k-1),(k-1))
			for (ll in 1:(k-1)){
				if (ll==j) {delta[ll,ll]=1}}
			for (ll in 1:t){
				delta11[[ll]]=delta
				}
			#Vakd hai que facela completa
			if (i==1){Vakd =as(Iq1%*%delta%*%t(Iq1),"sparseMatrix")}
			if (i==2){Vakd=bdiag(delta11)}

			for(l in 1:d){
				Va[[l]]=Vakd}
			if (i==1) {Lk[[o]]=(as(diag(rep(1,aa)),"sparseMatrix")-R1)%*%bdiag(Va)%*%V}
			if (i==2) {Lk[[o]]=(as(diag(rep(1,bb)),"sparseMatrix")-R2)%*%bdiag(Va)%*%V}
			o=o+1
				}}

A=list()
for(i in 1:D){
	n0=(k-1)*(D-i)
	if (n0<0) {n0=0}
	A[[i]]=cbind(matrix(0,nrow=k-1,ncol=(k-1)*(i-1)),diag(k-1),matrix(0,nrow=k-1,ncol=n0))}
aa=2*(k-1)

jj=1

gg3=list()
g3=as(matrix(0,(d*t*(k-1)),(d*t*(k-1))),"sparseMatrix")
for (i in 1:aa){
	for (j in 1:aa){
		g33=H%*%Lk[[i]]%*%V%*%t(as(Lk[[j]],"matrix"))%*%t(as(H,"matrix"))
		g3=g3+F[i,j]*g33
		jj=jj+1
		}}
for (kk in 1:D){
	gg3[[kk]]=A[[kk]]%*%g3%*%t(as(A[[kk]],"matrix"))}

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
return(mse)
}
