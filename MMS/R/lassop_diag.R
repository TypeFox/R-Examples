lassop_diag=function(data,Y,z,grp,mu,step,fix,rand,penalty.factor,alpha,showit)
{
	
	#-----------------------------------	
#	data=matrice data, the first column should be 1, for the intercept
#	Y=observation	Y=data*beta
#	z=random effects matrix, n*q, c'est la matrice à partager en trnasformer en matrice de structure Z
#	grp= groups, q*n, o autorise plusieurs structures, donc q lignes
#	mu=regularization parameter
#	fix= le nombre de variables qu'on ne veut pas séléctionner, ce sont les premières colonnes dans data
#	rand= si z contient des variables à la fois fixes et aléatoires, il est conseillé de ne pas les soumettre à séléction. un vecteur de longueur q, 0 si la variable n'est pas fixe+aléatoire, sa position dans data sinon.
#	penalty.factor=Separate penalty factors can be applied to each coefficient. This is a number that multiplies lambda to allow differential shrinkage. Can be 0 for some variables, which implies no shrinkage, and that variable is always included in the model. Default is 1 for all variables that are not in rand and in 1:fix.
# 	The elasticnet mixing parameter, with 0<α≤ 1. The penalty is defined as
#		(1-α)/2||β||_2^2+α||β||_1.
#alpha=1 is the lasso penalty.
#	showit=affiche les itérations de l'algorithme
#-----------------------------------------
	
ntot=nrow(data)	
p=ncol(data)
q=ncol(z)

	#safety
	if(length(Y)!=ntot){stop(" 'data' and 'Y' must have the same length ")}
	if(nrow(grp)!=q){stop("grp has to have q row(s)")}
	if(missing(z)||missing(grp)){stop("error")}
	if(missing(showit)){showit=FALSE}	
	if(missing(mu)){stop(" please, don't forget the regularization parameter..")}
	if(length(mu)>1){stop('mu has to be a number')}
	if(missing(rand)){rand=rep(0,q)}
	if(missing(alpha)){alpha=1}
	if(missing(step)){step=3000}

#verifier que rand est en accord entre z et data
if(length(rand)!=ncol(z)){stop("rand has to be of length nrow(grp)")}

I=which(rand!=0)
for(k in I)
{if(sum((z[,k]-data[,rand[k]])^2)>10^(-10)){stop(paste('Are you sure that random effect number',k,' is in position',rand[k],'of the data matrix?'))}
	}
	

#		-------------------------------------
#			on scale la matrice de départ
#		-------------------------------------

#si la matrice de depart ne contient pas l'intercept (colonne de 1) on la rajoute et on rajoute 1 dans var_nonselect s'il n'etait pas manquant
if(sum(data[,1])==ntot){data=cbind(data[,1],scale(data[,-1])*sqrt(ntot)/sqrt(ntot-1))
		data[which(is.na(data))]=0
	intercept=TRUE
			}else{
				data=scale(data)*sqrt(ntot)/sqrt(ntot-1)
				data[which(is.na(data))]=0
				data=cbind(1,data)
				intercept=FALSE
				p=p+1
			}
			
if(missing(fix)){fix=1
}else{if(!intercept){fix=fix+1}}			

if(fix<=0){stop("var_nonselect has to be positive, the intercept is not submitted to selection, not implemented in glmnet..")}

if(!intercept){rand[rand!=0]=rand[rand!=0]+1}

#on met a 0 les variables qui resteront non selectionnées car deja dans fix
for(i in rand)
{if(i<=fix)
	{rand[which(rand==i)]=0}}

p=ncol(data)

#construction des Z_i à partir de z et grp
grp=rbind(grp)
for(k in 1:q)
{
assign(paste("grp",k,sep="_"),factor(grp[k,]))
assign(paste("N",k,sep="_"),max(grp[k,]))

Z=matrix(0,nrow=ntot,ncol=get(paste("N",k,sep="_")))

for(i in 1:get(paste("N",k,sep="_")))
	{
		indic=which(get(paste("grp",k,sep="_"))==names(table(get(paste("grp",k,sep="_")))[i]))
		Z[indic,i]=z[indic,k]
		}
assign(paste("Z",k,sep="_"),Z)
}



#----------------------
# dans N_k:le nombre de groupe dans la structure k
# dans Z_k:la matrice de structure
#----------------------

#on initialise tout
sigma_uinit=rep(1,q)
sigma_u=sigma_uinit

#on cherche le nombre de variables fixes+aléatoires
a=table(rand)
b=names(a)[names(a)!=0]
set_random=b#[b!=1]#on enleve l'intercept, parce qu'il sera tjs dans les fixes, on ne le compte pas dans fixes+aleatoire

############################
# alg51bolm
###########################

	Iq=1:q #set des effets aléatoire non nul
	nonIq=numeric(0)	
	sigma_e=1
	
##############################################################################################################
	#################################### Initialisation beta #####################################################
	##############################################################################################################	
	if(missing(penalty.factor))
{vec=c(rep(0,fix),rep(1,p-fix))
}else{
	if(intercept){
		if(length(penalty.factor)!=p){stop('penalty.factor has to be of length p')}
		vec=penalty.factor}else{
			if(length(penalty.factor)!=(p-1)){stop('penalty.factor has to be of length p')}
			vec=c(0,penalty.factor)}
	}
	vec[as.numeric(set_random)]=0

	#on cherche beta en faisant un lasso a sigma_u et u fixé
	a=glmnet(data,Y,family="gaussian",lambda=mu*sigma_e,penalty.factor=vec,alpha=alpha)
	ind=which(as.vector(a$beta)!=0)
	ind=c(1,ind)
	
	beta_hat=as.vector(a$beta)
	beta_hat[1]=a$a0
	
	sum1=10
	sum2=rep(10,q)
	condition=1
	ft=numeric(0)
	uchap=NULL
	
	compteur_ijk=0
	fit=0
	converge=TRUE
	delete=1

	while(((sum1>10^-8)||(sum(sum2)>10^-8))||(abs(condition)>10^-6))
	{
		compteur_ijk=compteur_ijk+1
		if(showit==TRUE){print(paste("iteration=",compteur_ijk))}
	##############################################################################################################
#################################### E-step ##########################################################
##############################################################################################################
if(length(Iq)>0)
{
d=numeric(0)
for(k in Iq)
{d=c(d,rep(as.vector(sigma_u)[k],get(paste("N",k,sep="_"))))}
G=diag(d)
	
if(delete==1){
Z=numeric(0)
for(k in Iq)
{Z=cbind(Z,get(paste("Z",k,sep="_")))}
N=ncol(Z)}
delete=0
G1=solve(G)

yhat=Y-data%*%beta_hat
a=t(Z)%*%Z/sigma_e+G1
S1=solve(a)
uchap3=S1%*%t(Z)%*%yhat/sigma_e
if(length(uchap3)!=length(uchap)){uchap=matrix(0,N)	}
sum2=sum((uchap3-uchap)^2)
uchap=uchap3
yhat=as.vector(Y-Z%*%uchap)
}else{yhat=Y}

##############################################################################################################
#################################### M-step -1 ##########################################################
##############################################################################################################

#on construit un vecteur de length p qui est le penalty.factor du lasso
	if(missing(penalty.factor))
{	vec=c(rep(0,fix),rep(1,p-fix))
}else{
	if(intercept){vec=penalty.factor}else{vec=c(0,penalty.factor)}
	}
	vec[as.numeric(set_random)]=0

	
#on cherche beta en faisant un lasso a sigma_u et u fixé
a=glmnet(data,yhat,alpha=alpha,family="gaussian",lambda=mu*sigma_e,penalty.factor=vec)
ind=which(as.vector(a$beta)!=0)
ind=c(1,ind)

if(length(ind)>=min(ntot,p)){print("variables selected>min(ntot,p), please increase 'mu'")
	beta=rep(NA,p)
	sigma_u=rep(NA,q)
	sigma_e=NA
	yhatt=rep(NA,ntot)
	uchap=rep(NA,N)
	break}
if(showit){print(paste("number of selected variables:",length(ind)))}

beta_hat2=as.vector(a$beta)
beta_hat2[1]=a$a0
	 sum1=sum((beta_hat2-beta_hat)^2)
	 beta_hat=beta_hat2
##############################################################################################################
#################################### E-step -2 ##########################################################
##############################################################################################################

if(length(Iq)>0)
{
#on calcule les nouveaux uchapal[t+1] a partir de uchap[t]
	yhat2=Y-data%*%beta_hat#-sumj
	uchap3=S1%*%t(Z)%*%yhat2/sigma_e	
	if(length(uchap3)!=length(uchap)){uchap=matrix(0,N)	}
	uchap=uchap3

	a=1
	for(k in Iq)
	{uchapal=uchap[a:(a-1+get(paste("N",k,sep="_")))]
		a=a+get(paste("N",k,sep="_"))
		if(compteur_ijk>1)
		{sum2[k]=sum((uchapal-get(paste("uchap",k,sep="_")))^2)}else{sum2[k]=sum(uchapal^2)}
		assign(paste("uchap",k,sep="_"),uchapal)
	}

##############################################################################################################
#################################### M-step -2 ##########################################################
##############################################################################################################
for(k in Iq)
{	
	assign(paste("c",k,k,sep=""),t(get(paste("Z",k,sep="_")))%*%get(paste("Z",k,sep="_"))+diag(sigma_e/sigma_u[k],get(paste("N",k,sep="_"))))
	
	if(compteur_ijk==1)
	{for(j in Iq[-which(Iq<=k)])
		{assign(paste("c",k,j,sep=""),t(get(paste("Z",k,sep="_")))%*%get(paste("Z",j,sep="_")))}
	for(j in Iq[-which(Iq>=k)])
		{assign(paste("c",k,j,sep=""),t(get(paste("c",j,k,sep=""))))}
	}
}

A=numeric(0)
for(k in Iq)
{B=numeric(0)
	for(j in Iq)
	{B=cbind(B,get(paste("c",k,j,sep="")))}
	A=rbind(A,B) 
	}	
g1=solve(A)

a=0
for(k in Iq)
{
	assign(paste("C",k,k,sep=""),sum(diag(g1[((a+1):(a+get(paste("N",k,sep="_")))),((a+1):(a+get(paste("N",k,sep="_"))))])))
	a=a+get(paste("N",k,sep="_"))
	}
	
sigma_utemp=sigma_u
for(k in Iq)
{sigma_utemp[k]=sum((get(paste("uchap",k,sep="_")))^2/get(paste("N",k,sep="_")))+get(paste("C",k,k,sep=""))*sigma_e/get(paste("N",k,sep="_"))
}


#on calcule sigma_e sur les résidus du modele
sumk=Z%*%uchap
yhatt=data%*%beta_hat+sumk

a=0
for(k in Iq)
{a=a+get(paste("N",k,sep="_"))-sigma_e/sigma_utemp[k]*get(paste("C",k,k,sep=""))}
sigma_e=sum((Y-yhatt)^2)/ntot+(a*sigma_e/ntot)

sigma_u=sigma_utemp

if(showit){print(paste("sigma_u=",sigma_utemp))
	print(paste("sigma_e=",sigma_e))}
	
#on supprime les effets alétoire nul de Iq
for(k in Iq)
{
	if(var(get(paste("uchap",k,sep="_")))<10^-4*(min(sigma_e,1))){sigma_u[k]=0
		Iq=Iq[-which(Iq==k)]
		rand[k]=0
		delete=1

		a=table(rand)
		b=names(a)[names(a)!=0]
		set_random=b#[b!=1]#on enleve l'intercept, parce qu'il sera tjs dans les fixes, on ne le compte pas dans fixes+aleatoire

		var_nonselect=fix+length(set_random)
		nonIq=c(nonIq,k)
		sum2[k]=0
if(abs(mean(get(paste("uchap",k,sep="_"))))<10^-5){assign(paste("uchap",k,sep="_"),rep(0,get(paste("N",k,sep="_"))))}
}}


#fin length Iq 0
}else{yhatt=data%*%beta_hat
	sum1=0
	sum2=0
	condition=0}#fin length IQ>0


if((length(ind)==0)&&(length(Iq)==0)){break}
if(compteur_ijk>step){message("Algorithm can't converge..")
	converge=FALSE
	break}


suma=0
sumk=0
sumj=0
for(k in Iq)
{suma=suma+sum(get(paste("uchap",k,sep="_"))^2)/sigma_u[k]
	sumk=sumk+get(paste("Z",k,sep="_"))%*%get(paste("uchap",k,sep="_"))
	sumj=sumj+get(paste("N",k,sep="_"))*log(sigma_u[k])}

yhatt=data%*%beta_hat+sumk

ft=c(ft,ntot*log(sigma_e)+sumj+suma+sum((Y-yhatt)^2)/sigma_e)
if(length(ft)>2){condition=ft[length(ft)]-ft[length(ft)-1]}else{condition=10}


}#fin while

#########################
# fin alg51bolm
#########################
d=numeric(0)
for(k in 1:q)
{d=c(d,rep(as.vector(sigma_u)[k],get(paste("N",k,sep="_"))))}
G=diag(d)

d=numeric(0)
for(k in 1:q)
{d=c(d,as.vector(sigma_u)[k])}
Psi=diag(d,q)

names(beta_hat)=c("Intercept",paste("X", 2:p,sep=""))


out=list(data=list(X=data,Y=Y,z=z,grp=grp),beta=beta_hat,fitted.values=yhatt,Psi=Psi,sigma_e=sigma_e,it=compteur_ijk,converge=converge,u=uchap,call=match.call(),mu=mu)

out
structure(out,class="lassop")

}#fin function
