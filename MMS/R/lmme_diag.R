lmme_diag=function(data,Y,z,grp,step,showit)
{
	
	#-----------------------------------	
#	data=matrice data, the first column should be 1, for the intercept
#	Y=observation	Y=data*beta
#	z=random effects matrix, n*q, c'est la matrice à partager en trnasformer en matrice de structure Z
#	grp= groups, q*n, o autorise plusieurs structures, donc q lignes
#	diag= Psi est diagonale?
#	showit=affiche les itérations de l'algorithme
#	step=nombre max d'itérations de l'algorithme
#-----------------------------------------
	
ntot=nrow(data)	
p=ncol(data)
q=ncol(z)

	#safety
	if(ntot<p){stop('p has to be lower than n')}
	if(length(Y)!=ntot){stop(" 'data' and 'Y' must have the same length ")}
	if(nrow(grp)!=q){stop("grp has to have q row(s)")}
	if(missing(z)||missing(grp)){stop("error")}
	if(missing(showit)){showit=FALSE}	
	if(missing(step)){step=3000}

	
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


############################
# alg51bolm
###########################


	Iq=1:q #set des effets aléatoire non nul
	nonIq=numeric(0)
	
	sigma_e=1
	reg=lm(Y~data-1)
	beta_hat2=reg$coefficients
	beta_hat2[-which(beta_hat2!=0)]=0
	beta_hat=beta_hat2
#	beta_hat=matrix(0,p)
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
		if(showit){print(paste("iteration=",compteur_ijk))}
	##############################################################################################################
#################################### E-step ##########################################################
##############################################################################################################
if(length(Iq)>0)
{
d=numeric(0)
for(k in Iq)
{d=c(d,rep(as.vector(sigma_u)[k],get(paste("N",k,sep="_"))))}
G1=diag(1/d)
	
if(delete==1){
Z=numeric(0)
for(k in Iq)
{Z=cbind(Z,get(paste("Z",k,sep="_")))}
N=ncol(Z)}
delete=0

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

	reg=lm(yhat~data-1)
	beta_hat2=reg$coefficients
	beta_hat2[-which(beta_hat2!=0)]=0
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
	sum2=sum((uchap3-uchap)^2)
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
		if(var(get(paste("uchap",k,sep="_")))<(10^-4*sigma_e)){sigma_u[k]=0
#if (sigma_u[k]<10^-10){sigma_u[k]=0  #regarder si ce ne serait pas mieux de mettre a zero si le vecteur est a 0, parce que si =une constante?
		Iq=Iq[-which(Iq==k)]
		delete=1

		nonIq=c(nonIq,k)
		sum2[k]=0
if(abs(mean(get(paste("uchap",k,sep="_"))))<10^-5){assign(paste("uchap",k,sep="_"),rep(0,get(paste("N",k,sep="_"))))}
}}


#fin length Iq 0
}else{yhatt=data%*%beta_hat
	sum1=0
	sum2=0
	condition=0}#fin length IQ>0


if((length(Iq)==0)){break}
if(compteur_ijk>step){message("Algorithm did not converge..")
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

d=numeric(0)
for(k in 1:q)
{d=c(d,rep(as.vector(sigma_u)[k],get(paste("N",k,sep="_"))))}
G=diag(d)

d=numeric(0)
for(k in 1:q)
{d=c(d,as.vector(sigma_u)[k])}
Psi=diag(d,q)

names(beta_hat)=c("Intercept",paste("X", 2:p,sep=""))

#########################
# fin alg51bolm
#########################
out=list(data=list(X=data,Y=Y,z=z,grp=grp),beta=beta_hat,Psi=Psi,sigma_e=sigma_e,fitted.values=yhatt,it=compteur_ijk,converge=converge,u=uchap,call=match.call())
structure(out,class="lmme")


}#fin function
