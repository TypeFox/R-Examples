refit.mhtp_nodiag=function(object,Ynew,z,grp,fix,rand,alpha,step,num,ordre,m,show,IT,maxq,speed,...)
{

	#même Z et grp que object, seul le Y est différent
data=object$data$X
if(missing(z)) z=object$data$z
if(missing(grp)) grp=object$data$grp
ntot=nrow(data)	
p=ncol(data)
q=ncol(z)


	# #on reprend les mêmes argument que dans object$arg
if(missing(fix)) fix=object$arg$fix
if(missing(rand)) {rand_sauv=object$arg$rand_sauv}else{rand_sauv=rand}
if(missing(num)) num=object$arg$num
if(missing(ordre)) ordre=object$arg$ordre
	if(!ordre%in%c("bolasso","pval","pval_hd","FR")) stop('ordre has to be one of "bolasso","pval","pval_hd","FR" ')
if(missing(m)) m=object$arg$m
if(missing(IT)) IT=object$arg$IT
if(missing(maxq)) {maxqdep=object$arg$maxqdep}else{maxqdep=maxq}
if(missing(show)) show=object$arg$show
		showordre=show[1]
		showresult=show[2]
		showit=show[3]	
if(missing(step)) step=object$arg$step
if(missing(speed)){speed=1}
choix_ordre=ordre

	if(missing(Ynew)){stop('Ynew is missing')}
	if(length(Ynew)!=ntot){stop(paste(" 'data' and 'Ynew' must have the same length:",ntot))}
	
alpha=as.numeric(colnames(object$fitted.values))  	     

dec=decompbaseortho(data)
nonind=dec$nonind
U=dec$U
if(p>(ntot+length(nonind)))
{nonind2=c(nonind,(ntot+1+length(nonind)):p)	
	Uchap=U[,-nonind2]}else{if(length(nonind)==0){Uchap=U}else{Uchap=U[,-nonind]}}
dim_X=ncol(Uchap)			#nombre de variables utiles

if(nrow(object$ordrebeta)==1)#on a juste fait un procbol, donc les quantiles sont dans un tableaux de dim alpha*maxq*NBR
{
aV2=array(0,c(length(alpha),dim(object$quantile)[2],dim_X,1)) #on y met tous les quantiles
aV2[,,1:dim(object$quantile)[3],1]=object$quantile
}else{aV2=object$quantile}

	Y=Ynew
	ORDREBETA2=object$ordrebeta
	compteurordre=nrow(ORDREBETA2)
	indice2=NULL
	for(i in 1:nrow(object$ordrebeta))
	{j=1
	while(sum(aV2[,,j,i]^2)==0) 
	{j=j+1}
	k=j
	while(sum(aV2[,,k,i]^2)!=0) 
	{k=k+1}
	k=k-1
	indice2=cbind(indice2,c(j,k))
	}
	
#construction des Z_i à partir de z et grp
grp=rbind(grp)
Ni=max(grp[1,])
for(k in 1:q)
{
Z=matrix(0,nrow=ntot,ncol=Ni)

for(i in 1:Ni)
	{
		indic=which(grp[1,]==names(table(grp[1,])[i]))
		Z[indic,i]=z[indic,k]
		}
assign(paste("Z",k,sep="_"),Z)
}


#----------------------
# dans N_k:le nombre de groupe dans la structure k
# dans Z_k:la matrice de structure
#----------------------

BETA_HAT=matrix(0,nrow=p,ncol=length(alpha))
kchap=numeric(0)
SIGMAE=numeric(0)
indice=numeric(0)
compteur=0
rand=rand_sauv
		
##############################################################################################################
	#################################### Initialisation beta #####################################################
	##############################################################################################################	
	a=table(rand)
	b=names(a)[names(a)!=0]
	set_random=b#[b!=1]#on enleve l'intercept, parce qu'il sera tjs dans les fixes, on ne le compte pas dans fixes+aleatoire
	var_nonselect=fix+length(set_random)
	
	correspondance=1:p
	if(length(set_random)>0)
	{a=c(1:fix,as.numeric(set_random),(1:p)[-c(1:fix,as.numeric(set_random))])
	}else{a=c(1:fix,(1:p)[-c(1:fix)])}
		correspondance=rbind(correspondance,a)
	data2=data[,correspondance[2,]]		
		
a=matrix(runif(ntot*max(m,2),0,ntot),nrow=ntot)
random=ceiling(a)
		
UCHAP=NULL
y.fitted=NULL
COMPTEUR=NULL
PSI=array(0,c(q,q,length(alpha)))
for(alph in 1:length(alpha)) #boucle sur alpha
{
	
	#initialisation
	if(showit){print("initialization")}
mm2=MM2(data2=data2,yhat=Y,m=m,sigma=q,maxqdep=maxqdep,num=num,var_nonselect=var_nonselect,random=random,showit=showit,showresult=showresult,correspondance=correspondance,ORDREBETA2=ORDREBETA2,indice=indice2,aV2=aV2,choix_ordre=choix_ordre,alpha=alpha,alph=alph,IT=IT)
	
aV2=mm2$aV2
ORDREBETA2=mm2$ORDREBETA2
ind=mm2$ind
indice2=mm2$indice
compteurordre=mm2$compteurordre	
		if(showresult){print(paste("number of selected variables:",length(ind)))}

	
	#on cherche le nombre de variables fixes+aléatoires
rand=rand_sauv
a=table(rand)
b=names(a)[names(a)!=0]
set_random=b#[b!=1]#on enleve l'intercept, parce qu'il sera tjs dans les fixes, on ne le compte pas dans fixes+aleatoire

var_nonselect=fix+length(set_random)

		#on initialise tout
Psi=matrix(0.02,nrow=q,ncol=q)
diag(Psi)=1:q
Psi1=solve(Psi)
	##############################################################################################################
	#################################### Initialisation beta #####################################################
	##############################################################################################################	
	
	#ind=ind_sauv[alph,]
	beta_hat=matrix(0,p)
	if(length(ind)>0){
		reg=lm(Y~data[,ind]-1)
		
		beta_hat[ind]=reg$coefficients
		beta_hat[-which(beta_hat!=0)]=0}
		sigma_e=as.vector(var(Y-data%*%beta_hat))

	sum1=10
	sum2=rep(10,q)
	condition=1
	ft=numeric(0)
	uchap=NULL
	
	compteur_ijk=0
	converge=TRUE
	delete=0
Z=numeric(0)
for(k in 1:q)
{Z=cbind(Z,get(paste("Z",k,sep="_")))}
N=ncol(Z)#}
	while(((sum1>10^-8)||(sum(sum2)>10^-8))||(abs(condition)>10^-6))
	{
		compteur_ijk=compteur_ijk+1
		if(showit){print(paste("iteration=",compteur_ijk))}
		compteur=compteur+1
		#on cherche beta en faisant un lasso a sigma_u et u fixé
		
##############################################################################################################
#################################### E-step ##########################################################
##############################################################################################################
G=kronecker(Psi,diag(1,Ni))
G1=kronecker(Psi1,diag(1,Ni))

yhat=Y-data%*%beta_hat
a=t(Z)%*%Z+sigma_e*G1
S1=solve(a)
uchap3=S1%*%t(Z)%*%yhat
if(length(uchap3)!=length(uchap)){uchap=matrix(0,N)	}
sum2=sum((uchap3-uchap)^2)
uchap=uchap3
yhat=as.vector(Y-Z%*%uchap)
##############################################################################################################
#################################### M-step -1 ##########################################################
##############################################################################################################

if(((((sum1>10^-6)|(sum(sum2)>10^-6)))&(speed==1))|(speed==0)|(delete==1))
{	
correspondance=1:p
if(length(set_random)>0)
{a=c(1:fix,as.numeric(set_random),(1:p)[-c(1:fix,as.numeric(set_random))])
	}else{a=c(1:fix,(1:p)[-c(1:fix)])}
correspondance=rbind(correspondance,a)
data2=data[,correspondance[2,]]		#les var_nonselect premières colonnes sont les fixed et les random
#mean that data2[ordre]=data[correspondance[2,ordre]]=data[,bb]

	#on ne calcule que les manquants
mm2=MM2(data2=data2,yhat=yhat,m=m,sigma=sigma_e,maxqdep=maxqdep,num=num,var_nonselect=var_nonselect,random=random,showit=showit,showresult=showresult,correspondance=correspondance,ORDREBETA2=ORDREBETA2,indice=indice2,aV2=aV2,choix_ordre=choix_ordre,alpha=alpha,alph=alph,IT=IT)

aV2=mm2$aV2
ORDREBETA2=mm2$ORDREBETA2
ind=mm2$ind
indice2=mm2$indice
compteurordre=mm2$compteurordre	
}
	

if(showresult){print(paste("number of selected variables:",length(ind)))}
beta_hat2=matrix(0,p)

if(length(ind)>0){
	reg=lm(yhat~data[,ind]-1)

	beta_hat2[ind]=reg$coefficients
	beta_hat2[-which(beta_hat2!=0)]=0
	sum1=sum((beta_hat2-beta_hat)^2)
	#print(sum1)
	}else{sum1=10}
beta_hat=beta_hat2


##############################################################################################################
#################################### E-step -2 ##########################################################
##############################################################################################################


	yhat2=Y-data%*%beta_hat
	uchap3=S1%*%t(Z)%*%yhat2	
	if(length(uchap3)!=length(uchap)){uchap=matrix(0,N)	}
	sum2=sum((uchap3-uchap)^2)
	uchap=uchap3

	a=1
	for(k in 1:q)
	{uchapal=uchap[a:(a-1+Ni)]
		a=a+Ni
		assign(paste("uchap",k,sep="_"),uchapal)
	}

##############################################################################################################
#################################### M-step -2 ##########################################################
##############################################################################################################
a=t(Z)%*%Z+sigma_e*G1
Tkk=solve(a)

a=0
for(i in 1:q)
{
b=a
	for(j in i:q)
	{	
	assign(paste("T",i,j,sep=""),Tkk[(a+1):(a+Ni),(b+1):(b+Ni)])
	b=b+Ni
	}
	a=a+Ni}
U=NULL
for(k in 1:q)
U=cbind(U,get(paste("uchap",k,sep="_")))
E=t(U)%*%U #esperance condi de chaque u_i*u_j
#print(E)

#on estime G maintenant, Iq paramètre (q(q+1)/2)	
for(i in 1:(q))
{for(j in (i):q)
	{
	Psi[i,j]=(E[i,j]+sigma_e*sum(diag(get(paste("T",i,j,sep="")))))/Ni
}}
	
	
for(i in 1:(q))
{for(j in (i):q)
	{
	Psi[j,i]=Psi[i,j]#get(paste("s",i,j,sep=""))[1]
}}
Psi1=solve(Psi)
	
		
#on calcule sigma_e sur les résidus du modele
sumk=Z%*%uchap
yhatt=data%*%beta_hat+sumk

aa=Tkk%*%G1 #(T*G-1)
a=N-sigma_e*sum(diag(aa))
sigma_e=sum((Y-yhatt)^2)/ntot+(a*sigma_e/ntot)


if(compteur_ijk>step){message("Algorithm did not converge..")
		converge=FALSE
	break}

if(showit)
{print("Psi")
	print(Psi)}
	
suma=0
sumk=0
sumj=0
for(k in 1:q)
{sumk=sumk+get(paste("Z",k,sep="_"))%*%get(paste("uchap",k,sep="_"))}
suma=log(det(G))+t(uchap)%*%G1%*%uchap

yhatt=data%*%beta_hat+sumk

ft=c(ft,ntot*log(sigma_e)+suma+sum((Y-yhatt)^2)/sigma_e)

if(length(ft)>2){condition=ft[length(ft)]-ft[length(ft)-1]}else{condition=10}



}#fin while
	
uchap=NULL
for(k in 1:q)
{uchap=c(uchap,get(paste("uchap",k,sep="_")))}
		

PSI[,,alph]=Psi
SIGMAE=c(SIGMAE,sigma_e)
BETA_HAT[,alph]=beta_hat
kchap=c(kchap,length(ind))
UCHAP=cbind(UCHAP,uchap)
y.fitted=cbind(y.fitted,yhatt)
COMPTEUR=c(COMPTEUR,compteur_ijk)
}#fin alph

rownames(aV2)=paste("alpha=",alpha)
colnames(aV2)=paste("Hk,",0:(maxqdep-1))
dimnames(aV2)[[3]]=paste("ktest=",(var_nonselect!=0):(dim_X-(var_nonselect==0)))
dimnames(aV2)[[4]]=paste("ordre",1:compteurordre)

colnames(y.fitted)=alpha
#on doit faire une liste avec tous les quantiles utiles
rownames(ORDREBETA2)=paste("ordre",1:compteurordre)
colnames(BETA_HAT)=alpha

out=list(data=list(X=data,Y=Y,z=z,grp=grp),beta=BETA_HAT,fitted.values=y.fitted,u=UCHAP,Psi=PSI,sigma_e=SIGMAE,it=COMPTEUR,quantile=aV2,ordrebeta=ORDREBETA2,converge=converge,call=match.call(),arg=list(fix=fix,rand_sauv=rand_sauv,step=step,num=num,ordre=choix_ordre,m=m,random=random,show=show,IT=IT,maxqdep=maxqdep))

out
structure(out,class=c("mhtp","mhtp_nodiag"))
		
}#fin refit procbol
