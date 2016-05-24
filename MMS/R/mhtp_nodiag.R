mhtp_nodiag=function(data,Y,z,grp,fix,rand,alpha,step,num,ordre,m,show,IT,maxq,speed)
{
#-----------------------------------	
#	data=matrice data, the first column should be 1, for the intercept
#	Y=observation	Y=data*beta
#	z=random effects matrix, n*q, c'est la matrice à partager en transformer en matrice de structure Z. Si des effets aléatoires sont aussi des colonnes de data, elles doivent etre en premier dans z et il est conseillé de les dénombrer dans rand
#	grp= groups, q*n, on autorise plusieurs structures, donc q lignes
#	fix= le nombre de variables qu'on ne veut pas sélectionner, ce sont les premières colonnes dans data
#	rand= si z contient des variables à la fois fixes et aléatoires, il est conseillé de ne pas les soumettre à sélection. un vecteur de longueur q, 0 si la variable n'est pas fixe+aléatoire, sa position dans data sinon.

#	alpha=erreur de première espèce du test
#	step=number of max steps of the algorithm
#	num= nombre de variables max qu'on ordonne
#	ordre=avec quelle méthode on ordonne (pval, pval_hd ou bolasso)
#	m=nombre d'iteration bootstrap pour le bolasso
#	IT=nombre de simulations pour le calcul du quantile
#	maxq=nombre max d'hypotheses alternative testees
#	show=c(showordre,showresult,showit). 
	#	showordre=affiche l'ordre au fur et à mesure
	#	showresult=affiche les resultats des tests
	#	showit=affiche les itérations de l'algorithme
#	speed=c(0,1). Should the algorithm be speed up once everything but the lokelihood has "converged"?
#-----------------------------------------
	
	
	#probleme: si random=1:3,donc ordre= (1 2 3) 34 etc et qu'on supprime le 3 des random mais qu'il reste en 3eme position une fois ordonné, on zap le calcul du quantile pour l'ordre (1,2) 3 etc, a vérifier
	
ntot=nrow(data)	
p=ncol(data)
q=ncol(z)

	#safety
	if(length(Y)!=ntot){stop(" 'data' and 'Y' must have the same length ")}
	if(nrow(grp)!=q){stop("grp has to have q row(s)")}
	if(missing(z)||missing(grp)){stop("error, missing grp")}
	if(missing(m)){m=10}
	if(missing(num)){num=min((ntot-1)/2,(p-1)/2,30)}
	if(missing(alpha)){alpha=c(0.1)}
	alpha=sort(alpha,decreasing=TRUE)
	if(missing(ordre)){ordre="bolasso"}
	if(!ordre%in%c("bolasso","pval","pval_hd")) stop('ordre has to be one of "bolasso","pval","pval_hd" ')
	if(missing(IT)){IT=1000}
	if(missing(maxq)){maxq=min(log(min(ntot,p)-1,2),5)}
	if(missing(show)){show=c(0,0,0)}
		showordre=show[1]
		showresult=show[2]
		showit=show[3]	
	if(missing(rand)){rand=rep(0,q)}
	if(missing(step)){step=3000}
	if(missing(speed)){speed=1}
choix_ordre=ordre

#verifier que rand est en accord entre z et data
if(length(rand)!=ncol(z)){stop("rand has to be of length nrow(grp)")}

I=which(rand!=0)
for(k in I)
{if(sum((z[,k]-data[,rand[k]])^2)>10^(-10)){stop(paste('Are you sure that random effect number',k,' is in position',rand[k],'of the data matrix?'))}
	}


#		-------------------------------------
#			on scale la matrice de départ
#		-------------------------------------

#si la matrice de depart ne contient pas l'intercept (colonne de 1) on la rajoute et on rajoute 1 dans fix s'il n'etait pas manquant
if(sum(data[,1])==ntot){data=cbind(data[,1],scale(data[,-1])*sqrt(ntot)/sqrt(ntot-1))
		data[which(is.na(data))]=0
	intercept=TRUE
			}else{
				message("intercept has been added")
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


num=min(ntot-1,p-1,num+!intercept)
maxqdep=min(maxq,log(min(ntot,p)-1,2))


ntot=nrow(data)
p=ncol(data)


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
ORDREBETA2=numeric(0)
indice=numeric(0)
compteur=0
compteurordre=0
rand_sauv=rand

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
		
		if(showit){print("initialization")}

	a=mht(data2,Y,sigma=q,m=m,maxordre=num,ordre=ordre,var_nonselect=var_nonselect,maxq=maxq,alpha=alpha,show=c(showordre,0,showresult))
	compteurordre=compteurordre+1
	b=a$ordre
	bb=correspondance[2,b]

	if(showit){print(bb)}
		ORDREBETA2=rbind(ORDREBETA2,bb[1:num]) 
		
	#on cherche dans quel espace vit data	
	XI=data2[,a$ordrebeta]
	dec=decompbaseortho(XI)
	#on rajoute dans nonind2 les dernieres variables, celle qui n'ont pas d'utilitÈs puisque dans Rn
	nonind=dec$nonind
	U=dec$U
	if(p>(ntot+length(nonind)))
	{nonind2=c(nonind,(ntot+1+length(nonind)):p)	
		Uchap=U[,-nonind2]
	}else{
		if(length(nonind)==0){Uchap=U}else{Uchap=U[,-nonind]}
			}
	dim_X=ncol(Uchap)			#nombre de variables utiles
	aV=array(0,c(length(alpha),maxqdep,dim_X)) #on y met tous les quantiles
	
	NBR_effect=length(which(a$coefficients!=0))
	aV[,,1:NBR_effect]=array(a$quantile,c(length(alpha),maxqdep,NBR_effect))#}
	assign(paste("aV",compteurordre,sep="_"),aV)
	indice=rbind(var_nonselect-sum(nonind<=var_nonselect),NBR_effect)
	
	ind_sauv2=correspondance[2,a$relevant_var[1,]]
	ind_sauv=ind_sauv2
	if(length(alpha>1))
	{	
		for(alph in 1:length(alpha))
		{indi=correspondance[2,a$relevant_var[alph,]]
		indi2=c(indi,rep(0,length(ind_sauv2)-length(indi)))
		ind_sauv=rbind(ind_sauv,indi2)}
	}
a=matrix(runif(ntot*max(m,2),0,ntot),nrow=ntot)
random=ceiling(a)
		
UCHAP=NULL
y.fitted=NULL
COMPTEUR=NULL
PSI=array(0,c(q,q,length(alpha)))
for(alph in 1:length(alpha)) #boucle sur alpha
{
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
	
	ind=ind_sauv[alph,]
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
aV2=array(0,c(dim(aV)[1],dim(aV)[2],dim(aV)[3],nrow(ORDREBETA2)))
aV2[,,,1:nrow(ORDREBETA2)]=aV
	
mm2=MM2(data2=data2,yhat=yhat,m=m,sigma=sigma_e,maxqdep=maxqdep,num=num,var_nonselect=var_nonselect,random=random,showit=showit,showresult=showresult,correspondance=correspondance,ORDREBETA2=ORDREBETA2,indice=indice,aV2=aV2,choix_ordre=choix_ordre,alpha=alpha,alph=alph,IT=IT)

aV=mm2$aV2
ORDREBETA2=mm2$ORDREBETA2
ind=mm2$ind
indice=mm2$indice
compteurordre=mm2$compteurordre	
}

if(showresult){print(paste("number of selected variables:",length(ind)))}
beta_hat2=matrix(0,p)

if(length(ind)>0){
	reg=lm(yhat~data[,ind]-1)

	beta_hat2[ind]=reg$coefficients
	beta_hat2[-which(beta_hat2!=0)]=0
	sum1=sum((beta_hat2-beta_hat)^2)
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

rownames(aV)=paste("alpha=",alpha)
colnames(aV)=paste("Hk,",0:(maxqdep-1))
dimnames(aV)[[3]]=paste("ktest=",(var_nonselect!=0):(dim_X-(var_nonselect==0)))
dimnames(aV)[[4]]=paste("ordre",1:compteurordre)

colnames(y.fitted)=alpha
rownames(ORDREBETA2)=paste("ordre",1:compteurordre)
colnames(BETA_HAT)=alpha

show=c(showordre,showresult,showit)	
		
out=list(data=list(X=data,Y=Y,z=z,grp=grp),beta=BETA_HAT,fitted.values=y.fitted,u=UCHAP,Psi=PSI,sigma_e=SIGMAE,it=COMPTEUR,quantile=aV,ordrebeta=ORDREBETA2,converge=converge,call=match.call(),arg=list(fix=fix,rand_sauv=rand_sauv,step=step,num=num,ordre=choix_ordre,m=m,random=random,show=show,IT=IT,maxqdep=maxqdep))

out
structure(out,class=c("mhtp","mhtp_nodiag"))

}