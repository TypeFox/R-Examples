MM2=function(data2,yhat,m,sigma,maxqdep,num,var_nonselect,random,showit,showresult,correspondance,ORDREBETA2,indice,aV2,choix_ordre,alpha,alph,IT)
{#####################
# MM2 commence ###
#####################
compteurordre=nrow(ORDREBETA2)
ntot=nrow(data2)
p=ncol(data2)
#on ordonne les variables pour avoir les var_nonselect_fixed en premier et ensuite les var_nonselect_random!=0, pour pouvoir faire les test correctement sans selection var_nonselect=sum des deux

if(choix_ordre=="pval")
{
	if(p<ntot) 	#on calcule pval avec une seule regression comportant toutes les variables
	{reg=lm(yhat~data2-1)
		a=summary(reg)$coefficients[,4]
		a[1:var_nonselect]=0	#on ne selectionne pas l'intercept
		b=sort(a,index.return=TRUE)
		ORDREBETA=b$ix[1:num]
		XI_ord=data2[,b$ix] #on a ainsi les XI ordonnÈes
		if(showit){print(b$ix)}
		}else{		#on calcule pval avec FDR2: une regression pour chaque variable
			print("p>n, the order 'pval' is not possible, 'pval_hd' is used instead")
			choix_ordre="pval_hd"
			}
}

if(choix_ordre=="pval_hd")
{a=numeric(0)
	for(i in 1:p)
	{
		reg=lm(yhat~data2[,i]-1)
		c=summary(reg)$coefficients[,4]
		a=c(a,c)
	}
	a[1:var_nonselect]=0	#on ne selectionne pas l'intercept
	b=sort(a,index.return=TRUE)
		ORDREBETA=b$ix[1:num]
	XI_ord=data2[,b$ix] #on a ainsi les XI ordonnÈes
	if(showit){print(b$ix)}
}

if(choix_ordre=="bolasso")
{
ordre=dyadiqueordre(data2,yhat,m,maxordre=num,var_nonselect=var_nonselect,showtest=FALSE,showordre=FALSE,random=random)# donne l'ordre (dans ordre) et le nombre de fois ou l'algo a redemarré (dans prob)	
	b=ordre$ordre
bb=correspondance[2,b]
#mean that data2[ordre]=data[correspondance[2,ordre]]=data[,bb]
if(showit){print(bb)}
	ORDREBETA=bb[1:num]
#on complete l'ordre par les variables restantes
	a=match(1:p,b)
	b=c(b,(1:p)[which(is.na(a))])
	
XI_ord=data2[,b] #on a ainsi les XI ordonnées dans XI_ord
}

if(choix_ordre=="FR")
{
    ordre=1
    ind.left=2:p
    while(length(ordre)<min(ntot-1,p,num))
    {
        res=NULL
        for(i in 1:length(ind.left))
        {
            X=data2[,c(ordre,ind.left[i])]
            I=diag(1,ntot)
            H=X%*%solve(t(X)%*%X)%*%t(X)
            res[i]=t(yhat)%*%(I-H)%*%yhat	#same as reg=lm(Y~XI[,VRAINDICE]-1); sum(reg$residuals^2)
        }
        a=which.min(res)
        ordre=c(ordre,ind.left[a])
        ind.left=ind.left[-a]
        
    }
    b=ORDRE=ordre
    for(i in 1:p)
    {if(sum(i==b)==0){b=c(b,i)}} #on complete l'ordre par les variables restantes
    ORDREBETA=matrix(b,nrow=1)
    XI_ord=data2[,b] #on a ainsi les XI ordonnées dans XI_ord
	if(showit){print(ordre)}
}


dec=decompbaseortho(XI_ord)
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


beta2=lm(yhat~Uchap-1)$coefficients #decomposition de Y2 dans la base orthonormal (X(1),..,X(p))
beta2[-which(beta2!=0)]=0

indice2=var_nonselect-1-sum(nonind<=var_nonselect)

ktest=var_nonselect-sum(nonind<=var_nonselect)
nbr_test=var_nonselect

calcul=numeric(0)
maxq=min(log(min(ntot,dim_X)-1,2),maxqdep)
aV=array(0,c(length(alpha),maxq,dim_X)) #on y met tous les quantiles
	T=1
	while((T>0)&&(dim_X>ktest+1))
		{
		maxq=min(log(min(ntot,dim_X)-ktest-1,2),maxqdep)
	#on est a ktest fixé, on regarde dans tout ce qu'on a fait avant ou indice>ktest si on a deja calculé le quantile

	#on a indice de longueur compteur-1, aV_compteur avec les quantiles, var_select de dim compteur-1, et ORDREBETA2 
if(ktest>indice2)
{abc=compteurordre#dim de indice/ORDREBETA2/var_select
	i=0
	I=numeric(0)
	TT=0
	#A=numeric(0)
	while((TT==0)&&(i<abc))#dim(ORDREBETA2)=abc*num
	{i=i+1
	a=numeric(0)
	K=numeric(0)
	for(j in 1:nbr_test)
	{a=c(a,sum(ORDREBETA[j]==ORDREBETA2[i,1:nbr_test]))}
	
	if((sum(a)==nbr_test)&&(indice[2,i]>=ktest)&&(indice[1,i])<=ktest){TT=1
		I=c(I,i)}
	}

	if(length(I)!=0)
	{	aV[,,ktest]=aV2[,,ktest,I]
		calcul=c(calcul,0)#on met 0 si le quantile est deja calculé
		}else{
		if(showresult){print(paste("ktest=",ktest))}
		quant=quantilemht(XI_ord,ktest,alpha,IT=IT,maxq=maxq,sigma=sigma)
		aV[,1:maxq,ktest]=quant$quantile	
		calcul=c(calcul,1)#on met 1 si on a calculé un quantile manquant
		if(showresult){print(aV[,,ktest])}
		}
indice2=indice2+1
}
		#test S	
		bb=numeric(0)
	
		for(m in 0:(maxq-1))
			{
			a=sum(beta2[(ktest+1):(ktest+2^m)]^2)/sigma
			b=a>aV[alph,m+1,ktest]#F #1 si on doit rejeter le test, 0 sinon
			bb=c(bb,b) #on met tous les tests de Hk
			}
		if(length(which(bb!=0))>0) #securitÈ a la base, inutile maintenant? a verifier
			{bb[-which(bb!=0)]=0}else{bb=matrix(0,1,m)}
		
		if(sum(bb)>0){T=1
		ktest=ktest+1
			
		nbr_test=nbr_test+1
		if(length(nonind)>0)
		{for(i in 1:length(nonind))
			if(sum(nonind==(nbr_test+1))==1){nbr_test=nbr_test+1}}
											
		}else{
			T=0
			#print(ktest)
		} #on rejete Hk s'il y a au moins un test qui rejete
	}#fin while
if(ktest==dim_X){k0=dim_X}else{k0=ktest}

NBR=nbr_test #rÈsultat contenant le nombre de variables sÈlectionnÈes
NBR_effect=k0


if(sum(calcul)!=0)
{#on rajoute que si on a calculer au moins un quantile
	aV3=array(0,c(dim(aV2)[1],dim(aV2)[2],dim(aV2)[3],nrow(ORDREBETA2)+1))
	aV3[,,,1:nrow(ORDREBETA2)]=aV2
	compteurordre=compteurordre+1
	ORDREBETA2=rbind(ORDREBETA2,ORDREBETA)
	aV3[,,,nrow(ORDREBETA2)]=aV
	aV2=aV3

	indice=cbind(indice,c(var_nonselect-sum(nonind<=var_nonselect),NBR_effect))#c(indice,NBR)
	
	}	
	ind=ORDREBETA[1:NBR]

out=list(ind=ind,indice=indice,ORDREBETA2=ORDREBETA2,aV2=aV2,compteurordre=compteurordre,ORDREBETA=ORDREBETA,aV=aV)
###########################
#  fin MM2 ###
##########################
}