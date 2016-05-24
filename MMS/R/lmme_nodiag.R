lmme_nodiag=function(data,Y,z,grp,step,showit)
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

#on initialise tout
Psi=matrix(0.02,nrow=q,ncol=q)
diag(Psi)=1:q
Psi1=solve(Psi)
#on met tous ca dans une matrice Psi de taille q*q
#on prend tous les groupes =Ni


############################
# alg51bolm
###########################

	sigma_e=1
	reg=lm(Y~data-1)
	beta_hat2=reg$coefficients
	beta_hat2[-which(beta_hat2!=0)]=0
	beta_hat=beta_hat2
	sum1=10
	condition=1
	ft=numeric(0)
	uchap=NULL
	
	compteur_ijk=0
	fit=0
	converge=TRUE
	delete=1

Z=numeric(0)
for(k in 1:q)
{Z=cbind(Z,get(paste("Z",k,sep="_")))}
N=ncol(Z)

	while(((sum1>10^-8)||(sum(sum2)>10^-8))||(abs(condition)>10^-6))
	{
		compteur_ijk=compteur_ijk+1
		if(showit){print(paste("iteration=",compteur_ijk))}
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

	reg=lm(yhat~data-1)
	beta_hat2=reg$coefficients
	beta_hat2[-which(beta_hat2!=0)]=0
	sum1=sum((beta_hat2-beta_hat)^2)
	#print(sum1)
	beta_hat=beta_hat2
##############################################################################################################
#################################### E-step -2 ##########################################################
##############################################################################################################

#on calcule les nouveaux uchapal[t+1] a partir de uchap[t]
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
	Psi[j,i]=Psi[i,j]
}}
Psi1=solve(Psi)
	
		
#on calcule sigma_e sur les résidus du modele
sumk=Z%*%uchap
yhatt=data%*%beta_hat+sumk

aa=Tkk%*%G1 #(T*G-1)
a=N-sigma_e*sum(diag(aa))
sigma_e=sum((Y-yhatt)^2)/ntot+(a*sigma_e/ntot)


if(showit)
{print("Psi")
	print(Psi)
print(paste("sigma_e",sigma_e))}

if(compteur_ijk>step){message("Algorithm did not converge..")
	converge=FALSE
break}


suma=0
sumk=0
sumj=0
for(k in 1:q)
{sumk=sumk+get(paste("Z",k,sep="_"))%*%get(paste("uchap",k,sep="_"))}
suma=Ni*log(det(Psi))+t(uchap)%*%G1%*%uchap

yhatt=data%*%beta_hat+sumk
#log(det(G))+u'G-1u

ft=c(ft,ntot*log(sigma_e)+suma+sum((Y-yhatt)^2)/sigma_e)
if(length(ft)>2){condition=ft[length(ft)]-ft[length(ft)-1]}else{condition=10}

}#fin while


names(beta_hat)=c("Intercept",paste("X", 2:p,sep=""))

#########################
# fin alg51bolm
#########################
out=list(data=list(X=data,Y=Y,z=z,grp=grp),beta=beta_hat,Psi=Psi,sigma_e=sigma_e,fitted.values=yhatt,it=compteur_ijk,converge=converge,u=uchap,call=match.call())

out
structure(out,class="lmme")


}#fin function
