##################### likelihood function for the optimization with nlminb ###############
likelihood_nlminb<-function(Q_breslow,D,SigmaInv,X,X_aktuell,Eta_tilde,Betadach,W,n,family)
{
if(is.matrix(X)==FALSE)
{
X<-as.vector(X)
X<-as.matrix(X)
}
Q_breslow<-Q_breslow^2

likeliparam<-1000
mitte<-rep(0,likeliparam)
randli<-rep(0,likeliparam-1)
randre<-rep(0,likeliparam-1)
krit<-TRUE

if(is.null(family$multivariate)){
  V<-(diag(1/(D*SigmaInv*D))+W%*%(t(W)*rep(Q_breslow,n)))
}else{
  V<-(D%*%SigmaInv%*%t(D))+W%*%(t(W)*rep(Q_breslow,n))
}

VVV_tilde<-det(V)
if(VVV_tilde==0 || VVV_tilde==Inf || VVV_tilde==-Inf)
{

max_V_tilde<-max(V)
V_tilde<-V/max_V_tilde
mean_V_tilde<-max_V_tilde

V_tilde1<-V_tilde
V_tilde2<-matrix(c(18,0,0,1),2,2)

if(det(V_tilde)==0 || det(V_tilde)==Inf || det(V_tilde)==-Inf)
{
min_V_tilde<-min(V[V!=0])
V_tilde<-V/min_V_tilde
mean_V_tilde<-min_V_tilde
V_tilde2<-V_tilde
}


if(det(V_tilde1)==0 & det(V_tilde2)==0)
min_V_tilde<-0.001

VVV_tilde<-det(V_tilde)
if(VVV_tilde==0 || VVV_tilde==Inf || VVV_tilde==-Inf)
{
 mitte[1]<-min_V_tilde+0.5*(max_V_tilde-min_V_tilde)
 V_tilde<-V/mitte[1]
 mean_V_tilde<-mitte[1]
}

VVV_tilde<-det(V_tilde)
if(VVV_tilde==0 || VVV_tilde==Inf || VVV_tilde==-Inf)
{
krit<-FALSE

   if(det(V/mitte[1])==det(V/min_V_tilde))
   {
   randli[1]<-mitte[1]
   randre[1]<-max_V_tilde
   mitte[2]<-mitte[1]+0.5*(randre[1]-randli[1])
   }
   if(det(V/mitte[1])==det(V/max_V_tilde))
   {
   randli[1]<-min_V_tilde
   randre[1]<-mitte[1]
   mitte[2]<-mitte[1]-0.5*(randre[1]-randli[1])
   }
   V_tilde<-V/mitte[2]
   mean_V_tilde<-mitte[2]
      VVV_tilde<-det(V_tilde)
      if(VVV_tilde==0 || VVV_tilde==Inf || VVV_tilde==-Inf)
      {
      krit<-FALSE
      }
      else
      {
      krit<-TRUE
      }
}

j<-2
while(!krit & j<likeliparam+1)
{
   if(det(V/mitte[j])==det(V/randli[j-1]))
   {
   randli[j]<-mitte[j]
   randre[j]<-randre[j-1]
   mitte[j+1]<-mitte[j]+0.5*(randre[j]-randli[j])
   }
   if(det(V/mitte[j])==det(V/randre[j-1]))
   {
   randli[j]<-randli[j-1]
   randre[j]<-mitte[j]
   mitte[j+1]<-mitte[j]-0.5*(randre[j]-randli[j])
   }
   V_tilde<-V/mitte[j+1]
   mean_V_tilde<-mitte[j+1]
      VVV_tilde<-det(V_tilde)
      if(VVV_tilde==0 || VVV_tilde==Inf || VVV_tilde==-Inf)
      {
      krit<-FALSE
      }
      else
      {
      krit<-TRUE
      }
j<-j+1
}

likeli<- -0.5*(ncol(V_tilde)*log(mean_V_tilde)+log(VVV_tilde)
               +log(det(t(X_aktuell)%*%chol2inv(chol(V_tilde))%*%X_aktuell))
               +t(Eta_tilde-X%*%Betadach)%*%chol2inv(chol(V))%*%(Eta_tilde-X%*%Betadach))
}else{
likeli<- -0.5*(log(det(V))+log(det(t(X_aktuell)%*%chol2inv(chol(V))%*%X_aktuell))
               +t(Eta_tilde-X%*%Betadach)%*%chol2inv(chol(V))%*%(Eta_tilde-X%*%Betadach))
}

ret.obj<- -likeli
return(ret.obj)
}

###########

likelihood<-function(q_vec,D,SigmaInv,X,X_aktuell,Eta_tilde,Betadach,W,n,s,k,family)
{
if(is.matrix(X)==FALSE)
{
X<-as.vector(X)
X<-as.matrix(X)
}
Q_breslow1<-matrix(0,s,s)
Q_breslow1[lower.tri(Q_breslow1)]<-q_vec[(s+1):(s*(s+1)*0.5)]
Q_breslow1<-Q_breslow1+t(Q_breslow1)
diag(Q_breslow1)<-(q_vec[1:s])

if(all (eigen(Q_breslow1)$values>0))
{
Q_breslow<-matrix(0,n*s,n*s)
for (i in 1:n)
Q_breslow[((i-1)*s+1):(i*s),((i-1)*s+1):(i*s)]<-Q_breslow1

if(is.null(family$multivariate)){
  V<-(diag(1/(D*SigmaInv*D))+W%*%Q_breslow%*%t(W))
}else{
  V<-(D%*%SigmaInv%*%t(D))+W%*%Q_breslow%*%t(W)
}

VVV_tilde<-det(V)
if(VVV_tilde==0 || VVV_tilde==Inf || VVV_tilde==-Inf)
{


determ<-rep(0,n)
determ[1]<-det(V[(1:k[1]),(1:k[1])])

for (i in 2:n)
{
if((sum(k[1:(i-1)])+1)<(sum(k[1:i])))
{
determ[i]<-det(V[(sum(k[1:(i-1)])+1):(sum(k[1:i])),(sum(k[1:(i-1)])+1):(sum(k[1:i]))])
}else{
determ[i]<-V[(sum(k[1:i])),(sum(k[1:i]))]
}}

norm_determ<-min(determ[determ!=0])
determ_gross<-prod((1/norm_determ)*determ)
deti_max<-determ_gross
##########
mitte<-rep(0,1000)
randli<-rep(0,999)
randre<-rep(0,999)
deti<-rep(0,1000)
krit<-TRUE

if (determ_gross==0 || determ_gross==Inf || determ_gross==-Inf)
{
norm_determ<-max(determ)
determ_gross<-prod((1/norm_determ)*determ)
deti_min<-determ_gross
}


if (determ_gross==0 || determ_gross==Inf || determ_gross==-Inf)
{
mitte[1]<-min(determ)+0.5*(max(determ)-min(determ))
norm_determ<-mitte[1]
deti[1]<-prod((1/norm_determ)*determ)
determ_gross<-deti[1]
}



if(determ_gross==0 || determ_gross==Inf || determ_gross==-Inf)
{
krit<-FALSE

   if(deti[1]==deti_min)
   {
   randli[1]<-min(determ)
   randre[1]<-mitte[1]
   mitte[2]<-mitte[1]-0.5*(randre[1]-randli[1])
   deti_min<-deti[1]
   }
   if(deti[1]==deti_max)
   {
   randli[1]<-mitte[1]
   randre[1]<-max(determ)
   mitte[2]<-mitte[1]+0.5*(randre[1]-randli[1])
   deti_max<-deti[1]
   }
   norm_determ<-mitte[2]
   deti[2]<-prod((1/norm_determ)*determ)
   determ_gross<-deti[2]
   if(deti[2]==0 || deti[2]==Inf || deti[2]==-Inf)
      {
      krit<-FALSE
      }else{
      krit<-TRUE
      }
}

j<-2
while(!krit & j<1001)
{
  if(deti[j]==deti_min)
   {
   randli[j]<-randli[j-1]
   randre[j]<-mitte[j]
   mitte[j+1]<-mitte[j]-0.5*(randre[j]-randli[j])
   deti_min<-deti[j]
   }
   if(deti[j]==deti_max)
   {
   randli[j]<-mitte[j]
   randre[j]<-randre[j-1]
   mitte[j+1]<-mitte[j]+0.5*(randre[j]-randli[j])
   deti_max<-deti[j]
   }
   norm_determ<-mitte[j+1]
   deti[j+1]<-prod((1/norm_determ)*determ)
   determ_gross<-deti[j+1]
   if(deti[j+1]==0 || deti[j+1]==Inf || deti[j+1]==-Inf)
      {
      krit<-FALSE
      }else{
      krit<-TRUE
      }
j<-j+1
}

likeli<- -0.5*(n*log(norm_determ)+log(determ_gross)+(log(det(t(X_aktuell)%*%chol2inv(chol(V))%*%X_aktuell))+t(Eta_tilde-X%*%Betadach)%*%chol2inv(chol(V))%*%(Eta_tilde-X%*%Betadach)))
}else{
likeli<- -0.5*(log(det(V))+log(det(t(X_aktuell)%*%chol2inv(chol(V))%*%X_aktuell))+t(Eta_tilde-X%*%Betadach)%*%chol2inv(chol(V))%*%(Eta_tilde-X%*%Betadach))
}}else{
likeli<- -1e+20
}
ret.obj<- -likeli
ret.obj
return(ret.obj)
}







########### likelihood for diagonal


likelihood_diag<-function(q_vec,D,SigmaInv,X,X_aktuell,Eta_tilde,Betadach,W,n,s,k,rnd.len,family)
{
if(is.matrix(X)==FALSE)
{
X<-as.vector(X)
X<-as.matrix(X)
}
Q_breslow1<-diag(q_vec)

if(all (eigen(Q_breslow1)$values>0))
{
Q_breslow<-matrix(0,n%*%s,n%*%s)



for(i in 1:n[1])
Q_breslow[((i-1)*s[1]+1):(i*s[1]),((i-1)*s[1]+1):(i*s[1])]<-Q_breslow1[1:s[1],1:s[1]]

for(ie in 2:rnd.len)
{
for(i in 1:n[ie])
Q_breslow[((i-1)*s[ie]+1):(i*s[ie]),((i-1)*s[ie]+1):(i*s[ie])]<-Q_breslow1[(sum(s[1:(ie-1)])+1):sum(s[1:ie]),(sum(s[1:(ie-1)])+1):sum(s[1:ie])]
}
                                                                                                                                            
if(is.null(family$multivariate)){
  V<-(diag(1/(D*SigmaInv*D))+W%*%Q_breslow%*%t(W))
}else{
  V<-(D%*%SigmaInv%*%t(D))+W%*%Q_breslow%*%t(W)
}


VVV_tilde<-det(V)
if(VVV_tilde==0 || VVV_tilde==Inf || VVV_tilde==-Inf)
{


determ<-rep(0,sum(n))
determ[1]<-det(V[(1:k[1]),(1:k[1])])

for (i in 2:sum(n))
{
if((sum(k[1:(i-1)])+1)<(sum(k[1:i])))
{
determ[i]<-det(V[(sum(k[1:(i-1)])+1):(sum(k[1:i])),(sum(k[1:(i-1)])+1):(sum(k[1:i]))])
}else{
determ[i]<-V[(sum(k[1:i])),(sum(k[1:i]))]
}}

norm_determ<-min(determ[determ!=0])
determ_gross<-prod((1/norm_determ)*determ)
deti_max<-determ_gross
##########
mitte<-rep(0,1000)
randli<-rep(0,999)
randre<-rep(0,999)
deti<-rep(0,1000)
krit<-TRUE

if (determ_gross==0 || determ_gross==Inf || determ_gross==-Inf)
{
norm_determ<-max(determ)
determ_gross<-prod((1/norm_determ)*determ)
deti_min<-determ_gross
}


if (determ_gross==0 || determ_gross==Inf || determ_gross==-Inf)
{
mitte[1]<-min(determ)+0.5*(max(determ)-min(determ))
norm_determ<-mitte[1]
deti[1]<-prod((1/norm_determ)*determ)
determ_gross<-deti[1]
}



if(determ_gross==0 || determ_gross==Inf || determ_gross==-Inf)
{
krit<-FALSE

   if(deti[1]==deti_min)
   {
   randli[1]<-min(determ)
   randre[1]<-mitte[1]
   mitte[2]<-mitte[1]-0.5*(randre[1]-randli[1])
   deti_min<-deti[1]
   }
   if(deti[1]==deti_max)
   {
   randli[1]<-mitte[1]
   randre[1]<-max(determ)
   mitte[2]<-mitte[1]+0.5*(randre[1]-randli[1])
   deti_max<-deti[1]
   }
   norm_determ<-mitte[2]
   deti[2]<-prod((1/norm_determ)*determ)
   determ_gross<-deti[2]
   if(deti[2]==0 || deti[2]==Inf || deti[2]==-Inf)
      {
      krit<-FALSE
      }else{
      krit<-TRUE
      }
}

j<-2
while(!krit & j<1001)
{
  if(deti[j]==deti_min)
   {
   randli[j]<-randli[j-1]
   randre[j]<-mitte[j]
   mitte[j+1]<-mitte[j]-0.5*(randre[j]-randli[j])
   deti_min<-deti[j]
   }
   if(deti[j]==deti_max)
   {
   randli[j]<-mitte[j]
   randre[j]<-randre[j-1]
   mitte[j+1]<-mitte[j]+0.5*(randre[j]-randli[j])
   deti_max<-deti[j]
   }
   norm_determ<-mitte[j+1]
   deti[j+1]<-prod((1/norm_determ)*determ)
   determ_gross<-deti[j+1]
   if(deti[j+1]==0 || deti[j+1]==Inf || deti[j+1]==-Inf)
      {
      krit<-FALSE
      }else{
      krit<-TRUE
      }
j<-j+1
}

likeli<- -0.5*(n*log(norm_determ)+log(determ_gross)+(log(det(t(X_aktuell)%*%chol2inv(chol(V))%*%X_aktuell))+t(Eta_tilde-X%*%Betadach)%*%chol2inv(chol(V))%*%(Eta_tilde-X%*%Betadach)))
}else{
likeli<- -0.5*(log(det(V))+log(det(t(X_aktuell)%*%chol2inv(chol(V))%*%X_aktuell))+t(Eta_tilde-X%*%Betadach)%*%chol2inv(chol(V))%*%(Eta_tilde-X%*%Betadach))
}}else{
likeli<- -1e+20
}
ret.obj<- -likeli
ret.obj
return(ret.obj)
}




###########    likelihood for block

likelihood_block<-function(q_vec,D,SigmaInv,X,X_aktuell,Eta_tilde,Betadach,W,n,s,k,rnd.len,family)
{
if(is.matrix(X)==FALSE)
{
X<-as.vector(X)
X<-as.matrix(X)
}
Q_breslow1<-list()
eig.value<-rep(FALSE,rnd.len)


Q_breslow1[[1]]<-matrix(0,s[1],s[1])
if(s[1]>1)
Q_breslow1[[1]][lower.tri(Q_breslow1[[1]])]<-q_vec[(s[1]+1):(s[1]*(s[1]+1)*0.5)]
Q_breslow1[[1]]<-Q_breslow1[[1]]+t(Q_breslow1[[1]])
diag(Q_breslow1[[1]])<-(q_vec[1:s[1]])
eig.value[rnd.len]<-all (eigen(Q_breslow1[[1]])$values>0)
q_vec<-q_vec[-c(1:(s[1]*(s[1]+1)*0.5))]

for (zu in 2:rnd.len)
{
Q_breslow1[[zu]]<-matrix(0,s[zu],s[zu])
if(s[zu]>1)
Q_breslow1[[zu]][lower.tri(Q_breslow1[[zu]])]<-q_vec[(s[zu]+1):(s[zu]*(s[zu]+1)*0.5)]
Q_breslow1[[zu]]<-Q_breslow1[[zu]]+t(Q_breslow1[[zu]])
diag(Q_breslow1[[zu]])<-(q_vec[1:s[zu]])
eig.value[rnd.len]<-all (eigen(Q_breslow1[[zu]])$values>0)
q_vec<-q_vec[-c(1:(s[zu]*(s[zu]+1)*0.5))]
}

if(all (eig.value))
{

Q_breslow<-matrix(0,n%*%s,n%*%s)
for(i in 1:n[1])
Q_breslow[((i-1)*s[1]+1):(i*s[1]),((i-1)*s[1]+1):(i*s[1])]<-Q_breslow1[1:s[1],1:s[1]]

for(ie in 2:rnd.len)
{
for(i in 1:n[ie])
Q_breslow[((i-1)*s[ie]+1):(i*s[ie]),((i-1)*s[ie]+1):(i*s[ie])]<-Q_breslow1[(sum(s[1:(ie-1)])+1):sum(s[1:ie]),(sum(s[1:(ie-1)])+1):sum(s[1:ie])]
}


if(is.null(family$multivariate)){
  V<-(diag(1/(D*SigmaInv*D))+W%*%Q_breslow%*%t(W))
}else{
  V<-(D%*%SigmaInv%*%t(D))+W%*%Q_breslow%*%t(W)
}

VVV_tilde<-det(V)
if(VVV_tilde==0 || VVV_tilde==Inf || VVV_tilde==-Inf)
{


determ<-rep(0,sum(n))
determ[1]<-det(V[(1:k[1]),(1:k[1])])

for (i in 2:sum(n))
{
if((sum(k[1:(i-1)])+1)<(sum(k[1:i])))
{
determ[i]<-det(V[(sum(k[1:(i-1)])+1):(sum(k[1:i])),(sum(k[1:(i-1)])+1):(sum(k[1:i]))])
}else{
determ[i]<-V[(sum(k[1:i])),(sum(k[1:i]))]
}}

norm_determ<-min(determ[determ!=0])
determ_gross<-prod((1/norm_determ)*determ)
deti_max<-determ_gross
##########
mitte<-rep(0,1000)
randli<-rep(0,999)
randre<-rep(0,999)
deti<-rep(0,1000)
krit<-TRUE

if (determ_gross==0 || determ_gross==Inf || determ_gross==-Inf)
{
norm_determ<-max(determ)
determ_gross<-prod((1/norm_determ)*determ)
deti_min<-determ_gross
}


if (determ_gross==0 || determ_gross==Inf || determ_gross==-Inf)
{
mitte[1]<-min(determ)+0.5*(max(determ)-min(determ))
norm_determ<-mitte[1]
deti[1]<-prod((1/norm_determ)*determ)
determ_gross<-deti[1]
}



if(determ_gross==0 || determ_gross==Inf || determ_gross==-Inf)
{
krit<-FALSE

   if(deti[1]==deti_min)
   {
   randli[1]<-min(determ)
   randre[1]<-mitte[1]
   mitte[2]<-mitte[1]-0.5*(randre[1]-randli[1])
   deti_min<-deti[1]
   }
   if(deti[1]==deti_max)
   {
   randli[1]<-mitte[1]
   randre[1]<-max(determ)
   mitte[2]<-mitte[1]+0.5*(randre[1]-randli[1])
   deti_max<-deti[1]
   }
   norm_determ<-mitte[2]
   deti[2]<-prod((1/norm_determ)*determ)
   determ_gross<-deti[2]
   if(deti[2]==0 || deti[2]==Inf || deti[2]==-Inf)
      {
      krit<-FALSE
      }else{
      krit<-TRUE
      }
}

j<-2
while(!krit & j<1001)
{
  if(deti[j]==deti_min)
   {
   randli[j]<-randli[j-1]
   randre[j]<-mitte[j]
   mitte[j+1]<-mitte[j]-0.5*(randre[j]-randli[j])
   deti_min<-deti[j]
   }
   if(deti[j]==deti_max)
   {
   randli[j]<-mitte[j]
   randre[j]<-randre[j-1]
   mitte[j+1]<-mitte[j]+0.5*(randre[j]-randli[j])
   deti_max<-deti[j]
   }
   norm_determ<-mitte[j+1]
   deti[j+1]<-prod((1/norm_determ)*determ)
   determ_gross<-deti[j+1]
   if(deti[j+1]==0 || deti[j+1]==Inf || deti[j+1]==-Inf)
      {
      krit<-FALSE
      }else{
      krit<-TRUE
      }
j<-j+1
}

likeli<- -0.5*(n*log(norm_determ)+log(determ_gross)+(log(det(t(X_aktuell)%*%chol2inv(chol(V))%*%X_aktuell))+t(Eta_tilde-X%*%Betadach)%*%chol2inv(chol(V))%*%(Eta_tilde-X%*%Betadach)))
}else{
likeli<- -0.5*(log(det(V))+log(det(t(X_aktuell)%*%chol2inv(chol(V))%*%X_aktuell))+t(Eta_tilde-X%*%Betadach)%*%chol2inv(chol(V))%*%(Eta_tilde-X%*%Betadach))
}}else{
likeli<- -1e+20
}
ret.obj<- -likeli
ret.obj
return(ret.obj)
}

