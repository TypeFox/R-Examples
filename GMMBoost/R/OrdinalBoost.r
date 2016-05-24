##################### likelihood function for the optimization with bobyqa ###############
likelihood_bobyqa_cumul<-function(Q_breslow,D,Sigma,X,X_aktuell,Eta_tilde,Betadach,W,n)
{
Q_breslow<-Q_breslow^2

likeliparam<-1000
mitte<-rep(0,likeliparam)
randli<-rep(0,likeliparam-1)
randre<-rep(0,likeliparam-1)
krit<-TRUE

V<-chol2inv(chol(t(t(chol2inv(chol(Sigma))*D)*D)))+W%*%(t(W)*rep(Q_breslow,n))

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

likeli<- -0.5*(ncol(V_tilde)*log(mean_V_tilde)+log(VVV_tilde)+log(det(t(X_aktuell)%*%chol2inv(chol(V_tilde))%*%X_aktuell))+t(Eta_tilde-X%*%Betadach)%*%chol2inv(chol(V))%*%(Eta_tilde-X%*%Betadach))
}else{
likeli<- -0.5*(log(det(V))+log(det(t(X_aktuell)%*%chol2inv(chol(V))%*%X_aktuell))+t(Eta_tilde-X%*%Betadach)%*%chol2inv(chol(V))%*%(Eta_tilde-X%*%Betadach))
}

ret.obj<- -likeli
return(ret.obj)
}

##################### likelihood function for the optimization with bobyqa ###############
likelihood_bobyqa_seq<-function(Q_breslow,D,Sigma,X,X_aktuell,Eta_tilde,Betadach,W,n)
{
  Q_breslow<-Q_breslow^2
  
  likeliparam<-1000
  mitte<-rep(0,likeliparam)
  randli<-rep(0,likeliparam-1)
  randre<-rep(0,likeliparam-1)
  krit<-TRUE
  
  V<-(diag(1/(D*1/Sigma*D))+W%*%(t(W)*rep(Q_breslow,n)))
    
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
    
    likeli<- -0.5*(ncol(V_tilde)*log(mean_V_tilde)+log(VVV_tilde)+log(det(t(X_aktuell)%*%chol2inv(chol(V_tilde))%*%X_aktuell))+t(Eta_tilde-X%*%Betadach)%*%chol2inv(chol(V))%*%(Eta_tilde-X%*%Betadach))
  }else{
    likeli<- -0.5*(log(det(V))+log(det(t(X_aktuell)%*%chol2inv(chol(V))%*%X_aktuell))+t(Eta_tilde-X%*%Betadach)%*%chol2inv(chol(V))%*%(Eta_tilde-X%*%Betadach))
  }
  
  ret.obj<- -likeli
  return(ret.obj)
}

###########

likelihood_cumul<-function(q_vec,D,Sigma,X,X_aktuell,Eta_tilde,Betadach,W,n,s,k)
{
Q_breslow1<-matrix(0,s,s)
Q_breslow1[lower.tri(Q_breslow1)]<-q_vec[(s+1):(s*(s+1)*0.5)]
Q_breslow1<-Q_breslow1+t(Q_breslow1)
diag(Q_breslow1)<-(q_vec[1:s])

if(all (eigen(Q_breslow1)$values>0))
{
Q_breslow<-matrix(0,n*s,n*s)
for (i in 1:n)
Q_breslow[((i-1)*s+1):(i*s),((i-1)*s+1):(i*s)]<-Q_breslow1

V<-chol2inv(chol(t(t(chol2inv(chol(Sigma))*D)*D)))+W%*%Q_breslow%*%t(W)


VVV_tilde<-det(V)
if(VVV_tilde==0 || VVV_tilde==Inf || VVV_tilde==-Inf)
{


determ<-rep(0,n)
determ[1]<-det(V[(1:k[1]),(1:k[1])])

for (i in 2:n)
determ[i]<-det(V[(sum(k[1:(i-1)])+1):(sum(k[1:i])),(sum(k[1:(i-1)])+1):(sum(k[1:i]))])


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


###########

likelihood_seq<-function(q_vec,D,Sigma,X,X_aktuell,Eta_tilde,Betadach,W,n,s,k)
{
  Q_breslow1<-matrix(0,s,s)
  Q_breslow1[lower.tri(Q_breslow1)]<-q_vec[(s+1):(s*(s+1)*0.5)]
  Q_breslow1<-Q_breslow1+t(Q_breslow1)
  diag(Q_breslow1)<-(q_vec[1:s])
  
  if(all (eigen(Q_breslow1)$values>0))
  {
    Q_breslow<-matrix(0,n*s,n*s)
    for (i in 1:n)
      Q_breslow[((i-1)*s+1):(i*s),((i-1)*s+1):(i*s)]<-Q_breslow1
    
      V<-(diag(1/(D*1/Sigma*D))+W%*%Q_breslow%*%t(W))
    
    
    VVV_tilde<-det(V)
    if(VVV_tilde==0 || VVV_tilde==Inf || VVV_tilde==-Inf)
    {
      
      
      determ<-rep(0,n)
      determ[1]<-det(V[(1:k[1]),(1:k[1])])
      
      for (i in 2:n)
        determ[i]<-det(V[(sum(k[1:(i-1)])+1):(sum(k[1:i])),(sum(k[1:(i-1)])+1):(sum(k[1:i]))])
      
      
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

#################################################################################################################################

sequentialize<-function(y,X,katvar)
{
lev<-as.numeric(levels(as.factor(y)))

kat<-length(lev)

if(!is.null(katvar))
{
X_kat<-cbind(1,X[,katvar])
colnames(X_kat)<-c("Icept",katvar)
namesvec<-colnames(X)[!is.element(colnames(X),katvar)]
X<-X[,namesvec]
}

y_bin<-numeric()
X_neu<-numeric()

for (kf in 1:length(y))
{
if(is.null(katvar))
{
if(as.factor(y[kf])==lev[1])
{
y_bin<-c(y_bin,1)
X_neu<-rbind(X_neu,c(1,rep(0,kat-2),X[kf,]))
}

for (ij in 2:(kat-1))
{
if(as.factor(y[kf])==lev[ij])
y_bin<-c(y_bin,rep(0,(ij-1)),1)
}

if(kat>4)
{
for (ij in 2:(kat-3))
{
if(as.factor(y[kf])==lev[ij])
X_neu<-rbind(X_neu,cbind(diag(ij),matrix(0,ij,kat-ij-1),t(matrix(X[kf,],dim(X)[2],ij))))
}}

if(kat>3)
{
if(as.factor(y[kf])==lev[kat-2])
X_neu<-rbind(X_neu,cbind(diag(kat-2),rep(0,kat-2),t(matrix(X[kf,],dim(X)[2],kat-2))))
}

if(as.factor(y[kf])==lev[kat-1])
X_neu<-rbind(X_neu,cbind(diag(kat-1),t(matrix(X[kf,],dim(X)[2],kat-1))))

if(as.factor(y[kf])==lev[kat])
{
y_bin<-c(y_bin,rep(0,(kat-1)))
X_neu<-rbind(X_neu,cbind(diag(kat-1),t(matrix(X[kf,],dim(X)[2],kat-1))))
}
}else{
if(as.factor(y[kf])==lev[1])
{
y_bin<-c(y_bin,1)
X_neu<-rbind(X_neu,c(X_kat[kf,],rep(0,(kat-2)*dim(X_kat)[2]),X[kf,]))
}

for (ij in 2:(kat-1))
{
if(as.factor(y[kf])==lev[ij])
y_bin<-c(y_bin,rep(0,(ij-1)),1)
}

if(kat>4)
{
for (ij in 2:(kat-3))
{
if(as.factor(y[kf])==lev[ij])
{
X_kat_akt<-adiag(t(X_kat[kf,]),t(X_kat[kf,]))
if (ij>2)
{
for (jf in 3:ij)
X_kat_akt<-adiag(X_kat_akt,t(X_kat[kf,]))
}
X_neu<-rbind(X_neu,cbind(X_kat_akt,matrix(0,ij,(kat-ij-1)*dim(X_kat)[2]),t(matrix(X[kf,],dim(X)[2],ij))))
}
}}

if(kat>3)
{
if(as.factor(y[kf])==lev[kat-2])
{
X_kat_akt<-adiag(t(X_kat[kf,]),t(X_kat[kf,]))

if(kat>4)
{
for (jf in 3:(kat-2))
X_kat_akt<-adiag(X_kat_akt,t(X_kat[kf,]))
}
X_neu<-rbind(X_neu,cbind(X_kat_akt,matrix(0,kat-2,dim(X_kat)[2]),t(matrix(X[kf,],dim(X)[2],kat-2))))
}}

if(as.factor(y[kf])==lev[kat-1])
{
X_kat_akt<-adiag(t(X_kat[kf,]),t(X_kat[kf,]))

if(kat>3)
{
for (jf in 3:(kat-1))
X_kat_akt<-adiag(X_kat_akt,t(X_kat[kf,]))
}
X_neu<-rbind(X_neu,cbind(X_kat_akt,t(matrix(X[kf,],dim(X)[2],kat-1))))
}

if(as.factor(y[kf])==lev[kat])
{
y_bin<-c(y_bin,rep(0,(kat-1)))
X_kat_akt<-adiag(t(X_kat[kf,]),t(X_kat[kf,]))

if(kat>3)
{
for (jf in 3:(kat-1))
{
X_kat_akt<-adiag(X_kat_akt,t(X_kat[kf,]))
}}
X_neu<-rbind(X_neu,cbind(X_kat_akt,t(matrix(X[kf,],dim(X)[2],kat-1))))
}

}}
Data_new<-cbind(y_bin,X_neu)
if(is.null(katvar))
{colnames(X_neu)[1:(kat-1)]<-paste("Icept",c(1:(kat-1)),sep=".")
}else{
colnames(X_neu)[1:(((kat-1)*dim(X_kat)[2]))]<-paste(colnames(X_kat),rep(1:(kat-1), each = dim(X_kat)[2]),sep=".")
colnames(X_neu)[(((kat-1)*dim(X_kat)[2])+1):dim(X_neu)[2]]<-namesvec
}
list(y_bin=y_bin,X_new=X_neu,kat=kat)
}
#################################################################################################################################

kumulize<-function(y,X,katvar)
{
lev<-as.numeric(levels(as.factor(y)))

kat<-length(lev)

if(!is.null(katvar))
{
X_kat<-cbind(1,X[,katvar])
colnames(X_kat)<-c("Icept",katvar)
namesvec<-colnames(X)[!is.element(colnames(X),katvar)]
X<-X[,namesvec]
}

y_bin<-numeric()
X_neu<-numeric()

for (kf in 1:length(y))
{
   for (ih in 1:kat)
   {
   if(y[kf]==lev[ih])
   y_bin<-c(y_bin,rep(0,ih-1),rep(1,kat-ih))
}}

if(is.null(katvar))
{
   for (kf in 1:length(y))
   {
   XX<-numeric()
   for (i in 1:(kat-1))
   XX<-rbind(XX,X[kf,])
   X_neu<-rbind(X_neu,cbind(diag(kat-1),XX))
   }
}else{
   for (kf in 1:length(y))
   {
   XX<-numeric()
   for (i in 1:(kat-1))
   XX<-rbind(XX,X[kf,])
   X_kat_akt<-t(X_kat[kf,])
   if((kat-1)>1)
   {
   for (vf in 2:(kat-1))
   X_kat_akt<-adiag(X_kat_akt,t(X_kat[kf,]))
   }
   X_neu<-rbind(X_neu,cbind(X_kat_akt,XX))
   }
}
if(is.null(katvar))
{colnames(X_neu)[1:(kat-1)]<-paste("Icept",c(1:(kat-1)),sep=".")
}else{
colnames(X_neu)[1:(((kat-1)*dim(X_kat)[2]))]<-paste(colnames(X_kat),rep(1:(kat-1), each = dim(X_kat)[2]),sep=".")
colnames(X_neu)[(((kat-1)*dim(X_kat)[2])+1):dim(X_neu)[2]]<-namesvec
}
list(y_bin=y_bin,X_new=X_neu,kat=kat)
}

###############################################

glmm_final_cumul<-function(y,X,W,k,q_start,Delta_start,s,steps=1000,
                           method,kat,katvar,print.iter.final=FALSE,eps.final=1e-5)
{
family=binomial(link="logit")
N<-length(y)
lin<-ncol(X)
n<-length(k)
q<-kat-1
Eta<-cbind(X,W)%*%Delta_start
Mu<-as.vector(family$linkinv(Eta))

Sigma<-matrix(0,N,N)
Sig1<-list()
  for (rrr in 1:(N/q))
  {
  Sig1<-as.matrix(Mu[(((rrr-1)*q)+1):(rrr*q)]%*%t(1-Mu[(((rrr-1)*q)+1):(rrr*q)]))
  Sig1[lower.tri(Sig1)]<-t(Sig1)[lower.tri(Sig1)]
  Sigma[(((rrr-1)*q)+1):(rrr*q),(((rrr-1)*q)+1):(rrr*q)]<-Sig1
}
D<-as.vector(family$mu.eta(Eta))

Z_alles<-cbind(X,W)

if(s==1)
{
P1<-c(rep(0,lin),rep((q_start^(-1)),n*s))
P1<-diag(P1)
}else{
P1<-matrix(0,lin+n*s,lin+n*s)
for(jf in 1:n)
P1[(lin+(jf-1)*s+1):(lin+jf*s),(lin+(jf-1)*s+1):(lin+jf*s)]<-chol2inv(chol(q_start))
}

Delta<-matrix(0,steps,(lin+s*n))
Delta[1,]<-Delta_start

Q<-list()
Q[[1]]<-q_start

score_vec<-rep(0,(lin+s*n))
D<-as.vector(family$mu.eta(Eta))
Delta_r<-rep(0,lin+s*n)

l=1

if(print.iter.final)
  print(paste("Final Re-estimation Iteration ", 1,sep=""))

opt<-steps

Sigmainv<-chol2inv(chol(Sigma))
score_part<-t(Z_alles)%*%diag(D)%*%Sigmainv
score_vec<-score_part%*%(y-Mu)-P1%*%Delta[1,]
F_gross<-(score_part%*%t(diag(D))%*%Z_alles)+P1

InvFisher<-chol2inv(chol(F_gross))
Delta_r<-InvFisher%*%score_vec
Halbierungsindex<-0

index<-seq(from=1,to=(kat-1)*(length(katvar)+1),by=length(katvar)+1)
lin.test<-(kat-1)*(length(katvar)+1)
for (ij in 0:100)
{
Test_delta<-Delta[1,1:lin.test]+(0.5)^(ij)*Delta_r[1:lin.test]
if(is.null(katvar))
{
Test_delta_diff<-diff(diag(length(index)),differences=1)%*%Test_delta
}else{
Test_matrix<-matrix(0,dim(X)[1]/q,kat-1)

for (jt in 1:(dim(X)[1]/q))
Test_matrix[jt,]<-X[c((jt-1)*(kat-1)+1):(jt*(kat-1)),1:lin.test]%*%Test_delta
Test_delta_diff<-diff(diag(length(index)),differences=1)%*%t(Test_matrix)
}
Halbierungsindex<-ij

if(all (Test_delta_diff>0))
         break
}

Delta[1,]<-Delta[1,]+((0.5)^Halbierungsindex)*Delta_r

Eta<-Z_alles%*%Delta[1,]

Mu<-as.vector(family$linkinv(Eta))

Sigma<-matrix(0,N,N)
Sig1<-list()
  for (rrr in 1:(N/q))
  {
  Sig1<-as.matrix(Mu[(((rrr-1)*q)+1):(rrr*q)]%*%t(1-Mu[(((rrr-1)*q)+1):(rrr*q)]))
  Sig1[lower.tri(Sig1)]<-t(Sig1)[lower.tri(Sig1)]
  Sigma[(((rrr-1)*q)+1):(rrr*q),(((rrr-1)*q)+1):(rrr*q)]<-Sig1
}
D<-as.vector(family$mu.eta(Eta))

if(method=="EM")
{
F_gross<-t(Z_alles)%*%diag(D)%*%chol2inv(chol(Sigma))%*%t(diag(D))%*%Z_alles+P1

InvFisher<-chol2inv(chol(F_gross))
############################# Q updaten ################
Q1<-InvFisher[(lin+1):(lin+s),(lin+1):(lin+s)]+Delta[1,(lin+1):(lin+s)]%*%t(Delta[1,(lin+1):(lin+s)])

for (i in 2:n)
Q1<-Q1+InvFisher[(lin+(i-1)*s+1):(lin+i*s),(lin+(i-1)*s+1):(lin+i*s)]+Delta[1,(lin+(i-1)*s+1):(lin+i*s)]%*%t(Delta[1,(lin+(i-1)*s+1):(lin+i*s)])

Q1<-1/n*Q1
}else{
Eta_tilde<-Eta+(y-Mu)*1/D

Betadach<-Delta[1,1:lin]

if(s==1)
{
optim.obj<-nlminb(sqrt(q_start),likelihood_bobyqa_cumul,D=D,Sigma=Sigma,X=X,X_aktuell=X,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W,lower = 1e-14, upper=20)
Q1<-as.matrix(optim.obj$par)^2
}else{
q_start_vec<-c(diag(q_start),q_start[lower.tri(q_start)])
up1<-min(20,50*max(q_start_vec))#[(s+1):(s*(s+1)*0.5)]))
upp<-rep(up1,length(q_start_vec))
low<-c(rep(0,s),rep(-up1,0.5*(s^2-s)))
kkk_vec<-c(rep(-1,s),rep(0.5,0.5*(s^2-s)))
optim.obj<-bobyqa(q_start_vec,likelihood_cumul,D=D,Sigma=Sigma,X=X,X_aktuell=X,Eta_tilde=Eta_tilde,Betadach=Betadach,W=W,n=n,s=s,k=k,lower=low,upper=upp)
Q1<-matrix(0,s,s)
Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
Q1<-Q1+t(Q1)
diag(Q1)<-(optim.obj$par[1:s])

#### Check for positive definitness ########
      for (ttt in 0:100)
      {
      Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
      Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
      Q_solvetest<-try(solve(Q1))
         if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
         break
      }
}}

Q[[2]]<-Q1
###############################################################################################################################################
################################################################### Boost ###################################################################
eps<-eps.final*sqrt(length(Delta_r))

for (l in 2:steps)
{
  
  if(print.iter.final)
    print(paste("Final Re-estimation Iteration ", l,sep=""))
  
  if(s==1)
  {
  P1<-c(rep(0,lin),rep((Q1^(-1)),n*s))
  P1<-diag(P1)
  }else{
  P1<-matrix(0,lin+n*s,lin+n*s)
  for(jf in 1:n)
  P1[(lin+(jf-1)*s+1):(lin+jf*s),(lin+(jf-1)*s+1):(lin+jf*s)]<-chol2inv(chol(Q1))
  }

  Sigmainv<-chol2inv(chol(Sigma))
  score_part<-t(Z_alles)%*%diag(D)%*%Sigmainv
  score_vec<-score_part%*%(y-Mu)-P1%*%Delta[l-1,]
  F_gross<-(score_part%*%t(diag(D))%*%Z_alles)+P1


InvFisher<-chol2inv(chol(F_gross))
Delta_r<-InvFisher%*%score_vec
Halbierungsindex<-0


########## Check order restriction
for (ij in 0:1000)
{
Test_delta<-Delta[l-1,1:lin.test]+(0.5)^(ij)*Delta_r[1:lin.test]
if(is.null(katvar))
{
Test_delta_diff<-diff(diag(length(index)),differences=1)%*%Test_delta
}else{
Test_matrix<-matrix(0,dim(X)[1]/q,kat-1)

for (jt in 1:(dim(X)[1]/q))
Test_matrix[jt,]<-X[c((jt-1)*(kat-1)+1):(jt*(kat-1)),1:lin.test]%*%Test_delta
Test_delta_diff<-diff(diag(length(index)),differences=1)%*%t(Test_matrix)
}
Halbierungsindex<-ij

if(all (Test_delta_diff>0))
         break
}

Delta[l,]<-Delta[l-1,]+((0.5)^Halbierungsindex)*Delta_r

Eta<-Z_alles%*%Delta[l,]

Mu<-as.vector(family$linkinv(Eta))

Sigma<-matrix(0,N,N)
Sig1<-list()
  for (rrr in 1:(N/q))
  {
  Sig1<-as.matrix(Mu[(((rrr-1)*q)+1):(rrr*q)]%*%t(1-Mu[(((rrr-1)*q)+1):(rrr*q)]))
  Sig1[lower.tri(Sig1)]<-t(Sig1)[lower.tri(Sig1)]
  Sigma[(((rrr-1)*q)+1):(rrr*q),(((rrr-1)*q)+1):(rrr*q)]<-Sig1
}
  
Solve_test<-try(chol2inv(chol(Sigma)))

kek<-0
while(class(Solve_test)=="try-error")
{
  Delta[l,]<-Delta[l-1,]+(0.9^(kek+1))*Delta_r

Eta<-Z_alles%*%Delta[l,]

Mu<-as.vector(family$linkinv(Eta))

Sigma<-matrix(0,N,N)
Sig1<-list()
  for (rrr in 1:(N/q))
  {
  Sig1<-as.matrix(Mu[(((rrr-1)*q)+1):(rrr*q)]%*%t(1-Mu[(((rrr-1)*q)+1):(rrr*q)]))
  Sig1[lower.tri(Sig1)]<-t(Sig1)[lower.tri(Sig1)]
  Sigma[(((rrr-1)*q)+1):(rrr*q),(((rrr-1)*q)+1):(rrr*q)]<-Sig1
}
Solve_test<-try(chol2inv(chol(Sigma)))
kek<-kek+1
}

D<-as.vector(family$mu.eta(Eta))

if (method=="EM")
{

F_gross<-t(Z_alles)%*%diag(D)%*%chol2inv(chol(Sigma))%*%t(diag(D))%*%Z_alles+P1

InvFisher<-chol2inv(chol(F_gross))
############################# Q update ################
Q1<-InvFisher[(lin++1):(lin++s),(lin++1):(lin++s)]+Delta[l,(lin++1):(lin++s)]%*%t(Delta[l,(lin++1):(lin++s)])
for (i in 2:n)
Q1<-Q1+InvFisher[(lin++(i-1)*s+1):(lin++i*s),(lin++(i-1)*s+1):(lin++i*s)]+Delta[l,(lin++(i-1)*s+1):(lin++i*s)]%*%t(Delta[l,(lin++(i-1)*s+1):(lin++i*s)])

Q1<-1/n*Q1
}else{
Eta_tilde<-Eta+(y-Mu)*1/D

Betadach<-Delta[l,1:lin]

if(s==1)
{
optim.obj<-nlminb(sqrt(Q1),likelihood_bobyqa_cumul,D=D,Sigma=Sigma,X=X,X_aktuell=X,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = 1e-12, upper = 20)
Q1<-as.matrix(optim.obj$par)^2
}else{
Q1_vec<-c(diag(Q1),Q1[lower.tri(Q1)])
optim.obj<-bobyqa(Q1_vec,likelihood_cumul,D=D,Sigma=Sigma,X=X,X_aktuell=X,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W,lower=low,upper=upp)

Q1<-matrix(0,s,s)
Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
Q1<-Q1+t(Q1)
diag(Q1)<-(optim.obj$par[1:s])

#### Check for positiv definitness ########
for (ttt in 0:100)
      {
      Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
      Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
       Q_solvetest<-try(solve(Q1))
         if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
         break
      }
}}

Q[[l+1]]<-Q1

kritval<-sqrt(sum((Delta[l-1,]-Delta[l,])^2))/sqrt(sum(Delta[l-1,]^2))
if(kritval<eps)
break

if(l>2)
{
kritval2<-sqrt(sum((Delta[l-2,]-Delta[l,])^2))/sqrt(sum(Delta[l-2,]^2))
if(kritval2<eps)
break
}}

opt<-l
Deltafinal<-Delta[l,]
Q_final<-Q[[l+1]]

  if(s==1)
  {
  P1<-c(rep(0,lin),rep((Q_final^(-1)),n*s))
  P1<-diag(P1)
  }else{
  P1<-matrix(0,lin+n*s,lin+n*s)
  for(jf in 1:n)
  P1[(lin+(jf-1)*s+1):(lin+jf*s),(lin+(jf-1)*s+1):(lin+jf*s)]<-chol2inv(chol(Q_final))
  }

F_gross<-t(Z_alles)%*%diag(D)%*%chol2inv(chol(Sigma))%*%t(diag(D))%*%Z_alles+P1

Inv_F_opt<-chol2inv(chol(F_gross))

Standard_errors<-sqrt(diag(Inv_F_opt))

ret.obj=list()
ret.obj$opt<-opt
ret.obj$Delta<-Deltafinal
ret.obj$Q<-Q_final
ret.obj$Standard_errors<-Standard_errors
return(ret.obj)
}


###############################################

glmm_final_seq<-function(y,X,W,k,q_start,Delta_start,s,steps=1000,
                         method,kat,katvar,print.iter.final=FALSE,eps.final=1e-5)
{
  family=binomial(link="logit")
  N<-length(y)
  lin<-ncol(X)
  n<-length(k)
  q<-kat-1
  Eta<-cbind(X,W)%*%Delta_start
  Mu<-as.vector(family$linkinv(Eta))
  
  Sigma<-as.vector(family$variance(Mu))
  
  D<-as.vector(family$mu.eta(Eta))
  
  Z_alles<-cbind(X,W)
  
  if(s==1)
  {
    P1<-c(rep(0,lin),rep((q_start^(-1)),n*s))
    P1<-diag(P1)
  }else{
    P1<-matrix(0,lin+n*s,lin+n*s)
    for(jf in 1:n)
      P1[(lin+(jf-1)*s+1):(lin+jf*s),(lin+(jf-1)*s+1):(lin+jf*s)]<-chol2inv(chol(q_start))
  }
  
  Delta<-matrix(0,steps,(lin+s*n))
  Delta[1,]<-Delta_start
  
  Q<-list()
  Q[[1]]<-q_start
  
  score_vec<-rep(0,(lin+s*n))
  D<-as.vector(family$mu.eta(Eta))
  Delta_r<-rep(0,lin+s*n)
  
  l=1
  opt<-steps
  if(print.iter.final)
    print(paste("Final Re-estimation Iteration ", 1,sep=""))
  

  score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[1,]
  F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1

  
  InvFisher<-chol2inv(chol(F_gross))
  Delta_r<-InvFisher%*%score_vec
  Halbierungsindex<-0
    
  Delta[1,]<-Delta[1,]+((0.5)^Halbierungsindex)*Delta_r
  
  Eta<-Z_alles%*%Delta[1,]
  
  Mu<-as.vector(family$linkinv(Eta))
  
  Sigma<-as.vector(family$variance(Mu))

  
  D<-as.vector(family$mu.eta(Eta))
  
  if(method=="EM")
  {
    
    F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
  
    InvFisher<-chol2inv(chol(F_gross))
    ############################# Q updaten ################
    Q1<-InvFisher[(lin+1):(lin+s),(lin+1):(lin+s)]+Delta[1,(lin+1):(lin+s)]%*%t(Delta[1,(lin+1):(lin+s)])
    
    for (i in 2:n)
      Q1<-Q1+InvFisher[(lin+(i-1)*s+1):(lin+i*s),(lin+(i-1)*s+1):(lin+i*s)]+Delta[1,(lin+(i-1)*s+1):(lin+i*s)]%*%t(Delta[1,(lin+(i-1)*s+1):(lin+i*s)])
    
    Q1<-1/n*Q1
  }else{
    Eta_tilde<-Eta+(y-Mu)*1/D
    
    Betadach<-Delta[1,1:lin]
    
    if(s==1)
    {
      optim.obj<-nlminb(sqrt(q_start),likelihood_bobyqa_seq,D=D,Sigma=Sigma,X=X,X_aktuell=X,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W,lower = 1e-14, upper=20)
      Q1<-as.matrix(optim.obj$par)^2
    }else{
      q_start_vec<-c(diag(q_start),q_start[lower.tri(q_start)])
      up1<-min(20,50*max(q_start_vec))#[(s+1):(s*(s+1)*0.5)]))
      upp<-rep(up1,length(q_start_vec))
      low<-c(rep(0,s),rep(-up1,0.5*(s^2-s)))
      kkk_vec<-c(rep(-1,s),rep(0.5,0.5*(s^2-s)))
      optim.obj<-bobyqa(q_start_vec,likelihood_seq,D=D,Sigma=Sigma,X=X,X_aktuell=X,Eta_tilde=Eta_tilde,Betadach=Betadach,W=W,n=n,s=s,k=k,lower=low,upper=upp)
      Q1<-matrix(0,s,s)
      Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
      Q1<-Q1+t(Q1)
      diag(Q1)<-(optim.obj$par[1:s])
      
      #### Check for positive definitness ########
      for (ttt in 0:100)
      {
        Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
        Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
        Q_solvetest<-try(solve(Q1))
        if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
          break
      }
    }}
  
  Q[[2]]<-Q1
  ###############################################################################################################################################
  ################################################################### Boost ###################################################################
  eps<-eps.final*sqrt(length(Delta_r))
  
  for (l in 2:steps)
  {
    
    if(print.iter.final)
      print(paste("Final Re-estimation Iteration ", l,sep=""))
    
    if(s==1)
    {
      P1<-c(rep(0,lin),rep((Q1^(-1)),n*s))
      P1<-diag(P1)
    }else{
      P1<-matrix(0,lin+n*s,lin+n*s)
      for(jf in 1:n)
        P1[(lin+(jf-1)*s+1):(lin+jf*s),(lin+(jf-1)*s+1):(lin+jf*s)]<-chol2inv(chol(Q1))
    }
    
      score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[l-1,]
      F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
    
    InvFisher<-chol2inv(chol(F_gross))
    Delta_r<-InvFisher%*%score_vec
    Halbierungsindex<-0
    
    
    Delta[l,]<-Delta[l-1,]+((0.5)^Halbierungsindex)*Delta_r
    
    Eta<-Z_alles%*%Delta[l,]
    
    Mu<-as.vector(family$linkinv(Eta))
    
    Sigma<-as.vector(family$variance(Mu))
    
            
    D<-as.vector(family$mu.eta(Eta))
    
    if (method=="EM")
    {
      
      F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
      
      InvFisher<-chol2inv(chol(F_gross))
      ############################# Q update ################
      Q1<-InvFisher[(lin++1):(lin++s),(lin++1):(lin++s)]+Delta[l,(lin++1):(lin++s)]%*%t(Delta[l,(lin++1):(lin++s)])
      for (i in 2:n)
        Q1<-Q1+InvFisher[(lin++(i-1)*s+1):(lin++i*s),(lin++(i-1)*s+1):(lin++i*s)]+Delta[l,(lin++(i-1)*s+1):(lin++i*s)]%*%t(Delta[l,(lin++(i-1)*s+1):(lin++i*s)])
      
      Q1<-1/n*Q1
    }else{
      Eta_tilde<-Eta+(y-Mu)*1/D
      
      Betadach<-Delta[l,1:lin]
      
      if(s==1)
      {
        optim.obj<-nlminb(sqrt(Q1),likelihood_bobyqa_seq,D=D,Sigma=Sigma,X=X,X_aktuell=X,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W,lower = 1e-12, upper = 20)
        Q1<-as.matrix(optim.obj$par)^2
      }else{
        Q1_vec<-c(diag(Q1),Q1[lower.tri(Q1)])
        optim.obj<-bobyqa(Q1_vec,likelihood_seq,D=D,Sigma=Sigma,X=X,X_aktuell=X,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W,lower=low,upper=upp)
        
        Q1<-matrix(0,s,s)
        Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
        Q1<-Q1+t(Q1)
        diag(Q1)<-(optim.obj$par[1:s])
        
        #### Check for positiv definitness ########
        for (ttt in 0:100)
        {
          Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
          Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
          Q_solvetest<-try(solve(Q1))
          if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
            break
        }
      }}
    
    Q[[l+1]]<-Q1
    
    kritval<-sqrt(sum((Delta[l-1,]-Delta[l,])^2))/sqrt(sum(Delta[l-1,]^2))
    if(kritval<eps)
      break
    
    if(l>2)
    {
      kritval2<-sqrt(sum((Delta[l-2,]-Delta[l,])^2))/sqrt(sum(Delta[l-2,]^2))
      if(kritval2<eps)
        break
    }}
  
  opt<-l
  Deltafinal<-Delta[l,]
  Q_final<-Q[[l+1]]
  
  if(s==1)
  {
    P1<-c(rep(0,lin),rep((Q_final^(-1)),n*s))
    P1<-diag(P1)
  }else{
    P1<-matrix(0,lin+n*s,lin+n*s)
    for(jf in 1:n)
      P1[(lin+(jf-1)*s+1):(lin+jf*s),(lin+(jf-1)*s+1):(lin+jf*s)]<-chol2inv(chol(Q_final))
  }
  
  F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
  
  Inv_F_opt<-chol2inv(chol(F_gross))
  
  Standard_errors<-sqrt(diag(Inv_F_opt))
  
  ret.obj=list()
  ret.obj$opt<-opt
  ret.obj$Delta<-Deltafinal
  ret.obj$Q<-Q_final
  ret.obj$Standard_errors<-Standard_errors
  return(ret.obj)
}

#############################################################################################################################################
###################################################### main boosting function #################################################################

est.OrdinalBoost<-function(fix,rnd,data,model="sequential",control=list())
{
  if(model=="cumulative")
  {
  ret.obj<-est.OrdinalBoost.cumul(fix=fix,rnd=rnd,data=data,control=control)
  }else{
  ret.obj<-est.OrdinalBoost.seq(fix=fix,rnd=rnd,data=data,control=control)
  }
  ret.obj$model<-model
  return(ret.obj)
}  

#############################################################################################################################################
###################################################### main boosting function #################################################################
#fix=formula(y_katvec~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10);rnd=list(Person2=~1);data=H_kat;model="cumulative";control=list(steps=10,katvar=c("X3","X5"))
#fix=formula(y_katvec~V1+V2+V3+V4+V5+V6+V7+V8+as.factor(X9)+as.factor(X10));rnd=list(Person2=~1);data=H_kat;model="sequential";control=list(steps=10,lin=c("V6","as.factor(X10)"))
#fix=formula(pain~time+I(time^2)+th+age+sex);rnd=list(id=~1+time);data=knee; model="sequential";control=list(steps=10,katvar="age",lin="time")
est.OrdinalBoost.cumul<-function(fix,rnd,data,control=list())
{
family=binomial(link="logit")
variables<-attr(terms(fix),"variables")

if(attr(terms(fix), "intercept")==0)
{
variables<-attr(terms(fix),"term.labels")
fix<- paste(rownames((attr(terms(fix),"factors")))[1],"~ +1",sep="")
for (ir in 1:length(variables))
fix <- paste(fix, variables[ir], sep="+")
fix <-formula(fix)
}

y <- model.response(model.frame(fix, data))

X <- model.matrix(fix, data)
X<-X[,-1]

very.old.names<-attr(terms(fix),"term.labels")
very.old.names2<-attr(terms(fix),"term.labels")

old.names<-attr(X,"dimnames")[[2]]
old.names2<-attr(X,"dimnames")[[2]]

factor.names<-character()
for (i in 1:length(very.old.names))
{
if (substr(very.old.names[i],1,9)=="as.factor")
factor.names<-c(factor.names,very.old.names[i])
}

factor.list<-list()
if(length(factor.names)>0)
{
spl<-strsplit(factor.names,"\\(")

categ.names<-character()
for(uz in 1:length(spl))
categ.names<-c(categ.names,spl[[uz]][2])

spl2<-strsplit(categ.names,"\\)")
categ.names2<-character()
for(uz in 1:length(spl2))
categ.names2<-c(categ.names2,spl2[[uz]])
factor.names<-categ.names2

for (i in 1:length(factor.names))
factor.list[[i]]<-levels(as.factor(data[,factor.names[i]]))
}    





rndformula <- as.character(rnd)

trmsrnd <- terms(rnd[[1]])
newrndfrml <- "~ -1"
newrndfrml <- paste(newrndfrml,  if(attr(trmsrnd, "intercept")) names(rnd)[1] else "", sep=" + ")

if(length(attr(trmsrnd, "variables"))>1)
{
newrndfrml <- paste(newrndfrml,
         paste(sapply(attr(trmsrnd,"term.labels"), function(lbl){
                     paste(lbl, names(rnd)[1], sep=":")
                 }), collapse=" + "), sep="+") }

W_start <- model.matrix(formula(newrndfrml), data)

rnlabels<-terms(formula(newrndfrml))
random.labels<-attr(rnlabels,"term.labels")
k<-table(data[,colnames(data)==(names(rnd)[1])])
n<-length(k)
s<-dim(W_start)[2]/n

if(s>1)
{
W<-W_start[,seq(from=1,to=1+(s-1)*n,by=n)]
for (i in 2:n)
W<-cbind(W,W_start[,seq(from=i,to=i+(s-1)*n,by=n)])
}else{
W<-W_start
}

print(paste("Iteration ", 1,sep=""))
control<-do.call(OrdinalBoostControl, control)

katvar<-control$katvar
method<-control$method

old.control.lin<-control$lin

if(sum(substr(control$lin,1,9)=="as.factor")>0)
{
group0<-substr(control$lin,1,9)=="as.factor"
spl0<-strsplit(control$lin[group0],"\\)")

control$lin<-control$lin[!group0]
for(ur in 1:length(spl0))
control$lin<-c(control$lin,old.names[is.element(substr(old.names,1,nchar(spl0[[ur]])),spl0[[ur]])])
}

lin.out<-NULL
lin.out2<-NULL

if(!is.null(control$lin) || !is.null(control$katvar))
{
lin.out<-!is.element(very.old.names,c(control$katvar,old.control.lin[1:length(old.control.lin)]))
very.old.names<-very.old.names[lin.out]
lin.out2<-!is.element(old.names,c(control$katvar,control$lin[1:length(control$lin)]))
old.names<-old.names[lin.out2]
}

group<-substr(very.old.names,1,9)=="as.factor"

if(sum(group)>0)
{
spl<-strsplit(very.old.names[group],"\\(")

categ.names<-character()
for(uz in 1:length(spl))
categ.names<-c(categ.names,spl[[uz]][2])

spl2<-strsplit(categ.names,"\\)")
categ.names2<-character()
for(uz in 1:length(spl2))
categ.names2<-c(categ.names2,spl2[[uz]])

block<-numeric()
posi<-1
for(ip in 1:length(group))
{
  if(!group[ip])
  {
  block<-c(block,1)
  }else{
  block<-c(block,length(levels(as.factor(data[,categ.names2[posi]])))-1)
  posi<-posi+1
  }
}}else{
block<-rep(1,length(group))
}

SQ<-kumulize(y,cbind(X,W),katvar)


y<-SQ$y_bin
X<-cbind(SQ$X_new[,1:((SQ$kat-1)*(1+length(control$katvar)))],SQ$X_new[,control$lin])
colnames(X)<-c(colnames(SQ$X_new[,1:((SQ$kat-1)*(1+length(control$katvar)))]),control$lin)
U<-SQ$X_new[,old.names]
W<-SQ$X_new[,colnames(W)]

kat<-SQ$kat

if(length(control$q_start)==0)
{
control$q_start<-rep(0.1,s)
if(s>1)
control$q_start<-diag(control$q_start)
}

index<-seq(from=1,to=(kat-1)*(length(katvar)+1),by=length(katvar)+1)
if(is.null(control$start))
{
control$start<-rep(0,(dim(X)[2]+dim(W)[2]))

control$start[index]<-c(1:length(index))
control$start[index]<-control$start[index]-mean(control$start[index])
}

beta_null<-control$start[1:dim(X)[2]]
ranef_null<-control$start[(dim(X)[2]+1):(dim(X)[2]+dim(W)[2])]
q_start<-control$q_start

N<-length(y)

Z_fastalles<-cbind(X,U)

lin<-length(beta_null)

m<-length(block)
m_alt<-ncol(U)

lin0<-sum(beta_null!=0)
if(s==1)
{
Q_start<-diag(q_start^2,s)
p_start<-c(rep(0,lin0),rep((q_start^2)^(-1),n*s))
P_start<-diag(p_start)
}else{
Q_start<-q_start
P_start<-matrix(0,lin0+n*s,lin0+n*s)
for(jf in 1:n)
P_start[(lin0+(jf-1)*s+1):(lin0+jf*s),(lin0+(jf-1)*s+1):(lin0+jf*s)]<-chol2inv(chol(q_start))
}

Eta_start<-X%*%beta_null+W%*%ranef_null
D_start<-as.vector(family$mu.eta(Eta_start))
Mu_start<-as.vector(family$linkinv(Eta_start))

q<-kat-1

Sigma_start<-matrix(0,N,N)
  for (rrr in 1:(N/q))
  {
  Sig1<-as.matrix(Mu_start[(((rrr-1)*q)+1):(rrr*q)]%*%t(1-Mu_start[(((rrr-1)*q)+1):(rrr*q)]))
  Sig1[lower.tri(Sig1)]<-t(Sig1)[lower.tri(Sig1)]
  Sigma_start[(((rrr-1)*q)+1):(rrr*q),(((rrr-1)*q)+1):(rrr*q)]<-Sig1
}
Z_start<-cbind(X[,beta_null!=0],W)

W0_inv<-diag(D_start)%*%chol2inv(chol(Sigma_start))%*%diag(D_start)
M0<-Z_start%*%chol2inv(chol((t(Z_start)%*%W0_inv%*%Z_start)+P_start))%*%t(Z_start)%*%W0_inv

blocksum<-rep(0,m)
for (jg in 1:m)
blocksum[jg]<-sum(block[1:jg])

Mu<-as.vector(family$linkinv(Eta_start))

Z_alles<-cbind(X,U,W)

#########################  some definitions ################
Z_r<-list()

Z_r[[1]]<-cbind(X,as.matrix(U[,1:block[1]]))
for (r in 2:m)
Z_r[[r]]<-cbind(X,as.matrix(U[,(blocksum[r-1]+1):blocksum[r]]))

Z_r_gross<-list()
for (r in 1:m)
Z_r_gross[[r]]<-cbind(Z_r[[r]],W)

Sigma<-matrix(0,N,N)
Sig1<-list()
  for (rrr in 1:(N/q))
  {
  Sig1<-as.matrix(Mu[(((rrr-1)*q)+1):(rrr*q)]%*%t(1-Mu[(((rrr-1)*q)+1):(rrr*q)]))
  Sig1[lower.tri(Sig1)]<-t(Sig1)[lower.tri(Sig1)]
  Sigma[(((rrr-1)*q)+1):(rrr*q),(((rrr-1)*q)+1):(rrr*q)]<-Sig1
}
########################################################## some definitions ################################################
Delta<-matrix(0,control$steps,(lin+m_alt+s*n))

Standard<-matrix(0,control$steps,(lin+m_alt+s*n))

Standard_ma<-list()

for (ijk in 1:m)
Standard_ma[[ijk]]<-as.numeric()

Delta[1,1:lin]<-beta_null
Delta[1,(lin+m_alt+1):(lin+m_alt+s*n)]<-t(ranef_null)

Delta_start<-Delta[1,]

IC<-matrix(0,nrow=control$steps,ncol=m)
IC_stop<-rep(0,control$steps)
komp<-rep(0,control$steps)

E<-diag(N)

M1<-list()
diag_H1<-rep(0,m)
FinalHat<-list()
Halbierungsindex<-rep(0,m)

Q<-array(0,c(s,s,control$steps+1))
Q[,,1]<-Q_start
######### derive D and score_r and F_r  ###########
D<-as.vector(family$mu.eta(Eta_start))

Delta_r<-list()
for(r in 1:m)
Delta_r[[r]]<-rep(0,lin+block[r]+s*n)

Eta_r<-matrix(0,N,m)
Mu_r<-matrix(0,N,m)

grp<-rep(1,m_alt)
for (r in 2:m_alt)
grp[((r-1)+1):(r)]<-r

if(s==1)
{
P_alles<-c(rep(0,lin+m_alt),rep((Q_start^(-1)),n*s))
P_alles<-diag(P_alles)
}else{
P_alles<-matrix(0,lin+m_alt+n*s,lin+m_alt+n*s)
for(jf in 1:n)
P_alles[(lin+m_alt+(jf-1)*s+1):(lin+m_alt+jf*s),(lin+m_alt+(jf-1)*s+1):(lin+m_alt+jf*s)]<-chol2inv(chol(Q_start))
}

Sigmainv<-chol2inv(chol(Sigma))
score_part<-t(Z_alles)%*%diag(D)%*%Sigmainv
score_alles<-score_part%*%(y-Mu)
F_alles<-(score_part%*%t(diag(D))%*%Z_alles)+P_alles

Inv_r<-list()

l=1
for (r in 1:m)
{
  if(r==1)
  {
  Inv_r[[1]]<-chol2inv(chol(F_alles[c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s)),c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s))]))
  Delta_r[[r]]<-Inv_r[[1]]%*%score_alles[c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s))]
  }else{
  Inv_r[[r]]<-chol2inv(chol(F_alles[c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s)),c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s))]))
  Delta_r[[r]]<-Inv_r[[r]]%*%score_alles[c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s))]
  }

########## Check order restriction
lin.test<-(kat-1)*(length(katvar)+1)
for (ij in 0:100)
{
Test_delta<-Delta[1,1:lin.test]+(0.5)^(ij)*Delta_r[[r]][1:lin.test]
if(is.null(katvar))
{
Test_delta_diff<-diff(diag(length(index)),differences=1)%*%Test_delta
}else{
Test_matrix<-matrix(0,dim(X)[1]/q,kat-1)

for (jt in 1:(dim(X)[1]/q))
Test_matrix[jt,]<-X[c((jt-1)*(kat-1)+1):(jt*(kat-1)),1:lin.test]%*%Test_delta
Test_delta_diff<-diff(diag(length(index)),differences=1)%*%t(Test_matrix)
}
Halbierungsindex[r]<-ij

if(all (Test_delta_diff>0))
         break
}

Eta_r[,r]<-Eta_start+((0.5)^Halbierungsindex[r])*Z_r_gross[[r]]%*%Delta_r[[r]]
############# new  Mu_r ##############
Mu_r[,r]<-as.vector(family$linkinv(Eta_r[,r]))

################### derive Hat-Matrix  ##############
Eigen1<-eigen(Sigma)

P.sigma<-Eigen1$vectors
L.sigma<-sqrt(diag(Eigen1$values))
Sigma.halbe<-P.sigma%*%L.sigma%*%t(P.sigma)

Sigmainv<-chol2inv(chol(Sigma))
Eigen2<-eigen(t(t(Sigmainv*D)*D))

P.w<-Eigen2$vectors
L.w<-sqrt(diag(Eigen2$values))
w.halbe<-P.w%*%L.w%*%t(P.w)

M1[[r]]<-Sigma.halbe%*%w.halbe%*%(Z_r_gross[[r]])%*%Inv_r[[r]]%*%t(Z_r_gross[[r]])%*%w.halbe%*%chol2inv(chol(Sigma.halbe))

diag_H1[r]<-sum(diag(M0+M1[[r]]%*%(E-M0)))

######## likelihood of information criterion   ########
    IC[1,r]<-sum(y*log(Mu_r[,r])+(1-y)*log(1-Mu_r[,r]))#-0.5*pena

   ############ AIC ##############
    if (control$sel.method=="aic")
    IC[1,r]<-(-IC[1,r])+diag_H1[r]

   ############ BIC ###############
    if (control$sel.method=="bic")
    IC[1,r]<-(-IC[1,r])+0.5*diag_H1[r]*log(n)
}
   mi_ic<-match(min(IC[1,]),IC[1,])

#################### update ###########################
if(mi_ic==1)
{
g<-c(rep(T,lin),is.element(grp,1:blocksum[1]),rep(T,s*n))
}else{
g<-c(rep(T,lin),is.element(grp,(blocksum[mi_ic-1]+1):blocksum[mi_ic]),rep(T,s*n))
}

Delta_r[[mi_ic]]<-((0.5)^Halbierungsindex[mi_ic])*Delta_r[[mi_ic]]
Delta_r[[mi_ic]][(lin+1):(lin+block[mi_ic])]<-control$nue*Delta_r[[mi_ic]][(lin+1):(lin+block[mi_ic])]
Delta[1,g]<-Delta[1,g]+t(Delta_r[[mi_ic]])

Nue_ma<-((0.5)^Halbierungsindex[mi_ic])*c(rep(1,lin),rep(control$nue,block[mi_ic]),rep(1,s*n))

Uu<-(E-(Sigma.halbe%*%w.halbe%*%(Z_r_gross[[mi_ic]])%*%(Inv_r[[mi_ic]]*Nue_ma)%*%t(Z_r_gross[[mi_ic]])%*%w.halbe%*%chol2inv(chol(Sigma.halbe))))%*%(E-M0)

FinalHat[[1]]<-E-Uu

defre<-sum(diag(FinalHat[[1]]))

komp[1]<-mi_ic

Eta<-Z_alles%*%Delta[1,]

Mu<-as.vector(family$linkinv(Eta))

Sigma<-matrix(0,N,N)
Sig1<-list()
  for (rrr in 1:(N/q))
  {
  Sig1<-as.matrix(Mu[(((rrr-1)*q)+1):(rrr*q)]%*%t(1-Mu[(((rrr-1)*q)+1):(rrr*q)]))
  Sig1[lower.tri(Sig1)]<-t(Sig1)[lower.tri(Sig1)]
  Sigma[(((rrr-1)*q)+1):(rrr*q),(((rrr-1)*q)+1):(rrr*q)]<-Sig1
}

D<-as.vector(family$mu.eta(Eta))

if (control$method=="EM")
{
aktuell_vec<-!is.element(Delta[1,],0)
X_aktuell<-Z_alles[,aktuell_vec]

m_aktuell<-sum(block[is.element(1:m,komp)])

if(s==1)
{
P_gross<-c(rep(0,lin+m_aktuell),rep((Q_start^(-1)),n*s))
P_gross<-diag(P_gross)
}else{
P_gross<-matrix(0,lin+m_aktuell+n*s,lin+m_aktuell+n*s)
for(jf in 1:n)
P_gross[(lin+m_aktuell+(jf-1)*s+1):(lin+m_aktuell+jf*s),(lin+m_aktuell+(jf-1)*s+1):(lin+m_aktuell+jf*s)]<-chol2inv(chol(Q_start))
}

F_gross<-t(X_aktuell)%*%diag(D)%*%chol2inv(chol(Sigma))%*%t(diag(D))%*%X_aktuell+P_gross

############################# update Q ################
InvFisher<-chol2inv(chol(F_gross))

Q1<-InvFisher[(lin+m_aktuell+1):(lin+m_aktuell+s),(lin+m_aktuell+1):(lin+m_aktuell+s)]+Delta[1,(lin+m+1):(lin+m+s)]%*%t(Delta[1,(lin+m+1):(lin+m+s)])

for (i in 2:n)
Q1<-Q1+InvFisher[(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s),(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s)]+Delta[1,(lin+m+(i-1)*s+1):(lin+m+i*s)]%*%t(Delta[1,(lin+m+(i-1)*s+1):(lin+m+i*s)])

Q1<-1/n*Q1
}else{
Eta_tilde<-Eta+(y-Mu)*1/D

Betadach<-Delta[1,1:(lin+m_alt)]

aktuell_vec<-!is.element(Delta[1,1:(lin+m_alt)],0)
X_aktuell<-Z_fastalles[,aktuell_vec]

if(s==1)
{
optim.obj<-nlminb(sqrt(Q_start),likelihood_bobyqa_cumul,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W,lower = 1e-14, upper=20)
Q1<-as.matrix(optim.obj$par)^2
}else{
q_start_vec<-c(diag(q_start),q_start[lower.tri(q_start)])
up1<-min(20,50*max(q_start_vec))#[(s+1):(s*(s+1)*0.5)]))
upp<-rep(up1,length(q_start_vec))
low<-c(rep(0,s),rep(-up1,0.5*(s^2-s)))
kkk_vec<-c(rep(-1,s),rep(0.5,0.5*(s^2-s)))
optim.obj<-bobyqa(q_start_vec,likelihood_cumul,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,Betadach=Betadach,W=W,n=n,s=s,k=k,lower=low,upper=upp)
Q1<-matrix(0,s,s)
Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
Q1<-Q1+t(Q1)
diag(Q1)<-(optim.obj$par[1:s])

#### Check for positive definitness ########
      for (ttt in 0:100)
      {
      Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
      Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
      Q_solvetest<-try(solve(Q1))
         if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
         break
      }
}
}
    IC_stop[1]<-sum(y*log(Mu)+(1-y)*log(1-Mu))
   ############ AIC ##############
    if (control$sel.method=="aic")
    IC_stop[1]<-(-IC_stop[1])+defre

   ############ BIC ###############
    if (control$sel.method=="bic")
    IC_stop[1]<-(-IC_stop[1])+0.5*defre*log(n)

Q[,,2]<-Q1
###############################################################################################################################################
################################################################### Boost ###################################################################
if(control$steps!=1)
{
for (l in 2:control$steps)
{
print(paste("Iteration ", l,sep=""))

if(s==1)
{
P_alles<-c(rep(0,lin+m_alt),rep((Q1^(-1)),n*s))
P_alles<-diag(P_alles)
}else{
P_alles<-matrix(0,lin+m_alt+n*s,lin+m_alt+n*s)
for(jf in 1:n)
P_alles[(lin+m_alt+(jf-1)*s+1):(lin+m_alt+jf*s),(lin+m_alt+(jf-1)*s+1):(lin+m_alt+jf*s)]<-chol2inv(chol(Q1))
}

Sigmainv<-chol2inv(chol(Sigma))
score_part<-t(Z_alles)%*%diag(D)%*%Sigmainv
score_alles<-score_part%*%(y-Mu)
F_alles<-(score_part%*%t(diag(D))%*%Z_alles)+P_alles

  for (r in 1:m)
  {
  if(r==1)
  {
  Inv_r[[1]]<-chol2inv(chol(F_alles[c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s)),c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s))]))
  Delta_r[[r]]<-Inv_r[[1]]%*%score_alles[c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s))]
  }else{
  Inv_r[[r]]<-chol2inv(chol(F_alles[c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s)),c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s))]))
  Delta_r[[r]]<-Inv_r[[r]]%*%score_alles[c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s))]
  }

   ########## Check order restriction
for (ij in 0:100)
{
Test_delta<-Delta[l-1,1:lin.test]+(0.5)^(ij)*Delta_r[[r]][1:lin.test]
if(is.null(katvar))
{
Test_delta_diff<-diff(diag(length(index)),differences=1)%*%Test_delta
}else{
Test_matrix<-matrix(0,dim(X)[1]/q,kat-1)

for (jt in 1:(dim(X)[1]/q))
Test_matrix[jt,]<-X[c((jt-1)*(kat-1)+1):(jt*(kat-1)),1:lin.test]%*%Test_delta
Test_delta_diff<-diff(diag(length(index)),differences=1)%*%t(Test_matrix)
}
Halbierungsindex[r]<-ij

if(all (Test_delta_diff>0))
         break
}

   Eta_r[,r]<-Eta+((0.5)^Halbierungsindex[r])*Z_r_gross[[r]]%*%Delta_r[[r]]

   Mu_r[,r]<-as.vector(family$linkinv(Eta_r[,r]))

     ################### Jetzt Hat-Matrix berechnen ##############
      Eigen1<-eigen(Sigma)

      P.sigma<-Eigen1$vectors
      L.sigma<-sqrt(diag(Eigen1$values))
      Sigma.halbe<-P.sigma%*%L.sigma%*%t(P.sigma)

      Sigmainv<-chol2inv(chol(Sigma))
      Eigen2<-eigen(t(t(Sigmainv*D)*D))

      P.w<-Eigen2$vectors
      L.w<-sqrt(diag(Eigen2$values))
      w.halbe<-P.w%*%L.w%*%t(P.w)

      M1[[r]]<-Sigma.halbe%*%w.halbe%*%(Z_r_gross[[r]])%*%Inv_r[[r]]%*%t(Z_r_gross[[r]])%*%w.halbe%*%chol2inv(chol(Sigma.halbe))

  diag_H1[r]<-sum(diag(E-((E-M1[[r]])%*%Uu)))

      IC[l,r]<-sum(y*log(Mu_r[,r])+(1-y)*log(1-Mu_r[,r]))#-0.5*pena

   ############ AIC ##############
      if (control$sel.method=="aic")
      IC[l,r]<-(-IC[l,r])+diag_H1[r]

   ############ BIC ###############
      if (control$sel.method=="bic")
      IC[l,r]<-(-IC[l,r])+0.5*diag_H1[r]*log(n)
   }

   mi_ic<-match(min(IC[l,]),IC[l,])

   if(mi_ic==1)
   {
   g<-c(rep(T,lin),is.element(grp,1:blocksum[1]),rep(T,s*n))
   }else{
   g<-c(rep(T,lin),is.element(grp,(blocksum[mi_ic-1]+1):blocksum[mi_ic]),rep(T,s*n))
   }

   Delta[l,]<-Delta[l-1,]
   Delta_r[[mi_ic]]<-((0.5)^Halbierungsindex[mi_ic])*Delta_r[[mi_ic]]
   Delta_r[[mi_ic]][(lin+1):(lin+block[mi_ic])]<-control$nue*Delta_r[[mi_ic]][(lin+1):(lin+block[mi_ic])]
   Delta[l,g]<-Delta[l,g]+t(Delta_r[[mi_ic]])

komp[l]<-mi_ic

Nue_ma<-((0.5)^Halbierungsindex[mi_ic])*c(rep(1,lin),rep(control$nue,block[mi_ic]),rep(1,s*n))

Uu<-(E-(Sigma.halbe%*%w.halbe%*%(Z_r_gross[[mi_ic]])%*%(Inv_r[[mi_ic]]*Nue_ma)%*%t(Z_r_gross[[mi_ic]])%*%w.halbe%*%chol2inv(chol(Sigma.halbe))))%*%Uu

FinalHat[[l]]<-E-Uu
defre<-sum(diag(FinalHat[[l]]))

Eta<-Z_alles%*%Delta[l,]

Mu<-as.vector(family$linkinv(Eta))

Sigma<-matrix(0,N,N)
Sig1<-list()
  for (rrr in 1:(N/q))
  {
  Sig1<-as.matrix(Mu[(((rrr-1)*q)+1):(rrr*q)]%*%t(1-Mu[(((rrr-1)*q)+1):(rrr*q)]))
  Sig1[lower.tri(Sig1)]<-t(Sig1)[lower.tri(Sig1)]
  Sigma[(((rrr-1)*q)+1):(rrr*q),(((rrr-1)*q)+1):(rrr*q)]<-Sig1
}
D<-as.vector(family$mu.eta(Eta))

if (control$method=="EM")
{
aktuell_vec<-!is.element(Delta[l,],0)
X_aktuell<-Z_alles[,aktuell_vec]

m_aktuell<-sum(block[is.element(1:m,komp)])


if(s==1)
{
P_gross<-c(rep(0,lin+m_aktuell),rep((Q1^(-1)),n*s))
P_gross<-diag(P_gross)
}else{
P_gross<-matrix(0,lin+m_aktuell+n*s,lin+m_aktuell+n*s)
for(jf in 1:n)
P_gross[(lin+m_aktuell+(jf-1)*s+1):(lin+m_aktuell+jf*s),(lin+m_aktuell+(jf-1)*s+1):(lin+m_aktuell+jf*s)]<-chol2inv(chol(Q1))
}


F_gross<-t(X_aktuell)%*%diag(D)%*%chol2inv(chol(Sigma))%*%t(diag(D))%*%X_aktuell+P_gross

############################# update Q ################
InvFisher<-chol2inv(chol(F_gross))

Q1<-InvFisher[(lin+m_aktuell+1):(lin+m_aktuell+s),(lin+m_aktuell+1):(lin+m_aktuell+s)]+Delta[l,(lin+m+1):(lin+m+s)]%*%t(Delta[l,(lin+m+1):(lin+m+s)])

for (i in 2:n)
Q1<-Q1+InvFisher[(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s),(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s)]+Delta[l,(lin+m+(i-1)*s+1):(lin+m+i*s)]%*%t(Delta[l,(lin+m+(i-1)*s+1):(lin+m+i*s)])


Q1<-1/n*Q1
}else{
Eta_tilde<-Eta+(y-Mu)*1/D

Betadach<-Delta[l,1:(lin+m_alt)]

aktuell_vec<-!is.element(Delta[l,1:(lin+m_alt)],0)
X_aktuell<-Z_fastalles[,aktuell_vec]

if(s==1)
{
optim.obj<-try(nlminb(sqrt(Q1),likelihood_bobyqa_cumul,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W,lower = 1e-12, upper = 20))
Q1<-as.matrix(optim.obj$par)^2
}else{
Q1_vec<-c(diag(Q1),Q1[lower.tri(Q1)])
optim.obj<-try(bobyqa(Q1_vec,likelihood_cumul,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W,lower=low,upper=upp))

Q1<-matrix(0,s,s)
Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
Q1<-Q1+t(Q1)
diag(Q1)<-(optim.obj$par[1:s])


#### Check for positiv definitness ########
for (ttt in 0:100)
      {
      Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
      Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
       Q_solvetest<-try(solve(Q1))
         if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
         break
      }
}
}

    IC_stop[l]<-sum(y*log(Mu)+(1-y)*log(1-Mu))

   ############ AIC ##############
    if (control$sel.method=="aic")
    IC_stop[l]<-(-IC_stop[l])+defre

   ############ BIC ###############
    if (control$sel.method=="bic")
    IC_stop[l]<-(-IC_stop[l])+0.5*defre*log(n)

Q[,,l+1]<-Q1

kritval<-sqrt(sum((Delta[l-1,]-Delta[l,])^2))/sqrt(sum(Delta[l,]^2))


if(kritval<1e-8)
break


kritval2<-abs(IC_stop[l-1]-IC_stop[l])/abs(IC_stop[l])

if(kritval2<1e-8)
break


}}
if (l<control$steps)
{
for (wj in (l+1):control$steps)
{
Delta[wj,]<-Delta[l,]
Q[,,wj+1]<-Q1
IC[wj,]<-IC[l,]
IC_stop[wj]<-IC_stop[l]
}}

##########
opt<-match(min(IC_stop),IC_stop)

if (control$OPT==TRUE)
{
Delta_neu<-Delta[opt,]
}else{
Delta_neu<-Delta[l,]
}

Qfinal<-Q[,,opt+1]

aaa<-!is.element(Delta_neu[1:(lin+m_alt)],0)

glmm_final<-try(glmm_final_cumul(y,Z_fastalles[,aaa],W,k,q_start=Qfinal,
                                 Delta_start=Delta_neu[c(aaa,rep(T,n*s))],s,
                                 steps=control$maxIter,method=method,kat=kat,
                                 katvar=katvar,print.iter.final=control$print.iter.final,
                                 eps.final=control$eps.final))

if(class(glmm_final)=="try-error" || glmm_final$opt>(control$maxIter-5))
{
glmm_final<-try(glmm_final_cumul(y,Z_fastalles[,aaa],W,k,q_start=q_start,
                                 Delta_start=Delta_start[c(aaa,rep(T,n*s))],s,
                                 steps=2000,method=method,kat=kat,katvar=katvar,
                                 print.iter.final=control$print.iter.final,eps.final=control$eps.final))

if(class(glmm_final)=="try-error" || glmm_final$opt>1990)
{
cat("Warning:\n")
cat("Final Fisher scoring reestimation did not converge!")
}}
Delta_neu[c(aaa,rep(T,n*s))]<-glmm_final$Delta

Standard_errors<-rep(0,length(Delta_neu))
Standard_errors[c(aaa,rep(T,n*s))]<-glmm_final$Standard_errors

Qfinal<-glmm_final$Q

Eta_opt<-Z_alles%*%Delta_neu
Mu_opt<-as.vector(family$linkinv(Eta_opt))

if(s==1)
Qfinal<-sqrt(Qfinal)


if(!is.matrix(Qfinal))
Qfinal<-as.matrix(Qfinal)
colnames(Qfinal)<-random.labels
rownames(Qfinal)<-random.labels

names(Delta_neu)[1:lin]<-colnames(X)[1:lin]
names(Standard_errors)[1:lin]<-colnames(X)[1:lin]


names(Delta_neu)[(lin+1):(lin+m_alt)]<-colnames(U)
names(Standard_errors)[(lin+1):(lin+m_alt)]<-colnames(U)


permut<-numeric()
for(ys in 1:length(old.names2[!is.element(old.names2,control$katvar)]))
permut<-c(permut,match(old.names2[!is.element(old.names2,control$katvar)][ys],names(Delta_neu)))

Delta_neu[((SQ$kat-1)*(1+length(control$katvar))+1):((SQ$kat-1)*(1+length(control$katvar))+length(old.names2[!is.element(old.names2,control$katvar)]))]<-Delta_neu[permut]
Standard_errors[((SQ$kat-1)*(1+length(control$katvar))+1):((SQ$kat-1)*(1+length(control$katvar))+length(old.names2[!is.element(old.names2,control$katvar)]))]<-Standard_errors[permut]

names(Delta_neu)[((SQ$kat-1)*(1+length(control$katvar))+1):((SQ$kat-1)*(1+length(control$katvar))+length(old.names2[!is.element(old.names2,control$katvar)]))]<-old.names2[!is.element(old.names2,control$katvar)]
names(Standard_errors)[((SQ$kat-1)*(1+length(control$katvar))+1):((SQ$kat-1)*(1+length(control$katvar))+length(old.names2[!is.element(old.names2,control$katvar)]))]<-old.names2[!is.element(old.names2,control$katvar)]


Delta[,((SQ$kat-1)*(1+length(control$katvar))+1):((SQ$kat-1)*(1+length(control$katvar))+length(old.names2[!is.element(old.names2,control$katvar)]))]<-Delta[,permut]

names(Delta_neu)[(lin+m_alt+1):(lin+m_alt+n*s)]<-colnames(W)
names(Standard_errors)[(lin+m_alt+1):(lin+m_alt+n*s)]<-colnames(W)
colnames(Delta)<-names(Delta_neu)

colnames(IC)<-very.old.names

comp<-character()
for(oi in 1:length(komp))
comp<-c(comp,very.old.names[komp[oi]])

ret.obj=list()
ret.obj$IC<-IC
ret.obj$IC_sel<-IC_stop
ret.obj$opt<-opt
ret.obj$Deltamatrix<-Delta
ret.obj$ranef<-Delta_neu[(lin+m_alt+1):(lin+m_alt+n*s)]
ret.obj$coefficients<-Delta_neu[1:(lin+m_alt)]
ret.obj$fixerror<-Standard_errors[1:(lin+m_alt)]
ret.obj$ranerror<-Standard_errors[(lin+m_alt+1):(lin+m_alt+n*s)]
ret.obj$components<-comp
ret.obj$Q_long<-Q
ret.obj$Q<-Qfinal
ret.obj$y_hat<-Mu_opt
ret.obj$HatMatrix<-FinalHat[[opt]]
ret.obj$family<-family
ret.obj$fix<-fix
ret.obj$newrndfrml<-newrndfrml
ret.obj$subject<-names(rnd)[1]
ret.obj$k<-k
ret.obj$kat<-kat
ret.obj$katvar<-katvar
ret.obj$data<-data
ret.obj$factor.names<-factor.names
ret.obj$factor.list<-factor.list
return(ret.obj)
}




#############################################################################################################################################
###################################################### main boosting function #################################################################
#fix=formula(y_katvec~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10);rnd=list(Person2=~1);data=H_kat;model="cumulative";control=list(steps=10,katvar=c("X3","X5"))
#fix=formula(y_katvec~V1+V2+V3+V4+V5+V6+V7+V8+as.factor(X9)+as.factor(X10));rnd=list(Person2=~1);data=H_kat;model="sequential";control=list(steps=10,lin=c("V6","as.factor(X10)"))
#fix=formula(pain~time+I(time^2)+th+age+sex);rnd=list(id=~1+time);data=knee; model="sequential";control=list(steps=10,katvar="age",lin="time")
est.OrdinalBoost.seq<-function(fix,rnd,data,control=list())
{
  family=binomial(link="logit")
  variables<-attr(terms(fix),"variables")
  
  if(attr(terms(fix), "intercept")==0)
  {
    variables<-attr(terms(fix),"term.labels")
    fix<- paste(rownames((attr(terms(fix),"factors")))[1],"~ +1",sep="")
    for (ir in 1:length(variables))
      fix <- paste(fix, variables[ir], sep="+")
    fix <-formula(fix)
  }
  
  y <- model.response(model.frame(fix, data))
  
  X <- model.matrix(fix, data)
  X<-X[,-1]
  
  very.old.names<-attr(terms(fix),"term.labels")
  very.old.names2<-attr(terms(fix),"term.labels")
  
  old.names<-attr(X,"dimnames")[[2]]
  old.names2<-attr(X,"dimnames")[[2]]
  
  factor.names<-character()
  for (i in 1:length(very.old.names))
  {
    if (substr(very.old.names[i],1,9)=="as.factor")
      factor.names<-c(factor.names,very.old.names[i])
  }
  
  factor.list<-list()
  if(length(factor.names)>0)
  {
    spl<-strsplit(factor.names,"\\(")
    
    categ.names<-character()
    for(uz in 1:length(spl))
      categ.names<-c(categ.names,spl[[uz]][2])
    
    spl2<-strsplit(categ.names,"\\)")
    categ.names2<-character()
    for(uz in 1:length(spl2))
      categ.names2<-c(categ.names2,spl2[[uz]])
    factor.names<-categ.names2
    
    for (i in 1:length(factor.names))
      factor.list[[i]]<-levels(as.factor(data[,factor.names[i]]))
  }    
  
  
  
  
  
  rndformula <- as.character(rnd)
  
  trmsrnd <- terms(rnd[[1]])
  newrndfrml <- "~ -1"
  newrndfrml <- paste(newrndfrml,  if(attr(trmsrnd, "intercept")) names(rnd)[1] else "", sep=" + ")
  
  if(length(attr(trmsrnd, "variables"))>1)
  {
    newrndfrml <- paste(newrndfrml,
                        paste(sapply(attr(trmsrnd,"term.labels"), function(lbl){
                          paste(lbl, names(rnd)[1], sep=":")
                        }), collapse=" + "), sep="+") }
  
  W_start <- model.matrix(formula(newrndfrml), data)
  
  rnlabels<-terms(formula(newrndfrml))
  random.labels<-attr(rnlabels,"term.labels")
  k<-table(data[,colnames(data)==(names(rnd)[1])])
  n<-length(k)
  s<-dim(W_start)[2]/n
  
  if(s>1)
  {
    W<-W_start[,seq(from=1,to=1+(s-1)*n,by=n)]
    for (i in 2:n)
      W<-cbind(W,W_start[,seq(from=i,to=i+(s-1)*n,by=n)])
  }else{
    W<-W_start
  }
  
  print(paste("Iteration ", 1,sep=""))
  control<-do.call(OrdinalBoostControl, control)
  
  katvar<-control$katvar
  method<-control$method
  
  old.control.lin<-control$lin
  
  if(sum(substr(control$lin,1,9)=="as.factor")>0)
  {
    group0<-substr(control$lin,1,9)=="as.factor"
    spl0<-strsplit(control$lin[group0],"\\)")
    
    control$lin<-control$lin[!group0]
    for(ur in 1:length(spl0))
      control$lin<-c(control$lin,old.names[is.element(substr(old.names,1,nchar(spl0[[ur]])),spl0[[ur]])])
  }
  
  lin.out<-NULL
  lin.out2<-NULL
  
  if(!is.null(control$lin) || !is.null(control$katvar))
  {
    lin.out<-!is.element(very.old.names,c(control$katvar,old.control.lin[1:length(old.control.lin)]))
    very.old.names<-very.old.names[lin.out]
    lin.out2<-!is.element(old.names,c(control$katvar,control$lin[1:length(control$lin)]))
    old.names<-old.names[lin.out2]
  }
  
  group<-substr(very.old.names,1,9)=="as.factor"
  
  if(sum(group)>0)
  {
    spl<-strsplit(very.old.names[group],"\\(")
    
    categ.names<-character()
    for(uz in 1:length(spl))
      categ.names<-c(categ.names,spl[[uz]][2])
    
    spl2<-strsplit(categ.names,"\\)")
    categ.names2<-character()
    for(uz in 1:length(spl2))
      categ.names2<-c(categ.names2,spl2[[uz]])
    
    block<-numeric()
    posi<-1
    for(ip in 1:length(group))
    {
      if(!group[ip])
      {
        block<-c(block,1)
      }else{
        block<-c(block,length(levels(as.factor(data[,categ.names2[posi]])))-1)
        posi<-posi+1
      }
    }}else{
      block<-rep(1,length(group))
    }
  
    SQ<-sequentialize(y,cbind(X,W),katvar)
  
  y<-SQ$y_bin
  X<-cbind(SQ$X_new[,1:((SQ$kat-1)*(1+length(control$katvar)))],SQ$X_new[,control$lin])
  colnames(X)<-c(colnames(SQ$X_new[,1:((SQ$kat-1)*(1+length(control$katvar)))]),control$lin)
  U<-SQ$X_new[,old.names]
  W<-SQ$X_new[,colnames(W)]
  
  kat<-SQ$kat
  
  if(length(control$q_start)==0)
  {
    control$q_start<-rep(0.1,s)
    if(s>1)
      control$q_start<-diag(control$q_start)
  }
  
  index<-seq(from=1,to=(kat-1)*(length(katvar)+1),by=length(katvar)+1)
  if(is.null(control$start))
  control$start<-rep(0,(dim(X)[2]+dim(W)[2]))
   
  beta_null<-control$start[1:dim(X)[2]]
  ranef_null<-control$start[(dim(X)[2]+1):(dim(X)[2]+dim(W)[2])]
  q_start<-control$q_start
  
  N<-length(y)
  
  Z_fastalles<-cbind(X,U)
  
  lin<-length(beta_null)
  
  m<-length(block)
  m_alt<-ncol(U)
  
  lin0<-sum(beta_null!=0)
  if(s==1)
  {
    Q_start<-diag(q_start^2,s)
    p_start<-c(rep(0,lin0),rep((q_start^2)^(-1),n*s))
    P_start<-diag(p_start)
  }else{
    Q_start<-q_start
    P_start<-matrix(0,lin0+n*s,lin0+n*s)
    for(jf in 1:n)
      P_start[(lin0+(jf-1)*s+1):(lin0+jf*s),(lin0+(jf-1)*s+1):(lin0+jf*s)]<-chol2inv(chol(q_start))
  }
  
  Eta_start<-X%*%beta_null+W%*%ranef_null
  D_start<-as.vector(family$mu.eta(Eta_start))
  Mu_start<-as.vector(family$linkinv(Eta_start))
  
  q<-kat-1
  
      Sigma_start<-as.vector(family$variance(Mu_start))
  
  Z_start<-cbind(X[,beta_null!=0],W)
  
    W0_inv<-D_start*1/Sigma_start*D_start
    M0<-Z_start%*%chol2inv(chol((t(Z_start)%*%(Z_start*W0_inv))+P_start))%*%t(Z_start*W0_inv)
  
  blocksum<-rep(0,m)
  for (jg in 1:m)
    blocksum[jg]<-sum(block[1:jg])
  
  Mu<-as.vector(family$linkinv(Eta_start))
  
  Z_alles<-cbind(X,U,W)
  
  #########################  some definitions ################
  Z_r<-list()
  
  Z_r[[1]]<-cbind(X,as.matrix(U[,1:block[1]]))
  for (r in 2:m)
    Z_r[[r]]<-cbind(X,as.matrix(U[,(blocksum[r-1]+1):blocksum[r]]))
  
  Z_r_gross<-list()
  for (r in 1:m)
    Z_r_gross[[r]]<-cbind(Z_r[[r]],W)
  
      Sigma<-as.vector(family$variance(Mu))
  
  ########################################################## some definitions ################################################
  Delta<-matrix(0,control$steps,(lin+m_alt+s*n))
  
  Standard<-matrix(0,control$steps,(lin+m_alt+s*n))
  
  Standard_ma<-list()
  
  for (ijk in 1:m)
    Standard_ma[[ijk]]<-as.numeric()
  
  Delta[1,1:lin]<-beta_null
  Delta[1,(lin+m_alt+1):(lin+m_alt+s*n)]<-t(ranef_null)
  
  Delta_start<-Delta[1,]
  
  IC<-matrix(0,nrow=control$steps,ncol=m)
  IC_stop<-rep(0,control$steps)
  komp<-rep(0,control$steps)
  
  E<-diag(N)
  
  M1<-list()
  diag_H1<-rep(0,m)
  FinalHat<-list()
  Halbierungsindex<-rep(0,m)
  
  Q<-array(0,c(s,s,control$steps+1))
  Q[,,1]<-Q_start
  ######### derive D and score_r and F_r  ###########
  D<-as.vector(family$mu.eta(Eta_start))
  
  Delta_r<-list()
  for(r in 1:m)
    Delta_r[[r]]<-rep(0,lin+block[r]+s*n)
  
  Eta_r<-matrix(0,N,m)
  Mu_r<-matrix(0,N,m)
  
  grp<-rep(1,m_alt)
  for (r in 2:m_alt)
    grp[((r-1)+1):(r)]<-r
  
  if(s==1)
  {
    P_alles<-c(rep(0,lin+m_alt),rep((Q_start^(-1)),n*s))
    P_alles<-diag(P_alles)
  }else{
    P_alles<-matrix(0,lin+m_alt+n*s,lin+m_alt+n*s)
    for(jf in 1:n)
      P_alles[(lin+m_alt+(jf-1)*s+1):(lin+m_alt+jf*s),(lin+m_alt+(jf-1)*s+1):(lin+m_alt+jf*s)]<-chol2inv(chol(Q_start))
  }
  
    score_alles<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
    F_alles<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P_alles
  
  Inv_r<-list()
  
  l=1
  for (r in 1:m)
  {
    if(r==1)
    {
      Inv_r[[1]]<-chol2inv(chol(F_alles[c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s)),c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s))]))
      Delta_r[[r]]<-Inv_r[[1]]%*%score_alles[c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s))]
    }else{
      Inv_r[[r]]<-chol2inv(chol(F_alles[c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s)),c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s))]))
      Delta_r[[r]]<-Inv_r[[r]]%*%score_alles[c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s))]
    }
    
    
    Eta_r[,r]<-Eta_start+((0.5)^Halbierungsindex[r])*Z_r_gross[[r]]%*%Delta_r[[r]]
    ############# new  Mu_r ##############
    Mu_r[,r]<-as.vector(family$linkinv(Eta_r[,r]))
    
    ################### derive Hat-Matrix  ##############
    M1[[r]]<-(Z_r_gross[[r]]*sqrt(Sigma*D*1/Sigma*D))%*%Inv_r[[r]]%*%t(Z_r_gross[[r]]*sqrt(D*1/Sigma*D*1/Sigma))
 
    diag_H1[r]<-sum(diag(M0+M1[[r]]%*%(E-M0)))
    
    ######## likelihood of information criterion   ########
    IC[1,r]<-sum(y*log(Mu_r[,r])+(1-y)*log(1-Mu_r[,r]))#-0.5*pena
    
    ############ AIC ##############
    if (control$sel.method=="aic")
      IC[1,r]<-(-IC[1,r])+diag_H1[r]
    
    ############ BIC ###############
    if (control$sel.method=="bic")
      IC[1,r]<-(-IC[1,r])+0.5*diag_H1[r]*log(n)
  }
  mi_ic<-match(min(IC[1,]),IC[1,])
  
  #################### update ###########################
  if(mi_ic==1)
  {
    g<-c(rep(T,lin),is.element(grp,1:blocksum[1]),rep(T,s*n))
  }else{
    g<-c(rep(T,lin),is.element(grp,(blocksum[mi_ic-1]+1):blocksum[mi_ic]),rep(T,s*n))
  }
  
  Delta_r[[mi_ic]]<-((0.5)^Halbierungsindex[mi_ic])*Delta_r[[mi_ic]]
  Delta_r[[mi_ic]][(lin+1):(lin+block[mi_ic])]<-control$nue*Delta_r[[mi_ic]][(lin+1):(lin+block[mi_ic])]
  Delta[1,g]<-Delta[1,g]+t(Delta_r[[mi_ic]])
  
  Nue_ma<-((0.5)^Halbierungsindex[mi_ic])*c(rep(1,lin),rep(control$nue,block[mi_ic]),rep(1,s*n))
  
    Uu<-(E-(cbind(Z_r[[mi_ic]],W)*sqrt(Sigma*D*1/Sigma*D))%*%(Inv_r[[mi_ic]]*Nue_ma)%*%t(cbind(Z_r[[mi_ic]],W)*sqrt(D*1/Sigma*D*1/Sigma)))%*%(E-M0)

  FinalHat[[1]]<-E-Uu
  
  defre<-sum(diag(FinalHat[[1]]))
  
  komp[1]<-mi_ic
  
  Eta<-Z_alles%*%Delta[1,]
  
  Mu<-as.vector(family$linkinv(Eta))
  
      Sigma<-as.vector(family$variance(Mu))
  
  D<-as.vector(family$mu.eta(Eta))
  
  if (control$method=="EM")
  {
    aktuell_vec<-!is.element(Delta[1,],0)
    X_aktuell<-Z_alles[,aktuell_vec]
    
    m_aktuell<-sum(block[is.element(1:m,komp)])
    
    if(s==1)
    {
      P_gross<-c(rep(0,lin+m_aktuell),rep((Q_start^(-1)),n*s))
      P_gross<-diag(P_gross)
    }else{
      P_gross<-matrix(0,lin+m_aktuell+n*s,lin+m_aktuell+n*s)
      for(jf in 1:n)
        P_gross[(lin+m_aktuell+(jf-1)*s+1):(lin+m_aktuell+jf*s),(lin+m_aktuell+(jf-1)*s+1):(lin+m_aktuell+jf*s)]<-chol2inv(chol(Q_start))
    }
    
      F_gross<-t(X_aktuell)%*%(X_aktuell*D*1/Sigma*D)+P_gross
    
    ############################# update Q ################
    InvFisher<-chol2inv(chol(F_gross))
    
    Q1<-InvFisher[(lin+m_aktuell+1):(lin+m_aktuell+s),(lin+m_aktuell+1):(lin+m_aktuell+s)]+Delta[1,(lin+m+1):(lin+m+s)]%*%t(Delta[1,(lin+m+1):(lin+m+s)])
    
    for (i in 2:n)
      Q1<-Q1+InvFisher[(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s),(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s)]+Delta[1,(lin+m+(i-1)*s+1):(lin+m+i*s)]%*%t(Delta[1,(lin+m+(i-1)*s+1):(lin+m+i*s)])
    
    Q1<-1/n*Q1
  }else{
    Eta_tilde<-Eta+(y-Mu)*1/D
    
    Betadach<-Delta[1,1:(lin+m_alt)]
    
    aktuell_vec<-!is.element(Delta[1,1:(lin+m_alt)],0)
    X_aktuell<-Z_fastalles[,aktuell_vec]
    
    if(s==1)
    {
      optim.obj<-nlminb(sqrt(Q_start),likelihood_bobyqa_seq,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W,lower = 1e-14, upper=20)
      Q1<-as.matrix(optim.obj$par)^2
    }else{
      q_start_vec<-c(diag(q_start),q_start[lower.tri(q_start)])
      up1<-min(20,50*max(q_start_vec))#[(s+1):(s*(s+1)*0.5)]))
      upp<-rep(up1,length(q_start_vec))
      low<-c(rep(0,s),rep(-up1,0.5*(s^2-s)))
      kkk_vec<-c(rep(-1,s),rep(0.5,0.5*(s^2-s)))
      optim.obj<-bobyqa(q_start_vec,likelihood_seq,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,Betadach=Betadach,W=W,n=n,s=s,k=k,lower=low,upper=upp)
      Q1<-matrix(0,s,s)
      Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
      Q1<-Q1+t(Q1)
      diag(Q1)<-(optim.obj$par[1:s])
      
      #### Check for positive definitness ########
      for (ttt in 0:100)
      {
        Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
        Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
        Q_solvetest<-try(solve(Q1))
        if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
          break
      }
    }
  }
  IC_stop[1]<-sum(y*log(Mu)+(1-y)*log(1-Mu))
  ############ AIC ##############
  if (control$sel.method=="aic")
    IC_stop[1]<-(-IC_stop[1])+defre
  
  ############ BIC ###############
  if (control$sel.method=="bic")
    IC_stop[1]<-(-IC_stop[1])+0.5*defre*log(n)
  
  Q[,,2]<-Q1
  ###############################################################################################################################################
  ################################################################### Boost ###################################################################
  if(control$steps!=1)
  {
    for (l in 2:control$steps)
    {
      print(paste("Iteration ", l,sep=""))
      
      if(s==1)
      {
        P_alles<-c(rep(0,lin+m_alt),rep((Q1^(-1)),n*s))
        P_alles<-diag(P_alles)
      }else{
        P_alles<-matrix(0,lin+m_alt+n*s,lin+m_alt+n*s)
        for(jf in 1:n)
          P_alles[(lin+m_alt+(jf-1)*s+1):(lin+m_alt+jf*s),(lin+m_alt+(jf-1)*s+1):(lin+m_alt+jf*s)]<-chol2inv(chol(Q1))
      }
      
        score_alles<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
        F_alles<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P_alles
      
      for (r in 1:m)
      {
        if(r==1)
        {
          Inv_r[[1]]<-chol2inv(chol(F_alles[c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s)),c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s))]))
          Delta_r[[r]]<-Inv_r[[1]]%*%score_alles[c(1:(lin+block[1]),(lin+m_alt+1):(lin+m_alt+n*s))]
        }else{
          Inv_r[[r]]<-chol2inv(chol(F_alles[c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s)),c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s))]))
          Delta_r[[r]]<-Inv_r[[r]]%*%score_alles[c(1:lin,(lin+sum(block[1:(r-1)])+1):(lin+sum(block[1:r])),(lin+m_alt+1):(lin+m_alt+n*s))]
        }
        
        ########## Check order restriction
        
        Eta_r[,r]<-Eta+((0.5)^Halbierungsindex[r])*Z_r_gross[[r]]%*%Delta_r[[r]]
        
        Mu_r[,r]<-as.vector(family$linkinv(Eta_r[,r]))
        
        ################### Jetzt Hat-Matrix berechnen ##############
          M1[[r]]<-(Z_r_gross[[r]]*sqrt(Sigma*D*1/Sigma*D))%*%Inv_r[[r]]%*%t(Z_r_gross[[r]]*sqrt(D*1/Sigma*D*1/Sigma))

        diag_H1[r]<-sum(diag(E-((E-M1[[r]])%*%Uu)))
        
        IC[l,r]<-sum(y*log(Mu_r[,r])+(1-y)*log(1-Mu_r[,r]))#-0.5*pena
        
        ############ AIC ##############
        if (control$sel.method=="aic")
          IC[l,r]<-(-IC[l,r])+diag_H1[r]
        
        ############ BIC ###############
        if (control$sel.method=="bic")
          IC[l,r]<-(-IC[l,r])+0.5*diag_H1[r]*log(n)
      }
      
      mi_ic<-match(min(IC[l,]),IC[l,])
      
      if(mi_ic==1)
      {
        g<-c(rep(T,lin),is.element(grp,1:blocksum[1]),rep(T,s*n))
      }else{
        g<-c(rep(T,lin),is.element(grp,(blocksum[mi_ic-1]+1):blocksum[mi_ic]),rep(T,s*n))
      }
      
      Delta[l,]<-Delta[l-1,]
      Delta_r[[mi_ic]]<-((0.5)^Halbierungsindex[mi_ic])*Delta_r[[mi_ic]]
      Delta_r[[mi_ic]][(lin+1):(lin+block[mi_ic])]<-control$nue*Delta_r[[mi_ic]][(lin+1):(lin+block[mi_ic])]
      Delta[l,g]<-Delta[l,g]+t(Delta_r[[mi_ic]])
      
      komp[l]<-mi_ic
      
      Nue_ma<-((0.5)^Halbierungsindex[mi_ic])*c(rep(1,lin),rep(control$nue,block[mi_ic]),rep(1,s*n))
      
        Uu<-(E-(cbind(Z_r[[mi_ic]],W)*sqrt(Sigma*D*1/Sigma*D))%*%(Inv_r[[mi_ic]]*Nue_ma)%*%t(cbind(Z_r[[mi_ic]],W)*sqrt(D*1/Sigma*D*1/Sigma)))%*%Uu
      
      FinalHat[[l]]<-E-Uu
      defre<-sum(diag(FinalHat[[l]]))
      
      Eta<-Z_alles%*%Delta[l,]
      
      Mu<-as.vector(family$linkinv(Eta))
      
          Sigma<-as.vector(family$variance(Mu))
      
      D<-as.vector(family$mu.eta(Eta))
      
      if (control$method=="EM")
      {
        aktuell_vec<-!is.element(Delta[l,],0)
        X_aktuell<-Z_alles[,aktuell_vec]
        
        m_aktuell<-sum(block[is.element(1:m,komp)])
        
        
        if(s==1)
        {
          P_gross<-c(rep(0,lin+m_aktuell),rep((Q1^(-1)),n*s))
          P_gross<-diag(P_gross)
        }else{
          P_gross<-matrix(0,lin+m_aktuell+n*s,lin+m_aktuell+n*s)
          for(jf in 1:n)
            P_gross[(lin+m_aktuell+(jf-1)*s+1):(lin+m_aktuell+jf*s),(lin+m_aktuell+(jf-1)*s+1):(lin+m_aktuell+jf*s)]<-chol2inv(chol(Q1))
        }
        
          F_gross<-t(X_aktuell)%*%(X_aktuell*D*1/Sigma*D)+P_gross
        
        ############################# update Q ################
        InvFisher<-chol2inv(chol(F_gross))
        
        Q1<-InvFisher[(lin+m_aktuell+1):(lin+m_aktuell+s),(lin+m_aktuell+1):(lin+m_aktuell+s)]+Delta[l,(lin+m+1):(lin+m+s)]%*%t(Delta[l,(lin+m+1):(lin+m+s)])
        
        for (i in 2:n)
          Q1<-Q1+InvFisher[(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s),(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s)]+Delta[l,(lin+m+(i-1)*s+1):(lin+m+i*s)]%*%t(Delta[l,(lin+m+(i-1)*s+1):(lin+m+i*s)])
        
        
        Q1<-1/n*Q1
      }else{
        Eta_tilde<-Eta+(y-Mu)*1/D
        
        Betadach<-Delta[l,1:(lin+m_alt)]
        
        aktuell_vec<-!is.element(Delta[l,1:(lin+m_alt)],0)
        X_aktuell<-Z_fastalles[,aktuell_vec]
        
        if(s==1)
        {
          optim.obj<-try(nlminb(sqrt(Q1),likelihood_bobyqa_seq,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W,lower = 1e-12, upper = 20))
          Q1<-as.matrix(optim.obj$par)^2
        }else{
          Q1_vec<-c(diag(Q1),Q1[lower.tri(Q1)])
          optim.obj<-try(bobyqa(Q1_vec,likelihood_seq,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W,lower=low,upper=upp))
          
          Q1<-matrix(0,s,s)
          Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
          Q1<-Q1+t(Q1)
          diag(Q1)<-(optim.obj$par[1:s])
          
          
          #### Check for positiv definitness ########
          for (ttt in 0:100)
          {
            Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
            Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
            Q_solvetest<-try(solve(Q1))
            if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
              break
          }
        }
      }
      
      IC_stop[l]<-sum(y*log(Mu)+(1-y)*log(1-Mu))
      
      ############ AIC ##############
      if (control$sel.method=="aic")
        IC_stop[l]<-(-IC_stop[l])+defre
      
      ############ BIC ###############
      if (control$sel.method=="bic")
        IC_stop[l]<-(-IC_stop[l])+0.5*defre*log(n)
      
      Q[,,l+1]<-Q1
      
      kritval<-sqrt(sum((Delta[l-1,]-Delta[l,])^2))/sqrt(sum(Delta[l,]^2))
      
      
      if(kritval<1e-8)
        break
      
      
      kritval2<-abs(IC_stop[l-1]-IC_stop[l])/abs(IC_stop[l])
      
      if(kritval2<1e-8)
        break
      
      
    }}
  if (l<control$steps)
  {
    for (wj in (l+1):control$steps)
    {
      Delta[wj,]<-Delta[l,]
      Q[,,wj+1]<-Q1
      IC[wj,]<-IC[l,]
      IC_stop[wj]<-IC_stop[l]
    }}
  
  ##########
  opt<-match(min(IC_stop),IC_stop)
  
  if (control$OPT==TRUE)
  {
    Delta_neu<-Delta[opt,]
  }else{
    Delta_neu<-Delta[l,]
  }
  
  Qfinal<-Q[,,opt+1]
  
  aaa<-!is.element(Delta_neu[1:(lin+m_alt)],0)
  
  glmm_final<-try(glmm_final_seq(y,Z_fastalles[,aaa],W,k,q_start=Qfinal,
                                 Delta_start=Delta_neu[c(aaa,rep(T,n*s))],s,
                                 steps=control$maxIter,method=method,kat=kat,
                                 katvar=katvar,print.iter.final=control$print.iter.final,
                                 eps.final=control$eps.final))
        
        if(class(glmm_final)=="try-error" || glmm_final$opt>(control$maxIter-5))
        {
          glmm_final<-try(glmm_final_seq(y,Z_fastalles[,aaa],W,k,q_start=q_start,
                                         Delta_start=Delta_start[c(aaa,rep(T,n*s))],
                                         s,steps=2000,method=method,kat=kat,
                                         katvar=katvar,print.iter.final=control$print.iter.final,
                                         eps.final=control$eps.final))
          
          if(class(glmm_final)=="try-error" || glmm_final$opt>1990)
          {
            cat("Warning:\n")
            cat("Final Fisher scoring reestimation did not converge!")
          }}  
  
  Delta_neu[c(aaa,rep(T,n*s))]<-glmm_final$Delta
  
  Standard_errors<-rep(0,length(Delta_neu))
  Standard_errors[c(aaa,rep(T,n*s))]<-glmm_final$Standard_errors
  
  Qfinal<-glmm_final$Q
  
  Eta_opt<-Z_alles%*%Delta_neu
  Mu_opt<-as.vector(family$linkinv(Eta_opt))
  
  if(s==1)
    Qfinal<-sqrt(Qfinal)
  
  
  if(!is.matrix(Qfinal))
    Qfinal<-as.matrix(Qfinal)
  colnames(Qfinal)<-random.labels
  rownames(Qfinal)<-random.labels
  
  names(Delta_neu)[1:lin]<-colnames(X)[1:lin]
  names(Standard_errors)[1:lin]<-colnames(X)[1:lin]
  
  
  names(Delta_neu)[(lin+1):(lin+m_alt)]<-colnames(U)
  names(Standard_errors)[(lin+1):(lin+m_alt)]<-colnames(U)
  
  
  permut<-numeric()
  for(ys in 1:length(old.names2[!is.element(old.names2,control$katvar)]))
    permut<-c(permut,match(old.names2[!is.element(old.names2,control$katvar)][ys],names(Delta_neu)))
  
  Delta_neu[((SQ$kat-1)*(1+length(control$katvar))+1):((SQ$kat-1)*(1+length(control$katvar))+length(old.names2[!is.element(old.names2,control$katvar)]))]<-Delta_neu[permut]
  Standard_errors[((SQ$kat-1)*(1+length(control$katvar))+1):((SQ$kat-1)*(1+length(control$katvar))+length(old.names2[!is.element(old.names2,control$katvar)]))]<-Standard_errors[permut]
  
  names(Delta_neu)[((SQ$kat-1)*(1+length(control$katvar))+1):((SQ$kat-1)*(1+length(control$katvar))+length(old.names2[!is.element(old.names2,control$katvar)]))]<-old.names2[!is.element(old.names2,control$katvar)]
  names(Standard_errors)[((SQ$kat-1)*(1+length(control$katvar))+1):((SQ$kat-1)*(1+length(control$katvar))+length(old.names2[!is.element(old.names2,control$katvar)]))]<-old.names2[!is.element(old.names2,control$katvar)]
  
  
  Delta[,((SQ$kat-1)*(1+length(control$katvar))+1):((SQ$kat-1)*(1+length(control$katvar))+length(old.names2[!is.element(old.names2,control$katvar)]))]<-Delta[,permut]
  
  names(Delta_neu)[(lin+m_alt+1):(lin+m_alt+n*s)]<-colnames(W)
  names(Standard_errors)[(lin+m_alt+1):(lin+m_alt+n*s)]<-colnames(W)
  colnames(Delta)<-names(Delta_neu)
  
  colnames(IC)<-very.old.names
  
  comp<-character()
  for(oi in 1:length(komp))
    comp<-c(comp,very.old.names[komp[oi]])
  
  ret.obj=list()
  ret.obj$IC<-IC
  ret.obj$IC_sel<-IC_stop
  ret.obj$opt<-opt
  ret.obj$Deltamatrix<-Delta
  ret.obj$ranef<-Delta_neu[(lin+m_alt+1):(lin+m_alt+n*s)]
  ret.obj$coefficients<-Delta_neu[1:(lin+m_alt)]
  ret.obj$fixerror<-Standard_errors[1:(lin+m_alt)]
  ret.obj$ranerror<-Standard_errors[(lin+m_alt+1):(lin+m_alt+n*s)]
  ret.obj$components<-comp
  ret.obj$Q_long<-Q
  ret.obj$Q<-Qfinal
  ret.obj$y_hat<-Mu_opt
  ret.obj$HatMatrix<-FinalHat[[opt]]
  ret.obj$family<-family
  ret.obj$fix<-fix
  ret.obj$newrndfrml<-newrndfrml
  ret.obj$subject<-names(rnd)[1]
  ret.obj$k<-k
  ret.obj$kat<-kat
  ret.obj$katvar<-katvar
  ret.obj$data<-data
  ret.obj$factor.names<-factor.names
  ret.obj$factor.list<-factor.list
  return(ret.obj)
}

####################################################################################
####################################################################################

OrdinalBoost <- function(fix=formula, rnd=formula, data,  model="sequential", control=list()) UseMethod("OrdinalBoost")

OrdinalBoost.formula <- function(fix,rnd,...)
{
est <- est.OrdinalBoost(fix,rnd,...)
est$fitted.values <- est$y_hat
est$StdDev <- est$Q
est$call <- match.call()
class(est) <- "OrdinalBoost"
est
}


print.OrdinalBoost <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\nFixed Effects:\n")
cat("\nCoefficients:\n")
print(x$coefficients)
cat("\nRandom Effects:\n")
cat("\nStdDev:\n")
print(x$StdDev)
}




summary.OrdinalBoost <- function(object, ...)
{
se <- object$fixerror
zval <- coefficients(object) / se
TAB <- cbind(Estimate = coefficients(object),
StdErr = se,
z.value = zval,
p.value = 2*pnorm(-abs(zval)))
res <- list(call=object$call,
coefficients=TAB,StdDev=object$StdDev)
class(res) <- "summary.OrdinalBoost"
res
}


print.summary.OrdinalBoost <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\n")
cat("\nFixed Effects:\n")
cat("\nCoefficients:\n")
printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
cat("\nRandom Effects:\n")
cat("\nStdDev:\n")
print(x$StdDev)
}


predict.OrdinalBoost <- function(object,newdata=NULL,...)
{
family=binomial(link="logit")
kat<-object$kat
katvar<-object$katvar
model<-object$model
krit.random<-TRUE
  
subj.new<-levels(as.factor(newdata[,object$subject]))
subj.old<-levels(as.factor(object$data[,object$subject]))
subj.test<-is.element(subj.new,subj.old)
subj.ok<-subj.new[subj.test]

if(!is.null(newdata))
{
rand.ok<-is.element(newdata[,object$subject],subj.ok)

  if(length(object$factor.names>0))
  {
  for (i in 1:length(object$factor.names))
  {
  newdata[,object$factor.names[i]]<-factor(newdata[,object$factor.names[i]],levels=object$factor.list[[i]])
  }}

X <- model.matrix(object$fix, newdata)

X<-cbind(X,rand.ok)

krit.random<-!all(!is.element(subj.new,subj.old))
  if(krit.random)
  {
  W_start <- model.matrix(formula(object$newrndfrml), newdata)
  }else{
  W_start <- NULL}
}else{
X <- model.matrix(object$fix, object$data)
W_start <- model.matrix(formula(object$newrndfrml), object$data)
}

X<-X[,-1]



rnlabels<-terms(formula(object$newrndfrml))
random.labels<-attr(rnlabels,"term.labels")
s<-length(random.labels)
if(!is.null(newdata))
{
k<-table(newdata[,colnames(newdata)==(object$subject)])
}else{
k<-object$k
}
n<-length(k)
N<-dim(X)[1]

if(s>1)
{for (i in 2:s)
subj.test<-cbind(subj.test,subj.test)
subj.test<-as.vector(t(subj.test))
}


if(krit.random)
{
if(s>1)
{
W<-W_start[,seq(from=1,to=1+(s-1)*n,by=n)]
for (i in 2:n)
W<-cbind(W,W_start[,seq(from=i,to=i+(s-1)*n,by=n)])
}else{
W<-W_start
}
if(length(subj.test)>0)
W<-W[,subj.test]
}else{
W<-NULL
}

X_neu<-numeric()

if(is.null(katvar))
{
   for (kf in 1:N)
   {
   XX<-numeric()
   for (i in 1:(kat-1))
   XX<-rbind(XX,c(X[kf,],W[kf,]))
   X_neu<-rbind(X_neu,cbind(diag(kat-1),XX))
   }
}else{
X_kat<-cbind(1,X[,katvar])
colnames(X_kat)<-c("Icept",katvar)
namesvec<-colnames(X)[!is.element(colnames(X),katvar)]
X<-X[,namesvec]

   for (kf in 1:N)
   {
   XX<-numeric()
   for (i in 1:(kat-1))
   XX<-rbind(XX,c(X[kf,],W[kf,]))
   X_kat_akt<-t(X_kat[kf,])
   if((kat-1)>1)
   {
   for (vf in 2:(kat-1))
   X_kat_akt<-adiag(X_kat_akt,t(X_kat[kf,]))
   }
   X_neu<-rbind(X_neu,cbind(X_kat_akt,XX))
   }
}

if(krit.random)
{
  if(!is.null(newdata))
  {
  rand.ok<-as.logical(X_neu[,is.element(attr(X_neu,"dimnames")[[2]],"rand.ok")])
  X_neu<-X_neu[,!is.element(attr(X_neu,"dimnames")[[2]],"rand.ok")]
  Mu<-family$linkinv(X_neu[,1:length(object$coef)]%*%object$coef)
  Mu[rand.ok]<-family$linkinv(X_neu[rand.ok,]%*%c(object$coef,object$ranef[match(colnames(W),names(object$ranef))]))
  }else{
  Mu<-family$linkinv(X_neu%*%c(object$coef,object$ranef))
  }
}else{
Mu<-family$linkinv(X_neu%*%object$coef)
}
y<-numeric()
y0<-numeric()


if(model=="sequential")
{
  for (zt in 1:N)
  {
  y0<-Mu[(zt-1)*(kat-1)+1]
  for(ru in 2:(kat-1))
  y0<-c(y0,Mu[(zt-1)*(kat-1)+ru]*prod(1-Mu[((zt-1)*(kat-1)+1):((zt-1)*(kat-1)+ru-1)]))
  y0<-c(y0,prod(1-Mu[((zt-1)*(kat-1)+1):((zt-1)*(kat-1)+kat-1)]))
  y<-rbind(y,y0)
  }
}else{
  for (zt in 1:N)
  {
  y0<-Mu[(zt-1)*(kat-1)+1]
  for(ru in 2:(kat-1))
  y0<-c(y0,Mu[(zt-1)*(kat-1)+ru]-Mu[(zt-1)*(kat-1)+ru-1])
  y0<-c(y0,1-Mu[(zt-1)*(kat-1)+kat-1])
  y<-rbind(y,y0)
  }
}
colnames(y)<-paste("Category",c(1:kat),sep="")
rownames(y)<-c(1:N)
y
}


