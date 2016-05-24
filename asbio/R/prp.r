#############################################################################################################
################################################ prp ########################################################
#############################################################################################################

prp<-function(Time,S.X,S.Y,N.X,N.Y,habitat=NULL,near.angle=NULL, F.0.NA=TRUE){


pairw<-length(S.X)-1
  
  Euc<-function(x,y)sqrt((x[1]-y[1])^2+(x[2]-y[2])^2) ##Euclidean distance

if(!is.null(near.angle))
  {Find.quad<-function(x) 
    {for(i in 1:length(S.X))
      {if(x[i]>=0&x[i]<=90){x[i]=1}
      if(x[i]<=180&x[i]>90){x[i]=2}
      if(x[i]<(-90)&x[i]>=(-180)){x[i]=3}
      if(x[i]<0&x[i]>=(-90)){x[i]=4}
      }
    x
    }
  near.angle<-Find.quad(near.angle)
  }

######################################## Lines ###########################################

Deer<-cbind(S.X,S.Y)
Near<-cbind(N.X,N.Y)

A<-matrix(ncol=1,nrow=pairw)
B<-matrix(ncol=1,nrow=pairw)
C<-matrix(ncol=1,nrow=pairw)
D<-matrix(ncol=1,nrow=pairw)
F<-matrix(ncol=1,nrow=pairw)

for(i in 1:(pairw))
  {A[i]<-Euc(Deer[i,],Near[i,])
  B[i]<-Euc(Deer[i+1,],Near[i+1,])
  C[i]<-Euc(Deer[i,],Deer[i+1,])
  D[i]<-Euc(Deer[i+1,],Near[i,])
  F[i]<-Euc(Near[i,],Near[i+1,])
  }


######################################## Angles ##########################################

kappa<-acos((C^2+A^2-D^2)/(2*C*A))*180/pi
gamma<-acos((B^2+F^2-D^2)/(2*B*F))*180/pi

######################################## Index ##########################################

m.C<-matrix(ncol=1,nrow=pairw)
m.F<-matrix(ncol=1,nrow=pairw)
b.C<-matrix(ncol=1,nrow=pairw)
b.F<-matrix(ncol=1,nrow=pairw)
int.x.coord<-matrix(ncol=1,nrow=pairw)
int.y.coord<-matrix(ncol=1,nrow=pairw)

for(i in 1:pairw)
  {if(F[i]!=0&C[i]!=0){
  m.C[i]<-(S.Y[i+1]-S.Y[i])/(S.X[i+1]-S.X[i])
  m.F[i]<-(N.Y[i+1]-N.Y[i])/(N.X[i+1]-N.X[i])}
  if(C[i]==0&F[i]==0){ ## No movement
  m.C[i]=-9999
  m.F[i]=-9999}
  if(F[i]==0&C[i]!=0){
  m.C[i]=-9998
  m.F[i]=-9998}
  if(F[i]!=0&C[i]==0)stop("F!=0 but C=0")
  }

for(i in 1:pairw)#### Infinite slopes??  No problem
  {if(m.C[i]!=-9999&m.C[i]!=-9998)
    {if(m.C[i]==m.F[i]){##Parallel
    m.C[i]= -9997
    m.F[i]= -9997}
    if(m.C[i]!=m.F[i]){##Not Parallel
      {if((m.C[i]==Inf|m.C[i]==-Inf)&(m.F[i]!=Inf|m.F[i]!=-Inf)){
      int.x.coord[i]<-S.X[i]
      b.F[i]<-N.Y[i]-N.X[i]*m.F[i]
      int.y.coord[i]<-int.x.coord[i]*m.F[i]+b.F[i]}
      if((m.C[i]!=Inf|m.C[i]!=-Inf)&(m.F[i]==Inf|m.F[i]==-Inf)){
      int.x.coord[i]<-N.X[i]
      b.C[i]<-S.Y[i]-S.X[i]*m.C[i]
      int.y.coord[i]<-int.x.coord[i]*m.C[i]+b.C[i]}
      if((m.C[i]!=Inf&m.C[i]!=-Inf)&(m.F[i]!=Inf&m.C[i]!=-Inf)){##No infinite slopes
      b.C[i]<-S.Y[i]-S.X[i]*m.C[i]
      b.F[i]<-N.Y[i]-N.X[i]*m.F[i]
      int.x.coord[i]<-(b.C[i]-b.F[i])/(m.F[i]-m.C[i])
      int.y.coord[i]<-int.x.coord[i]*m.F[i]+b.F[i]
        }
      }
    }
  }
  if(m.C[i]==-9999){## No movement
  int.x.coord[i]=-9999
  int.y.coord[i]=-9999
  }
  if(m.C[i]==-9998){## Same nearest neighbor
  int.x.coord[i]=-9998
  int.y.coord[i]=-9998
  }
}
  
C.ext<-matrix(ncol=1,nrow=pairw)
F.ext<-matrix(ncol=1,nrow=pairw)
last<-matrix(ncol=1,nrow=pairw)

for(i in 1:pairw)
  {if(m.C[i]!=-9999&m.F[i]!=-9999){##Find which side delta is on and calculte S and T
  C.ext[i]<-max(c(Euc(c(int.x.coord[i],int.y.coord[i]),c(S.X[i],S.Y[i])),Euc(c(int.x.coord[i],int.y.coord[i]),
  c(S.X[i+1],S.Y[i+1]))))
  F.ext[i]<-max(c(Euc(c(int.x.coord[i],int.y.coord[i]),c(N.X[i],N.Y[i])),Euc(c(int.x.coord[i],int.y.coord[i]),
  c(N.X[i+1],N.Y[i+1]))))
  last[i]<-max(A[i],B[i])}##opposite side of delta
  if(m.C[i]==-9999|m.F[i]==-9999){
  C.ext[i]=-9999
  F.ext[i]=-9999
  last[i]=-9999}
  }



delta<-matrix(ncol=1,nrow=pairw)
for(i in 1:pairw)
  {if(F[i]==0&is.null(habitat)){
    if(F[i]==0&C[i]!=0&F.0.NA==FALSE)delta[i]=90 ## F = 0; let delta = 90
    if(F[i]==0&C[i]!=0&F.0.NA==TRUE)delta[i]=NA   ## F = 0; let delta = NA
  }
  if(F[i]!=0|(F[i]==0&!is.null(habitat))){
    if(F[i]==0&!is.null(habitat)){
      if(habitat[i]!=habitat[i+1])delta[i]=90
      if(habitat[i]==habitat[i+1]){
        if(C[i]!=0&F.0.NA==FALSE)delta[i]=90
        if(C[i]!=0&F.0.NA==TRUE)delta[i]=NA}
        }
    if(m.C[i]!=-9999&m.F[i]!=-9999&m.C[i]!=-9998){
    delta[i]<-acos((C.ext[i]^2+F.ext[i]^2-last[i]^2)/(2*C.ext[i]*F.ext[i]))*180/pi}
    if(m.C[i]==-9999|m.F[i]==-9999)delta[i]=NA
    if(m.C[i]==-9997|m.F[i]==-9997)delta[i]=0
    }
  }
  
perp<-matrix(ncol=1,nrow=pairw)

for(i in 1:pairw)
  {if(is.na(delta[i])){delta[i]=NA}
  if(!is.na(delta[i])){
  if(delta[i]<=90){perp[i]=delta[i]/90}
  if(delta[i]>90&delta[i]<=135){perp[i]=(90-(delta[i]-90))/90}
  if(delta[i]>135&delta[i]<=180){perp[i]=(delta[i]-90)/90}}
  }

##################################################################################################

if(!is.null(near.angle)|!is.null(habitat)){ 
  
  for(i in 1:pairw)## Adjustments for potential problems with switching boundaries in a patch
    {if(!is.na(perp[i])){
    perp[i]<-ifelse(near.angle[i]!=near.angle[i+1]&habitat[i]==habitat[i+1],NA,perp[i])
    delta[i]<-ifelse(near.angle[i]!=near.angle[i+1]&habitat[i]==habitat[i+1],NA,delta[i])}
    }
    
  crossing<-matrix(ncol=1,nrow=pairw)
  for(i in 1:pairw){
  crossing[i]<-ifelse(habitat[i]!=habitat[i+1],1,0)}
  
  hab.type<-matrix(ncol=1,nrow=pairw)
  for(i in 1:pairw){
  if(habitat[i]==habitat[i+1]){hab.type[i]=habitat[i]}
  if(habitat[i]!=habitat[i+1]){hab.type[i]="Crossing"}
  }
  

  border.change<-matrix(ncol=1,nrow=pairw)
  for(i in 1:pairw){
  border.change[i]<-ifelse(near.angle[i]!=near.angle[i+1],"Y","N")}
  }

############################################ Results ################################################

  end.time<-matrix(ncol=1,nrow=pairw)
  for(i in 1:pairw){
  end.time[i]<-round(Time[i+1],0)}
  
res<-list()

  angles<-round(cbind(kappa,gamma,delta),3)
  colnames(angles)<-c("Kappa","Gamma","Delta")
  rownames(angles)<-end.time
  res$angles<-angles
  
  lines<-round(cbind(A,B,C,D,F),3)
  colnames(lines)<-c("A","B","C","D","F")
  rownames(lines)<-end.time
  res$lines<-lines
  
  n.perp<-length(!is.na(perp))
  SD.perp<-sapply(perp,function(x)sd(x,na.rm=TRUE))
  SE.perp<-SD.perp/sqrt(n.perp)
  
if(!is.null(habitat))
  {res$moment.by.moment<-data.frame(matrix(nrow=pairw,ncol=5,data=c(end.time,round(perp,4),round(delta,4),hab.type,border.change),dimnames=list(seq(1,pairw),c("End.time","Eta Index","Delta","Habitat","Brdr chng"))))
  }
if(is.null(habitat))
  {res$moment.by.moment<-data.frame(matrix(nrow=pairw,ncol=3,data=c(end.time,round(perp,4),round(delta,4)),dimnames=list(seq(1,pairw),c("End.time","Eta Index","Delta"))))
  }
res$P.summary<-matrix(nrow=1,ncol=3,data=c(mean(perp,na.rm=TRUE),SE.perp,n.perp),dimnames=list(c(""),c("Mean.Eta","SE.Eta","n")))
if(!is.null(near.angle)|!is.null(habitat))
  {
  x<-length(crossing[crossing==1]) 
  p.wilson<-(x+2)/(pairw+4)
  res$crossing.summary<-matrix(nrow=2,ncol=2,data=c(mean(crossing),sqrt((mean(crossing)*(1-mean(crossing)))/pairw),p.wilson,sqrt((p.wilson*(1-p.wilson))/(pairw+4))),dimnames=list(c("p.hat","S.p.hat"),c("Binomial estimates","Wilson estimates")))#estimators for binomial dist
  }
if(is.null(near.angle)|is.null(habitat))res$crossing.summary<-NA
class(res)<-"prp.index"
res
}
                               
