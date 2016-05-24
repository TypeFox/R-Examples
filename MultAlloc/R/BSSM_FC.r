BSSM_FC<-function(Nh,Sh2j,Yj,cvt,nmin=2,ch=NULL,certain=FALSE)
{
produced_cvs<-NULL
H<-length(Nh)
z<-cumsum(Nh)
z=matrix(sort(c(1,z,z+1)[-(2*H+1)]),ncol=2,byrow=TRUE)
if (is.matrix(Sh2j)==FALSE) {Sh2j<-t(as.matrix(Sh2j))}
number_variables=dim(Sh2j)[1]
fobj=unlist(apply(as.matrix(Nh),1,function(x) seq(1,x,by=1)))
values<-fobj
if (length(ch)>0)
  {Ch<-cbind(Nh,ch)
   CC<-as.vector(unlist(apply(Ch,1,function(x) rep(x[2],times=x[1]))))
   fobj<-fobj*CC
  }
restriction1=matrix(rep(rep(0,sum(Nh)),H),nrow=H)
for(i in 1:H) {restriction1[i,z[i,1]:z[i,2]]=1}
restriction2=NULL
Cj=NULL
for(i in 1:number_variables) {Cj=rbind(Cj,Nh^2*Sh2j[i,])}
for(i in 1:number_variables) {restriction2=rbind(restriction2,rep(0,sum(Nh)))}
for(j in 1:number_variables)
  {for(i in 1:H)
   {restriction2[j,z[i,1]:z[i,2]]=(rep(Cj[j,i],Nh[i])*(1/c(1:Nh[i]))*(1-c(1:Nh[i])/Nh[i]))}
    restriction2[j,]=restriction2[j,]/(Yj[j]^2*cvt[j]^2)
   }
########################################################################################
desig=c(rep("==",H),rep("<=",number_variables))
b=c(rep(1,H),rep(1.0,number_variables))
restriction3=apply(restriction1,1,function(restriction1) which(restriction1>0))
restriction4=restriction1
if (is.matrix(restriction3)==FALSE)
{for(i in 1:H) {restriction4[i,restriction3[[i]]]=c(1:Nh[i])}
}
else
{for(i in 1:H) {restriction4[i,restriction3[,i]]=c(1:Nh[i])}
}
A=rbind(restriction1,restriction2,restriction4)

desig=c(desig,rep(">=",H))
nmin_modif<-ifelse(Nh<nmin,Nh,nmin)
b=c(b,nmin_modif)
if (certain==TRUE)
{desig[length(desig)]<-"=="
b[length(b)]<-Nh[H]
}
tempo<-proc.time()
x=Rglpk::Rglpk_solve_LP(fobj,A,desig,b,types=rep("B",length(fobj)),max=FALSE)
tempo<-proc.time()-tempo
nh=x$solution[1:sum(Nh)]*values
nh=nh[nh>0]
for(i in 1:number_variables)
   {produced_cvs=cbind(produced_cvs,sqrt(sum(Nh^2*Sh2j[i,]/nh*(1-nh/Nh)))/Yj[i])
   }
produced_cvs<-round(produced_cvs,5)
return(list(n=sum(nh),nh=nh,cvs=produced_cvs,time_cpu=tempo[3]))
}
