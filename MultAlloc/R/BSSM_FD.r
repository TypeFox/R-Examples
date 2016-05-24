BSSM_FD<-function(Nh,Sh2j,Yj,Cust=NULL,nmin=2,ch=NULL,w=NULL,certain=FALSE)
{
produced_cvs<-NULL
H<-length(Nh)
z<-cumsum(Nh)
z=matrix(sort(c(1,z,z+1)[-(2*H+1)]),ncol=2,byrow=TRUE)
if (is.matrix(Sh2j)==FALSE) {Sh2j<-t(as.matrix(Sh2j))}
number_variables=dim(Sh2j)[1]
CYj<-NULL
cxh<-NULL
if (length(w)==0) {w=rep(1/number_variables,number_variables)}
    for(j in 1:number_variables) {CYj<-c(CYj,w[j]*(1/(Yj[j]^2)))}
    cvw<-NULL
    for(i in 1:H) {cvw<-c(cvw,sum(CYj*Sh2j[,i]))}
    for(i in 1:H) {cxh<-c(cxh,cvw[i]*(Nh[i]^2/(1:Nh[i])))}
    fobj<-cxh
    restriction1=matrix(rep(rep(0,sum(Nh)),H),nrow=H)
    for(i in 1:H) {restriction1[i,z[i,1]:z[i,2]]=1}
    restriction2<-NULL
    values<-NULL
    if (length(ch)>0)
      {for(i in 1:H)
        {restriction2<-c(restriction2,ch[i]*c(1:Nh[i]))
         values<-c(values,c(1:Nh[i]))
        }
      }
    else
      {for(i in 1:H) {restriction2<-c(restriction2,c(1:Nh[i]))}
        values<-restriction2
      }
     restriction3<-restriction1
     for(i in 1:H) {restriction3[i,z[i,1]:z[i,2]]<-c(1:Nh[i])}
     A<-rbind(restriction1,restriction2,restriction3)
     b<-c(rep(1,H),Cust,rep(nmin,times=H))
     desig<-c(rep("==",H),"<=",rep(">=",times=H))
     if (certain==TRUE)
      {desig[length(desig)]<-"=="
       b[length(b)]<-Nh[H]
      }
     tempo<-proc.time()
     x=Rglpk::Rglpk_solve_LP(fobj,A,desig,b,types=rep("B",length(fobj)),max=FALSE)
     tempo<-proc.time()-tempo
     nh=x$solution*values
     nh=nh[nh>0]
     cvh<-NULL
     for(j in 1:number_variables)
        {cvh<-c(cvh,sum(Sh2j[j,]*(Nh^2/nh)*(1-nh/Nh)))
        }
     produced_cvs<-round(sqrt(cvh)/Yj,5)
return(list(nh=nh,n=sum(nh),cvs=produced_cvs,time_cpu=tempo[3]))
}
