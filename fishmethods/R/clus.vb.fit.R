clus.vb.fit<-function(len=NULL,age=NULL,cluster=NULL,nboot=1000,
                    control=list(maxiter=10000,
                    minFactor=1/1024,tol=1e-5)){

if(length(age)!=length(len)) stop ("Vectors of different lengths")
if(length(cluster)!=length(len)) stop ("Vectors of different lengths")
if(length(cluster)!=length(age)) stop ("Vectors of different lengths")
if(mode(cluster)!="numeric") stop ("cluster must be a numeric variable")

ins<-as.data.frame(cbind(len,age,cluster))
ins<-ins[!is.na(ins$len) & !is.na(ins$age) & !is.na(ins$cluster),]
ins<-ins[order(ins$cluster),]
   nclus<-length(as.character(as.data.frame(table(ins$cluster))$Var1))
   listclus<-as.character(as.data.frame(table(ins$cluster))$Var1)   
  parms<-data.frame(NULL)
  for(i in 1:nboot){
    #boot strap by cluster
    x<-data.frame(NULL)
    out<-NULL
    temp<-sample(1:nclus,replace=T)
    tempclus<-listclus[temp]
      for(j in 1:nclus){
        out<-ins[ins$cluster==tempclus[j],]
        x<-rbind(x,out)
      }
  #calculate starting values for parameters
  mbt<-NULL;mbxt<-NULL;g2<-NULL
  g2<-aggregate(x$len,list(round(x$age,0)),mean)
  mbt<-g2;mbt[,1]<-mbt[,1]-1
  mbxt<-merge(g2,mbt,by.x=c("Group.1"),by.y=c("Group.1"))
  outboth<-lm(mbxt[,3]~mbxt[,2])
  Kboth<-abs(log(coef(outboth)[2]));Lboth<--coef(outboth)[1]/(coef(outboth)[2]-1)
  dxb<-as.data.frame(cbind(Lboth-g2$x,g2[,1]));dxb<-dxb[dxb[,1]>0,]
  t0b<-(coef(lm(log(dxb[,1])~dxb[,2]))[1]-log(Lboth))/Kboth
  Linf1<-Lboth;K1<-Kboth;t01<-t0b
  #fit model
  vbmodel<-try(nls(len~Linf*(1-exp(-K*(age-t0))),data=x,
              start=c(Linf=Linf1[[1]],K=K1[[1]],t0=t01[[1]]),
              control=control),silent=TRUE)
  if(class(vbmodel)!="try-error"){ 
     parms[i,1]<-coef(vbmodel)[1]
     parms[i,2]<-coef(vbmodel)[2]
     parms[i,3]<-coef(vbmodel)[3] 
     cor<-summary(vbmodel,correlation=TRUE)$correlation
     parms[i,4]<-cor[1,2]
     parms[i,5]<-cor[1,3]
     parms[i,6]<-cor[2,3]
  }
  if(class(vbmodel)=="try-error"){ 
    parms[i,1]<-NA
    parms[i,2]<-NA
    parms[i,3]<-NA
    parms[i,4]<-NA
    parms[i,5]<-NA
    parms[i,6]<-NA
  }
} # boot loop
 #calculate statistics
 means<-colMeans(parms,na.rm=TRUE)
 sds<-apply(parms,2,sd,na.rm=TRUE)
 qLinf<-quantile(parms[,1],c(0.025,0.975),na.rm=TRUE)
 qK<-quantile(parms[,2],c(0.025,0.975),na.rm=TRUE)
 qt0<-quantile(parms[,3],c(0.025,0.975),na.rm=TRUE)
 # output
 nfails<-length(parms[is.na(parms[,1]),1]) 
 title<-c(paste("Parameter coefficients based on ",nboot-nfails," successful fits"))
 ans<-NULL 
 ans$summary<-paste(title)
 ans$results<-matrix(NA,3L,6L)
 ans$results<-rbind(cbind(round(means[1],2),round(sds[1],3),round(qLinf[1],3),round(qLinf[2],3)),
                   cbind(round(means[2],3),round(sds[2],4),round(qK[1],4),round(qK[2],4)),
                   cbind(round(means[3],3),round(sds[3],4),round(qt0[1],4),round(qt0[2],4)))                   
dimnames(ans$results)<-list(cbind("Linf","K","t0"),c("Estimate","SE","95% LCI","95% UCI"))
ans$correlations<-matrix(NA,3L,3L)
ans$correlations<-rbind(cbind(1,round(means[4],3),round(means[5],3)),
                   cbind(round(means[4],3),1,round(means[6],3)),
                   cbind(round(means[5],3),round(means[6],3),1))                   
dimnames(ans$correlations)<-list(cbind("Linf","K","t0"),c("Linf","K","t0"))
return(ans) 
}# function
