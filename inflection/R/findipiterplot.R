findipiterplot <-
  function(x,y,index)
  {
    ESE<-c();AESE<-c();BESE<-c();
    EDE<-c();AEDE<-c();BEDE<-c();
    xpes<-c();ypes<-c();npes<-c();
    xped<-c();yped<-c();nped<-c();
    A<-findiplist(x,y,index);
    B<-A;
    x1<-x;y1<-y;
    x2<-x;y2<-y;
    ESE<-c(ESE,A[1,3]);AESE<-c(AESE,A[1,3])
    EDE<-c(EDE,A[2,3]);AEDE<-c(AEDE,A[2,3])
    #First run:
    ans1<-matrix(c(A[1,1],A[1,2],A[1,3],A[2,1],A[2,2],A[2,3]),nrow=2,ncol=3,byrow=TRUE,dimnames=list(c("ESE","EDE"),c("i1","i2","chi_S,D")))
    #ESE iterations:
    i<-0
    while (is.na(A[1,3])==FALSE)  
    {
      if ((A[1,2]>A[1,1]+1)==TRUE)
      {
        i<-i+1
        xpes<-c(xpes,x1);ypes<-c(ypes,y1);npes<-c(npes,dim(x1)[1]);
        x1<-cbind(x1[A[1,1]:A[1,2],1])
        y1<-cbind(y1[A[1,1]:A[1,2],1])
        A<-findiplist(x1,y1,index)
        if ((A[1,2]>A[1,1]+1)==TRUE)
        {
          ESE<-c(ESE,A[1,3]);AESE<-c(AESE,A[1,3]);
          EDE<-c(EDE,A[2,3]);BEDE<-c(BEDE,A[2,3]);         
        }      
      }
      else
        break    
    }
    #EDE iterations:
    j<-0
    while (is.na(B[2,3])==FALSE)
    {
      if ((B[2,2]>B[2,1]+1)==TRUE)
      {
        j<-j+1
        xped<-c(xped,x2);yped<-c(yped,y2);nped<-c(nped,dim(x2)[1]);
        x2<-cbind(x2[B[2,1]:B[2,2],1])
        y2<-cbind(y2[B[2,1]:B[2,2],1])
        B<-findiplist(x2,y2,index)
        if ((B[2,2]>B[2,1]+1)==TRUE)
        {
          if(is.na(B[2,3])==FALSE)
          {
            ESE<-c(ESE,B[1,3]);BESE<-c(BESE,B[1,3]);
            EDE<-c(EDE,B[2,3]);AEDE<-c(AEDE,B[2,3]);
          }
          else
            break
        }      
      }
      else
        break    
    }
    #Building answer, descriptive statistics and 95% confidence intervals:
    ESE=rbind(ESE);AESE=rbind(AESE);BESE=rbind(BESE);
    EDE=rbind(EDE);AEDE=rbind(AEDE);BEDE=rbind(BEDE);
    ans=new.env()
    ans$first=ans1
    ans$aese=AESE;rownames(ans$aese)=c("ESE iterations");
    #Cross iterative run of ESE - 95% confidence interval:
{
  if(is.null(BESE)==FALSE)
  {ans$bese=BESE;
   rownames(ans$bese)=c("ESE results from EDE iters");
   bmesm=apply(BESE,1,mean,na.rm = TRUE);
   bsesm=apply(BESE,1,sd,na.rm = TRUE);
   beresm=qt(0.975,df=length(BESE)-1)*bsesm/sqrt(length(BESE));
   blciesm=bmesm-beresm;brciesm=bmesm+beresm;
   ans$besmout=cbind(bmesm,bsesm,blciesm,brciesm)
   colnames(ans$besmout)=c("mean","sdev","95%(l)","95%(r)")
   rownames(ans$besmout)=c("ESE from EDE iters")  
  }
  else{ans$bese=NA}
}
    ans$esm=ESE;rownames(ans$esm)=c("ESE all iterations");
    ans$aede=AEDE;rownames(ans$aede)=c("EDE iterations");
    #Cross iterative run of EDE - 95% confidence interval:
{
  if(is.null(BEDE)==FALSE)
  {ans$bede=BEDE;
   rownames(ans$bede)=c("EDE results from ESE iters");
   bmedm=apply(BEDE,1,mean,na.rm = TRUE);
   bsedm=apply(BEDE,1,sd,na.rm = TRUE);
   beredm=qt(0.975,df=length(BEDE)-1)*bsedm/sqrt(length(BEDE));
   blciedm=bmedm-beredm;brciedm=bmedm+beredm;
   ans$bedmout=cbind(bmedm,bsedm,blciedm,brciedm)
   colnames(ans$bedmout)=c("mean","sdev","95%(l)","95%(r)")
   rownames(ans$bedmout)=c("EDE from ESE iters")    
  }
  else
  {ans$bede=NA}  
}
    #Iterative run of methods - 95% confidence intervals:
    ans$edm=EDE;rownames(ans$edm)=c("EDE all iterations");
    amesm=apply(AESE,1,mean,na.rm = TRUE);
    amedm=apply(AEDE,1,mean,na.rm = TRUE);
    asesm=apply(AESE,1,sd,na.rm = TRUE);
    asedm=apply(AEDE,1,sd,na.rm = TRUE);
    aeresm=qt(0.975,df=length(AESE)-1)*asesm/sqrt(length(AESE));
    aeredm=qt(0.975,df=length(AEDE)-1)*asedm/sqrt(length(AEDE));
    alciesm=amesm-aeresm;arciesm=amesm+aeresm;
    alciedm=amedm-aeredm;arciedm=amedm+aeredm;
    ans$aesmout=cbind(amesm,asesm,alciesm,arciesm)
    colnames(ans$aesmout)=c("mean","sdev","95%(l)","95%(r)")
    rownames(ans$aesmout)=c("ESE method")
    ans$aedmout=cbind(amedm,asedm,alciedm,arciedm)
    colnames(ans$aedmout)=c("mean","sdev","95%(l)","95%(r)")
    rownames(ans$aedmout)=c("EDE method")
    mesm=apply(ESE,1,mean,na.rm = TRUE);
    medm=apply(EDE,1,mean,na.rm = TRUE);
    sesm=apply(ESE,1,sd,na.rm = TRUE);
    sedm=apply(EDE,1,sd,na.rm = TRUE);
    eresm=qt(0.975,df=length(ESE)-1)*sesm/sqrt(length(ESE));
    eredm=qt(0.975,df=length(EDE)-1)*sedm/sqrt(length(EDE));
    lciesm=mesm-eresm;rciesm=mesm+eresm;
    lciedm=medm-eredm;rciedm=medm+eredm;
    ans$esmout=cbind(mesm,sesm,lciesm,rciesm)
    colnames(ans$esmout)=c("mean","sdev","95%(l)","95%(r)")
    rownames(ans$esmout)=c("ESE method, all results")
    ans$edmout=cbind(medm,sedm,lciedm,rciedm)
    colnames(ans$edmout)=c("mean","sdev","95%(l)","95%(r)")
    rownames(ans$edmout)=c("EDE method, all results")
    #Estimation of inflection point using results from both ESE and EDE methods:
    IP=rbind(c(ESE,EDE));
    mip=apply(IP,1,mean,na.rm = TRUE);
    sip=apply(IP,1,sd,na.rm = TRUE);
    erip=qt(0.975,df=length(IP)-1)*sip/sqrt(length(IP));
    lcip=mip-erip;rcip=mip+erip;
    ans$ipall=cbind(mip,sip,lcip,rcip)
    colnames(ans$ipall)=c("mean","stdev","95%(l)","95%(r)")
    rownames(ans$ipall)=c("All methods, all results")
    nps<-length(npes);
    npd<-length(nped);
    xsm<-cbind(xpes);ysm<-cbind(ypes);
    xdm<-cbind(xped);ydm<-cbind(yped);
    ns=cbind(npes);
    nd=cbind(nped);
    dev.new()
    par(mfrow=c(2,ceiling(nps/2)))
    n0<-0
    dfs<-c();
    for (i in 1:nps)
    {
      x1<-cbind(xsm[(n0+1):(n0+ns[1])]);y1<-cbind(ysm[(n0+1):(n0+ns[1])]);
      x<-cbind(x1[!is.na(x1)]);y<-cbind(y1[!is.na(y1)]);dfs=c(dfs,data.frame(x=x,y=y));
      plot(x,y,xlab=paste("x",i),ylab=paste("y",i))
      abline(v=AESE[i],lty=2,col="blue")
      title(paste("ESE iter",i),sub=paste("ip=",AESE[i]))
      n0<-n0+ns[i];
    }
    dev.new()
    par(mfrow=c(2,ceiling(npd/2)))
    n0<-0
    dfd<-c();
    for (i in 1:npd)
    {
      x1<-cbind(xdm[(n0+1):(n0+nd[1])]);y1<-cbind(ydm[(n0+1):(n0+nd[1])]);
      x<-cbind(x1[!is.na(x1)]);y<-cbind(y1[!is.na(y1)]);dfd=c(dfd,data.frame(x=x,y=y));
      plot(x,y,xlab=paste("x",i),ylab=paste("y",i))
      abline(v=AEDE[i],lty=2,col="blue")
      title(paste("EDE iter",i),sub=paste("ip=",AEDE[i]))
      n0<-n0+nd[i];
    }
    ans$xysl=dfs;
    ans$xydl=dfd;
    print(matrix(c(mip,sip,lcip,rcip),nrow=1,ncol=4,byrow=TRUE,dimnames=list(c("ip all methods"),c("mean","sdev","95%(l)","95%(r)"))))
    ans
  }