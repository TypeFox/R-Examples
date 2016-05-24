
vblrt<-function(len=NULL,age=NULL,group=NULL,error=1,select=1,Linf=NULL,K=NULL,t0=NULL,plottype=0,
                control=list(maxiter=10000,minFactor=1/1024,tol=1e-5)){
   if(is.null(len)) 
         stop ("len is missing") 
   if(is.null(age)) 
         stop ("age is missing") 
   if(is.null(group)) 
         stop ("group is missing.") 
   if(length(age)!=length(len)) stop ("Vectors of different lengths")
   if(nlevels(as.factor(group))>2) stop("Only two groups allowed in data")
   if(select==2 & (is.null(Linf)|is.null(K)|is.null(t0))) stop("User-specified values of Linf, K, and t0 are required")
   
   cat<-as.integer(as.factor(group))-1 
     	x<-as.data.frame(cbind(len,age,cat))
      x<-x[!is.na(x$len) & !is.na(x$age) & !is.na(x$cat),]
    wgt<-NULL
  if(select==1){ 
      m1xt<-NULL;m2xt<-NULL;mbxt<-NULL
      g1<-aggregate(x$len,list(x$cat,trunc(x$age)),mean)
      m1<-g1[g1[,1]==0,];m2<-g1[g1[,1]==1,]
      m1t<-m1[,c(2:3)];m1t[,1]<-m1t[,1]-1; m2t<-m2[,c(2:3)];m2t[,1]<-m2t[,1]-1
      m1xt<-merge(m1,m1t,by.x=c("Group.2"),by.y=c("Group.2"))
      out1<-lm(m1xt[,4]~m1xt[,3])
      K1<-abs(log(coef(out1)[2]));L1<--coef(out1)[1]/(coef(out1)[2]-1)
      dx1<-as.data.frame(cbind(L1-m1$x,m1[,2]));dx1<-dx1[dx1[,1]>0,]
      t01<-(coef(lm(log(dx1[,1])~dx1[,2]))[1]-log(L1))/K1

      m2xt<-merge(m2,m2t,by.x=c("Group.2"),by.y=c("Group.2"))
      out2<-lm(m2xt[,4]~m2xt[,3])
 	K2<-abs(log(coef(out2)[2]));L2<--coef(out2)[1]/(coef(out2)[2]-1)
      dx2<-as.data.frame(cbind(L2-m2$x,m2[,2]));dx2<-dx2[dx2[,1]>0,]
      t02<-(coef(lm(log(dx2[,1])~dx2[,2]))[1]-log(L2))/K2

      g2<-aggregate(x$len,list(round(x$age,0)),mean)
      mbt<-g2;mbt[,1]<-mbt[,1]-1
      mbxt<-merge(g2,mbt,by.x=c("Group.1"),by.y=c("Group.1"))
      outboth<-lm(mbxt[,3]~mbxt[,2])
      Kboth<-abs(log(coef(outboth)[2]));Lboth<--coef(outboth)[1]/(coef(outboth)[2]-1)
 	dxb<-as.data.frame(cbind(Lboth-g2$x,g2[,1]));dxb<-dxb[dxb[,1]>0,]
      t0b<-(coef(lm(log(dxb[,1])~dxb[,2]))[1]-log(Lboth))/Kboth
      Ld<-L2-L1;Kd<-K2-K1;td<-t02-t01
    }
    if(select==2){
         L1<-L2<-Lboth<-Linf; K1<-K2<-Kboth<-K;t01<-t02<-t0b<-t0
         Ld<-L2-L1;Kd<-K2-K1;td<-t02-t01
      }
      if(error==1) x$wgt<-1
 	if(error==2){
         x<-aggregate(x$len,list(x$cat,x$age),mean)
         names(x)<-c("cat","age","len")
         x$wgt<-1
        }
      if(error==3){
         d4<-merge(aggregate(x$len,list(x$cat,x$age),mean),
           aggregate(x$len,list(x$cat,x$age),var),by.y=c("Group.1","Group.2"),
           by.x=c("Group.1","Group.2"))
	   d4<-merge(d4,aggregate(x$len,list(x$cat,x$age),length),by.y=c("Group.1","Group.2"),by.x=c("Group.1","Group.2"))
         names(d4)<-c("cat","age","len","s2","n")
         x<-d4
         x$wgt<-x$n/x$s2 
         if(any(is.na(x$wgt))) stop("At least one age has a single length observation. Need at least two observations to calculate variance." )
        }
  	 Ho<-nls(len~(Linf+ls*cat)*(1-exp(-(K+ks*cat)*(age-(t0+ts*cat)))),data=x,       
        	 weights=wgt,start=list(Linf=L1,ls=Ld,K=K1,ks=Kd,t0=t01,ts=td),
             control=control)
             resid0<-residuals(Ho)
 	 H1<-nls(len~Linf*(1-exp(-(K+ks*cat)*(age-(t0+ts*cat)))),data=x,        
		weights=wgt,start=list(Linf=Lboth,K=K1,ks=Kd,t0=t01,ts=td),
            control=control)
	    resid1<-residuals(H1)  
 	 H2<-nls(len~(Linf+ls*cat)*(1-exp(-K*(age-(t0+ts*cat)))),data=x,       
		weights=wgt,start=list(Linf=L1,ls=Ld,K=Kboth,t0=t01,ts=td),
             control=control)
         resid2<-residuals(H2)
    	 H3<-nls(len~(Linf+ls*cat)*(1-exp(-(K+ks*cat)*(age-t0))),data=x,       
         	weights=wgt,start=list(Linf=L1,ls=Ld,K=K1,ks=Kd,t0=t0b),
             control=control)
         resid3<-residuals(H3)
  	 H4<-nls(len~Linf*(1-exp(-K*(age-t0))),data=x ,      
       	  weights=wgt,start=list(Linf=Lboth,K=Kboth,t0=t0b),
              control=control)
         resid4<-residuals(H4)
  	 RSS<-c(sum(residuals(Ho)^2),sum(residuals(H1)^2),sum(residuals(H2)^2),
             sum(residuals(H3)^2),sum(residuals(H4)^2))

 	 N<-length(residuals(Ho))
  	 X<-round(c(-N*log(RSS[1]/RSS[2]),-N*log(RSS[1]/RSS[3]),-N*log(RSS[1]/RSS[4]),
             -N*log(RSS[1]/RSS[5])),2)
  	 df<-c(length(coef(Ho))-length(coef(H1)),length(coef(Ho))-length(coef(H2)),
       	length(coef(Ho))-length(coef(H3)),length(coef(Ho))-length(coef(H4)))
  	 p<-round(1-pchisq(X,df),3)

      labs<-c("Ho","H1","H2","H3","H4")
      hyp<-c("Linf1=Linf2","K1=K2","t01=t02","Linf1=Linf2,K1=K2,t01=t02")
      labels<-c("Ho vs H1","Ho vs H2","Ho vs H3","Ho vs H4")
      compout<-data.frame(tests=labels,hypothesis=hyp,chisq=X,df=df,p=p)
      rss<-as.data.frame(cbind(labs,RSS));names(rss)<-c("model","rss")
      residuals_all<-as.data.frame(cbind(resid0,resid1,resid2,resid3,resid4))
      nlsout<-list(compout,summary(Ho),summary(H1),summary(H2),summary(H3), summary(H4),
             rss,residuals_all)
      names(nlsout)<-c("results",c(paste("model",labs)),"rss","residuals")
    # Plot observed versus predicted
    if (plottype>0){
      if (plottype==1){
 		par(mfrow=c(3,2))
		plotages<-seq(min(x$age),max(x$age)+1,1)
          # Ho model
	     Linf1<-nlsout$'model Ho'$coefficients[1]
		K1<-nlsout$'model Ho'$coefficients[3]
		t01<-nlsout$'model Ho'$coefficients[5]
		Linf2<-Linf1+nlsout$'model Ho'$coefficients[2]
		K2<-K1+nlsout$'model Ho'$coefficients[4]
		t02<-t01+nlsout$'model Ho'$coefficients[6]
		plot(len~age,data=x[x$cat==0,],main=paste("Ho Model ",levels(group)[1],"=black ", levels(group)[2],"=red"),xlab="Age",ylab="Length",ylim=c(0,max(x$len)))
		points(len~age,data=x[x$cat==1,],col="red")
		lines(Linf1*(1-exp(-K1*(plotages-t01)))~plotages)
		lines(Linf2*(1-exp(-K2*(plotages-t02)))~plotages,col="red")

        # H1 model
	     Linf1<-nlsout$'model H1'$coefficients[1]
		K1<-nlsout$'model H1'$coefficients[2]
		t01<-nlsout$'model H1'$coefficients[4]
		Linf2<-Linf1
		K2<-K1+nlsout$'model H1'$coefficients[3]
		t02<-t01+nlsout$'model H1'$coefficients[5]
		plot(len~age,data=x[x$cat==0,],main=paste("H1 Model ",levels(group)[1],"=black ", levels(group)[2],"=red"),xlab="Age",ylab="Length",ylim=c(0,max(x$len)))
		points(len~age,data=x[x$cat==1,],col="red",xlab="Age",ylab="Length")
		lines(Linf1*(1-exp(-K1*(plotages-t01)))~plotages)
		lines(Linf2*(1-exp(-K2*(plotages-t02)))~plotages,col="red")
    # H2 model
	     Linf1<-nlsout$'model H2'$coefficients[1]
		K1<-nlsout$'model H2'$coefficients[3]
		t01<-nlsout$'model H2'$coefficients[4]
		Linf2<-Linf1+nlsout$'model H2'$coefficients[2]
		K2<-K1
		t02<-t01+nlsout$'model H2'$coefficients[5]
		plot(len~age,data=x[x$cat==0,],main=paste("H2 Model ",levels(group)[1],"=black ", levels(group)[2],"=red"),xlab="Age",ylab="Length",ylim=c(0,max(x$len)))
		points(len~age,data=x[x$cat==1,],col="red")
		lines(Linf1*(1-exp(-K1*(plotages-t01)))~plotages)
		lines(Linf2*(1-exp(-K2*(plotages-t02)))~plotages,col="red")
 # H3 model
	     Linf1<-nlsout$'model H3'$coefficients[1]
		K1<-nlsout$'model H3'$coefficients[3]
		t01<-nlsout$'model H3'$coefficients[5]
		Linf2<-Linf1+nlsout$'model H3'$coefficients[2]
		K2<-K1+nlsout$'model H3'$coefficients[4]
		t02<-t01
	     plot(len~age,data=x[x$cat==0,],main=paste("H3 Model ",levels(group)[1],"=black ", levels(group)[2],"=red"),xlab="Age",ylab="Length",ylim=c(0,max(x$len)))
		points(len~age,data=x[x$cat==1,],col="red")
		lines(Linf1*(1-exp(-K1*(plotages-t01)))~plotages)
		lines(Linf2*(1-exp(-K2*(plotages-t02)))~plotages,col="red")
 # H4 model
	     Linf1<-nlsout$'model H4'$coefficients[1]
		K1<-nlsout$'model H4'$coefficients[2]
		t01<-nlsout$'model H4'$coefficients[3]
		Linf2<-Linf1
		K2<-K1
		t02<-t01
	    plot(len~age,data=x[x$cat==0,],main=paste("H4 Model ",levels(group)[1],"=black ", levels(group)[2],"=red"),xlab="Age",ylab="Length",ylim=c(0,max(x$len)))
	    points(len~age,data=x[x$cat==1,],col="red")
		lines(Linf1*(1-exp(-K1*(plotages-t01)))~plotages)
		lines(Linf2*(1-exp(-K2*(plotages-t02)))~plotages,col="red")
      }
   if (plottype==2){
 		par(mfrow=c(3,2))
		plotages<-seq(min(x$age),max(x$age)+1,1)
          # Ho model
	     Linf1<-nlsout$'model Ho'$coefficients[1]
		K1<-nlsout$'model Ho'$coefficients[3]
		t01<-nlsout$'model Ho'$coefficients[5]
		Linf2<-Linf1+nlsout$'model Ho'$coefficients[2]
		K2<-K1+nlsout$'model Ho'$coefficients[4]
		t02<-t01+nlsout$'model Ho'$coefficients[6]
          obs1<-x[x$cat==0,]
          obs2<-x[x$cat==1,]
          res1<-obs1$len-Linf1*(1-exp(-K1*(x[x$cat==0,2]-t01)))
          res2<-obs2$len-Linf2*(1-exp(-K2*(x[x$cat==1,2]-t02)))
         	plot(res1~age,data=x[x$cat==0,],main=paste("Ho Model ",levels(group)[1],"=black ", levels(group)[2],"=red") ,xlab="Age",ylab="Residual",ylim=c(-max(abs(res1),abs(res2)),max(abs(res1),abs(res2))))
		points(res2~age,data=x[x$cat==1,],col="red")
          abline(h=0)
   # H1 model
	     Linf1<-nlsout$'model H1'$coefficients[1]
		K1<-nlsout$'model H1'$coefficients[2]
		t01<-nlsout$'model H1'$coefficients[4]
		Linf2<-Linf1
		K2<-K1+nlsout$'model H1'$coefficients[3]
		t02<-t01+nlsout$'model H1'$coefficients[5]
		  obs1<-x[x$cat==0,]
          obs2<-x[x$cat==1,]
          res1<-obs1$len-Linf1*(1-exp(-K1*(x[x$cat==0,2]-t01)))
          res2<-obs2$len-Linf2*(1-exp(-K2*(x[x$cat==1,2]-t02)))
         	plot(res1~age,data=x[x$cat==0,],main=paste("H1 Model ",levels(group)[1],"=black ", levels(group)[2],"=red") ,xlab="Age",ylab="Residual",ylim=c(-max(abs(res1),abs(res2)),max(abs(res1),abs(res2))))
		points(res2~age,data=x[x$cat==1,],col="red")
          abline(h=0)
    # H2 model
	     Linf1<-nlsout$'model H2'$coefficients[1]
		K1<-nlsout$'model H2'$coefficients[3]
		t01<-nlsout$'model H2'$coefficients[4]
		Linf2<-Linf1+nlsout$'model H2'$coefficients[2]
		K2<-K1
		t02<-t01+nlsout$'model H2'$coefficients[5]
		  obs1<-x[x$cat==0,]
          obs2<-x[x$cat==1,]
          res1<-obs1$len-Linf1*(1-exp(-K1*(x[x$cat==0,2]-t01)))
          res2<-obs2$len-Linf2*(1-exp(-K2*(x[x$cat==1,2]-t02)))
       	plot(res1~age,data=x[x$cat==0,],main=paste("H2 Model ",levels(group)[1],"=black ", levels(group)[2],"=red") ,xlab="Age",ylab="Residual",ylim=c(-max(abs(res1),abs(res2)),max(abs(res1),abs(res2))))
		points(res2~age,data=x[x$cat==1,],col="red")
          abline(h=0)
 # H3 model
	     Linf1<-nlsout$'model H3'$coefficients[1]
		K1<-nlsout$'model H3'$coefficients[3]
		t01<-nlsout$'model H3'$coefficients[5]
		Linf2<-Linf1+nlsout$'model H3'$coefficients[2]
		K2<-K1+nlsout$'model H3'$coefficients[4]
		t02<-t01
		  obs1<-x[x$cat==0,]
          obs2<-x[x$cat==1,]
          res1<-obs1$len-Linf1*(1-exp(-K1*(x[x$cat==0,2]-t01)))
          res2<-obs2$len-Linf2*(1-exp(-K2*(x[x$cat==1,2]-t02)))
        	plot(res1~age,data=x[x$cat==0,],main=paste("H3 Model ",levels(group)[1],"=black ", levels(group)[2],"=red") ,xlab="Age",ylab="Residual",ylim=c(-max(abs(res1),abs(res2)),max(abs(res1),abs(res2))))
		points(res2~age,data=x[x$cat==1,],col="red")
          abline(h=0)
 # H4 model
	     Linf1<-nlsout$'model H4'$coefficients[1]
		K1<-nlsout$'model H4'$coefficients[2]
		t01<-nlsout$'model H4'$coefficients[3]
		Linf2<-Linf1
		K2<-K1
		t02<-t01
		obs1<-x[x$cat==0,]
            obs2<-x[x$cat==1,]
            res1<-obs1$len-Linf1*(1-exp(-K1*(x[x$cat==0,2]-t01)))
            res2<-obs2$len-Linf2*(1-exp(-K2*(x[x$cat==1,2]-t02)))
        	plot(res1~age,data=x[x$cat==0,],main=paste("H4 Model ",levels(group)[1],"=black ", levels(group)[2],"=red") ,xlab="Age",ylab="Residual",ylim=c(-max(abs(res1),abs(res2)),max(abs(res1),abs(res2))))
		points(res2~age,data=x[x$cat==1,],col="red")
          abline(h=0)
      }
   }
  return(nlsout)
}

