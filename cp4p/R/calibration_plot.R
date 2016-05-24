
############################################################
# plot to verify the assumptions of FDR control procedures #
############################################################

calibration.plot=function(p, pi0.method="pounds",  nbins = 20, pz=0.05){
  #library(MESS);
  Fr=ecdf(1-p);
  abs=seq(0,1,by=0.001);
  fc=Fr(abs);
  AUC=NULL;
  AUC2=NULL;
  plot(abs,fc,ty="l",xlab="1-p.value",ylab="Cumulative Distribution Function of 1-p.value",lwd=2,ylim=c(0,1.05));
  if (is.numeric(pi0.method)==TRUE){
    if (pi0.method<=1 && pi0.method>=0){
      pi0=pi0.method;
      title(main=paste("Calibration Plot - pi0 =",pi0.method))

      #straight line y=pi0*x
      dr=as.numeric(pi0)*abs;
      dif=abs-fc;
      dif[dif<0]=0;
      dif[dif==0]=0.0001;
      f=fc[dif>0]-dr[dif>0];
      f[f<0]=0;
      f=f+dr[dif>0];
      #Grey area
      x.bons=c(abs[dif>0],abs[dif>0][1]);
      y.bons=c(f,f[1]);
      polygon(x = x.bons, y = y.bons, col = "gray93",border="white");
    
      dif=fc-dr;
      lfc=length(fc);
      i=lfc;
      while (sign(dif[i-1])==sign(dif[i])){i=i-1;}
      #Green area
      x.bons=c((i:lfc)/lfc,1,i/lfc,i/lfc);
      y.bons=c(fc[i:lfc],pi0,pi0*i/lfc,fc[i]);
      polygon(x = x.bons, y = y.bons, col = "darkseagreen1",border="darkseagreen1");  
      #Compute DA protein concentration
      if (i<lfc){auc1=auc(abs[i:lfc],dif[i:lfc]);}else{auc1=0;}
      if ((1-pi0)>2*auc1){
        AUC=((1-pi0)/2-auc1)/((1-pi0)/2);
      }else{AUC=0;}
      text(0,0.9,paste("DA protein concentration: ",floor(AUC*1000)/10,"%"),col="green4",pos=4);
      #Compute divergence to uniformity assumption
      dif2=dif;
      dif2[dif2<0]=0;
      fc2=dr+dif2;
      AUC2=auc(abs[1:(i-1)],dif2[1:(i-1)]);
      text(0,0.85,paste("Uniformity underestimation: ",floor(AUC2*100000)/1000),col="red",pos=4);
      #Red area
      x.bons=c((1:(i-1))/1000,0);
      y.bons=c(fc2[1:(i-1)],0);
      polygon(x = x.bons, y = y.bons, col = "red",border="white");
    
      lines(abs,dr,lwd=2,col="blue",lty=1);
      lines(abs,fc,lwd=2);
      text(0,0.95,paste("Non-DA protein proportion: ",floor(pi0*1000)/10, "%"),col="blue",pos=4);
    }else{warning("\n Error in input pi0.method:\n Please write a numeric value between 0 and 1 or the name of an estimation method among st.spline, st.boot, jiang, histo, langaas, pounds, abh, slim or ALL.\n");}
  }
  if (is.numeric(pi0.method)==FALSE){
    if (pi0.method=="ALL"||pi0.method=="st.spline"||pi0.method=="st.boot"||pi0.method=="jiang"||pi0.method=="histo"||pi0.method=="langaas"||pi0.method=="pounds"||pi0.method=="abh"||pi0.method=="slim"){
      r=estim.pi0(p=p, pi0.method = pi0.method,  nbins = nbins, pz=pz);
      pi0=r$pi0;
      if (pi0.method=="ALL"){
        for (i in 1:length(pi0)){
          if (i==1){lines(abs,as.numeric(pi0[i])*abs,ty="l",col="chartreuse1", lty=i,lwd=2);name="st.spline";}
          if (i==2){lines(abs,as.numeric(pi0[i])*abs,ty="l",col="chartreuse3", lty=i,lwd=2);name=c(name,"st.boot");}
          if (i==3){lines(abs,as.numeric(pi0[i])*abs,ty="l",col="darkgreen", lty=i,lwd=2);name=c(name,"jiang");}
          if (i==4){lines(abs,as.numeric(pi0[i])*abs,ty="l",col="darkorange", lty=i,lwd=2);name=c(name,"histo");}
          if (i==5){lines(abs,as.numeric(pi0[i])*abs,ty="l",col="darkorange3", lty=i,lwd=2);name=c(name,"langaas");}
          if (i==6){lines(abs,as.numeric(pi0[i])*abs,ty="l",col="red", lty=i,lwd=2);name=c(name,"pounds");}
          if (i==7){lines(abs,as.numeric(pi0[i])*abs,ty="l",col="deepskyblue", lty=i,lwd=2);name=c(name,"abh");}
          if (i==8){lines(abs,as.numeric(pi0[i])*abs,ty="l",col="blue", lty=i,lwd=2);name=c(name,"slim");}
        }
        title(main="Calibration Plot - All methods");
        legend("topleft",name,col=c("chartreuse1","chartreuse3","darkgreen","darkorange","darkorange3","red","deepskyblue","blue"),lty=1:8,lwd=rep(2,8),title="pi0.method");
      }else{
        title(main=paste("Calibration Plot -",pi0.method,"method"));
        
        #straight line y=pi0*x
        dr=as.numeric(pi0)*abs;
        dif=abs-fc;
        dif[dif<0]=0;
        dif[dif==0]=0.0001;
        f=fc[dif>0]-dr[dif>0];
        f[f<0]=0;
        f=f+dr[dif>0];
        #Grey area
        x.bons=c(abs[dif>0],abs[dif>0][1]);
        y.bons=c(f,f[1]);
        polygon(x = x.bons, y = y.bons, col = "gray93",border="white");
      
        dif=fc-dr;
        lfc=length(fc);
        i=lfc;
        while (sign(dif[i-1])==sign(dif[i])){i=i-1;}
        #Green area
        x.bons=c((i:lfc)/lfc,1,i/lfc,i/lfc);
        y.bons=c(fc[i:lfc],pi0,pi0*i/lfc,fc[i]);
        polygon(x = x.bons, y = y.bons, col = "darkseagreen1",border="darkseagreen1");  
        #Compute DA protein concentration
        if (i<lfc){auc1=auc(abs[i:lfc],dif[i:lfc]);}else{auc1=0;}
        if ((1-pi0)>2*auc1){
          AUC=((1-pi0)/2-auc1)/((1-pi0)/2);
        }else{AUC=0;}
        text(0,0.9,paste("DA protein concentration: ",floor(AUC*1000)/10,"%"),col="green4",pos=4);
        #Compute divergence to uniformity assumption
        dif2=dif;
        dif2[dif2<0]=0;
        fc2=dr+dif2;
        AUC2=auc(abs[1:(i-1)],dif2[1:(i-1)]);
        text(0,0.85,paste("Uniformity underestimation: ",floor(AUC2*100000)/1000),col="red",pos=4);
        #Red area
        x.bons=c((1:(i-1))/1000,0);
        y.bons=c(fc2[1:(i-1)],0);
        polygon(x = x.bons, y = y.bons, col = "red",border="white");
      
        lines(abs,dr,lwd=2,col="blue",lty=1);
        lines(abs,fc,lwd=2);
        text(0,0.95,paste("Non-DA protein proportion: ",floor(pi0*1000)/10, "%"),col="blue",pos=4);
      }
    }else{warning("\n Error in input pi0.method:\n Please write a numeric value between 0 and 1 or the name of an estimation method among st.spline, st.boot, jiang, histo, langaas, pounds, abh, slim or ALL.\n");}
  }
  return(list(pi0=pi0,h1.concentration=AUC,unif.under=AUC2*100));
}
