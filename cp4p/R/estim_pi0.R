
###################################################
# estimate the proportion of true null hypotheses #
###################################################

estim.pi0=function (p, pi0.method = "ALL",  nbins = 20, pz = 0.05){
  #library(pi0)
  #library(qvalue)
 if (pi0.method=="ALL"||pi0.method=="st.spline"||pi0.method=="st.boot"||pi0.method=="jiang"||pi0.method=="histo"||pi0.method=="langaas"||pi0.method=="pounds"||pi0.method=="abh"||pi0.method=="slim"){
    
  ######################################
  #PI0 from LSL method of Benjamini-Hochberg (2000)
  #Reference: On the adaptive control of the false discovery rate in multiple testing with independent statistics.
  pi0.abh=function(p){
    m=length(p)
    sortp=sort(p)
    s=sort(1-sortp)/(1:m)
    m0raw=m
    i=m
    while ((i>1)&&(s[i]<=s[i-1])){i=i-1;}
    if (i>1){ m0raw=1/s[i-1];}
    else{ m0raw=1/s[1];}
    m0=min(floor(1 + m0raw), m)
    pi0=m0/m
    return(pi0)
  }
  ######################################
  #PI0 from Nettleton (2006)
  #Reference: Estimating the number of true null hypotheses from a histogram of p values.
  pi0.histo=function(p, nbin){
    bin=c(-0.1, (1:nbins)/nbin)
    bin.counts=tabulate(cut(p,bin))
    tail.means=rev(cumsum(rev(bin.counts))/(1:nbin))
    index=which(tail.means >= bin.counts)[1]
    tail.means[index]/tail.means[1]
  }
  ######################################
  #PI0 from Jiang and Doerge (2008)
  #Reference: Estimating the proportion of true null hypotheses for multiple comparisons.
  pi0.jiang=function(p, nbin){
    m=length(p) 
    t=seq(0,1,length=nbin+1)
    NB=rep(0,nbin)
    NBaverage=rep(0,nbin)
    NS=rep(0,nbin)
    pi=rep(0,nbin)
    for(i in 1:nbin){
      NB[i]=length(p[p>=t[i]])
      NBaverage[i]=NB[i]/(nbin-(i-1))
      NS[i]=length(p[p>=t[i]]) - length(p[p>=t[i+1]])
      pi[i]=NB[i]/(1-t[i])/m
    }
    i=min(which(NS <= NBaverage))
    pi0=min(1, mean(pi[(i-1):nbin]))
    return(pi0)
  }
  ######################################
  #PI0 from SLIM (2011)
  #Adapted from https://github.com/al2na/methylKit/blob/master/R/diffMeth.R
  #Copyright by Tsai Lab of UGA, US, and Hong-Qiang Wang, IIM, CAS, China
  #Reference: SLIM: A Sliding Linear Model for Estimating the Proportion of True Null Hypotheses in Datasets With Dependence Structures
  #inputs: 
  # rawp:p-values
  # STA:lambda1
  # Divi:number of segments (n)
  # Pz:maximum of p-values of alternative tests (pmax)
  # B:number of quantile points
  #outputs: 
  # pi0: estimated value of pi0
  ###############
  pi0.slim=function(rawp,STA=0.1,Divi=10,Pz=0.05,B=100){
    pi0s_est_COM=NULL;
    P_pi1_mtx=NULL;
    pi1_act_mtx=NULL;
    PI0=NULL;
    ########f1 function
    f1<-function(cutoff,rawp){sum(rawp<cutoff)/length(rawp)};
    ########estimation
    pi0_mtx=NULL;
    itv=(1-STA)/Divi;
    for (i in 1:Divi){
      cutoff=STA+(i/Divi)*(1-STA);
      lambda=seq(cutoff-itv,cutoff,itv/10);
      gamma_mtx=sapply(lambda,f1,rawp=rawp);
      LModel=lm(gamma_mtx~lambda);
      pi0_mtx=c(pi0_mtx,coefficients(LModel)[2]);
    }
    ########QValuesfun function
    QValuesfun=function(rawp,pi0){
      order_rawp=sort(rawp);
      qvalues=pi0*length(order_rawp)*order_rawp/c(1:length(order_rawp));
      temp=cummin(qvalues[seq(length(qvalues),1,-1)])
      qvalues=temp[seq(length(temp),1,-1)];
      qvalues=qvalues[order(order(rawp))]
    }
    ########searching
    maxFDR_mtx=NULL;
    quapoint_mtx=seq(0.01,0.99,1/B);
    for (k in 1:length(quapoint_mtx)){
      qua_point=quapoint_mtx[k];
      pi0_combLR=min(quantile(pi0_mtx,qua_point),1);
      pi0_est=pi0_combLR;
      ########Calculate independent index of raw p vlaues
      PI0=rbind(PI0,pi0_mtx);
      pi0s_est_COM=c(pi0s_est_COM,pi0_est);
      ########Condition 1
      P_pi1=sort(rawp)[max(length(rawp)*(1-pi0_est),1)];
      P_pi1_mtx=c(P_pi1_mtx,P_pi1);
      pi0=pi0_est;
      maxFDR=Pz*pi0/(1-(1-Pz)*pi0);
      maxFDR_mtx=c(maxFDR_mtx,maxFDR);
      qvalues_combLR=QValuesfun(rawp,pi0);
      qvalue_cf=maxFDR;
      selected=which(qvalues_combLR<qvalue_cf);
      Sel_qvalues_combLR=selected;
      
      pi1_act_mtx=c(pi1_act_mtx,length(Sel_qvalues_combLR)/length(rawp));
    }
    ####doing judging
    ##by max FDR
    pi1s_est_COM=1-pi0s_est_COM;
    Diff=sum(rawp<=Pz)/length(rawp)-pi1_act_mtx;
    ###
    loc=which.min(abs(Diff));
    Diff.loc=Diff[loc];
    pi0_Est=min(1,pi0s_est_COM[loc]);
    
    return(list(pi0=pi0_Est));
  }
  ######################################
  #OUTCOMES
  
  if (pi0.method=="ALL"){
    if (!inherits(try(qvalue(p,pi0.method="smoother")$pi0,FALSE), "try-error")==TRUE){
      Storey.S=qvalue(p,pi0.method="smoother");
      pi0.Storey.Spline=Storey.S$pi0;
    }else{warning("Warning: st.spline does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.Storey.Spline=1;}
    if (!inherits(try(qvalue(p,pi0.method="bootstrap")$pi0,FALSE), "try-error")==TRUE){
      Storey.B=qvalue(p,pi0.method="bootstrap");
      pi0.Storey.Boot=Storey.B$pi0
    }else{warning("Warning: st.boot does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.Storey.Boot=1;}
    if (!inherits(try(convest(p),FALSE), "try-error")==TRUE){
      Langaas=convest(p);
      pi0.Langaas=Langaas[1];
    }else{warning("Warning: langaas does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.Langaas=1;}
    if (!inherits(try(pi0.histo(p,nbins),FALSE), "try-error")==TRUE){
      pi0.Histo=pi0.histo(p,nbins);
      if (pi0.Histo>1){warning("Warning: histo does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.Histo=1;}
    }else{warning("Warning: histo does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.Histo=1;}
    if (!inherits(try(pi0.jiang(p,nbins),FALSE), "try-error")==TRUE){
      pi0.Jiang=pi0.jiang(p,nbins);
    }else{warning("Warning: jiang does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.Jiang=1;}
    if (!inherits(try(min(1,2*mean(p)),FALSE), "try-error")==TRUE){
      pi0.Pounds=min(1,2*mean(p));
    }else{warning("Warning: pounds does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.Pounds=1;}
    if (!inherits(try(pi0.abh(p),FALSE), "try-error")==TRUE){
      pi0.ABH=pi0.abh(p);
    }else{warning("Warning: abh does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.ABH=1;}
    if (!inherits(try(pi0.slim(p, Pz=pz),FALSE), "try-error")==TRUE){
      slim=pi0.slim(p, Pz=pz);
      pi0.SLIM=slim[[1]];
    }else{warning("Warning: slim does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.SLIM=1;}
    res=data.frame(pi0.Storey.Spline, pi0.Storey.Boot, pi0.Jiang, pi0.Histo, pi0.Langaas, pi0.Pounds, pi0.ABH, pi0.SLIM);
    return(list(pi0.est=res)); 
  }
  if (pi0.method=="st.spline"){
    if (!inherits(try(qvalue(p,pi0.method="smoother")$pi0,FALSE), "try-error")==TRUE){
      Storey.S=qvalue(p,pi0.method="smoother");
      pi0.Storey.Spline=Storey.S$pi0;
    }else{warning("Warning: st.spline does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.Storey.Spline=1;}
    return(list(pi0.Storey.Spline=pi0.Storey.Spline)) ; 
  }
  if (pi0.method=="st.boot"){
    if (!inherits(try(qvalue(p,pi0.method="bootstrap")$pi0,FALSE), "try-error")==TRUE){
      Storey.B=qvalue(p,pi0.method="bootstrap");
      pi0.Storey.Boot=Storey.B$pi0
    }else{warning("Warning: st.boot does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.Storey.Boot=1;}
    return(list(pi0.Storey.Boot=pi0.Storey.Boot)) ; 
  }
  if (pi0.method=="langaas"){
    if (!inherits(try(convest(p),FALSE), "try-error")==TRUE){
      Langaas=convest(p);
      pi0.Langaas=Langaas[1];
    }else{warning("Warning: langaas does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.Langaas=1;}
    return(list(pi0.Langaas=pi0.Langaas)) ; 
  }
  if (pi0.method=="histo"){
    if (!inherits(try(pi0.histo(p,nbins),FALSE), "try-error")==TRUE){
      pi0.Histo=pi0.histo(p,nbins);
      if (pi0.Histo>1){warning("Warning: histo does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.Histo=1;}
    }else{warning("Warning: histo does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.Histo=1;}
    return(list(pi0.Histo=pi0.Histo)) ; 
  }
  if (pi0.method=="jiang"){
    if (!inherits(try(pi0.jiang(p,nbins),FALSE), "try-error")==TRUE){
      pi0.Jiang=pi0.jiang(p,nbins);
    }else{warning("Warning: jiang does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.Jiang=1;}
    return(list(pi0.Jiang=pi0.Jiang)) ; 
  }
  if (pi0.method=="pounds"){
    if (!inherits(try(min(1,2*mean(p)),FALSE), "try-error")==TRUE){
      pi0.Pounds=min(1,2*mean(p));
    }else{warning("Warning: pounds does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.Pounds=1;}
    return(list(pi0.Pounds=pi0.Pounds)) ; 
  }
  if (pi0.method=="abh"){
    if (!inherits(try(pi0.abh(p),FALSE), "try-error")==TRUE){
      pi0.ABH=pi0.abh(p);
    }else{warning("Warning: abh does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.ABH=1;}
    return(list(pi0.ABH=pi0.ABH)) ; 
  }
  if (pi0.method=="slim"){
    if (!inherits(try(pi0.slim(p, Pz=pz),FALSE), "try-error")==TRUE){
      slim=pi0.slim(p, Pz=pz);
      pi0.SLIM=slim[[1]];
    }else{warning("Warning: slim does not work correctly (pi0 is thus fixed to 1)!!\n\n");pi0.SLIM=1;}
    return(list(pi0.SLIM=pi0.SLIM)) ; 
  }
 }else{warning("\n Error in input pi0.method:\n Please write the name of an estimation method among st.spline, st.boot, jiang, histo, langaas, pounds, abh, slim or ALL.\n");}
}
