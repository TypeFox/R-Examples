############################
# compute ajusted p-values #
############################

adjust.p=function(p, pi0.method=1, alpha=0.05, nbins = 20, pz=0.05){
  #library(multtest)
  if ((pi0.method=="bky")==0){
    if (is.numeric(pi0.method)==TRUE){
      if (pi0.method<=1 && pi0.method>=0){
          pi0=pi0.method;
          if (pi0.method==1){cat("Procedure of Benjamini-Hochberg is used. pi0 is fixed to 1.\n");}
      }else{stop("\n Error in input pi0.method:\n Please write a numeric value between 0 and 1 or the name of an estimation method among st.spline, st.boot, jiang, histo, langaas, pounds, abh or slim.\n");}
    }
    if (is.numeric(pi0.method)==FALSE){
      if (pi0.method=="st.spline"||pi0.method=="st.boot"||pi0.method=="jiang"||pi0.method=="histo"||pi0.method=="langaas"||pi0.method=="pounds"||pi0.method=="abh"||pi0.method=="slim"){
          r=estim.pi0(p, pi0.method = pi0.method,  nbins = nbins, pz=pz);
          pi0=r$pi0
          cat(pi0.method,"method is used. pi0 is estimated to",as.numeric(pi0),".\n");
      }else{stop("\n Error in input pi0.method:\n Please write a numeric value between 0 and 1 or the name of an estimation method among st.spline, st.boot, jiang, histo, langaas, pounds, abh or slim.\n");}
    }
    qa=mt.rawp2adjp(p, proc = "BH");
    adjp=data.frame(qa$adjp[order(qa$index),1],qa$adjp[order(qa$index),2]*as.numeric(pi0));
    colnames(adjp)=c("rawp","adjusted.p");
    return(list(pi0=pi0,adjp=adjp));
  }
  else{
    qa=mt.rawp2adjp(p, proc = "TSBH", alpha=alpha);
    adjp=data.frame(qa$adjp[order(qa$index),1:2]);
    pi0=as.numeric(qa$h0.TSBH)/length(p);
    cat("Procedure of BKY with a FDR of",alpha,"is used. pi0 is estimated to",pi0,".\n");
    colnames(adjp)=c("rawp","adjusted.p");
    return(list(pi0=pi0,adjp=adjp));
  }
}
