
Cluster.RNASeq=function(data,model,centers=NULL,method=c("EM","DA","SA"),iter.max=30,TMP=NULL){ 
    n=data$Count
    s=data$Normalizer
    t=data$Treatment
    logFC=data$logFC
    n[n<=0]=1e-10
    nG=nrow(n)
   
    if( model=="poisson") nb.disp=rep(1e-10,nG)
    if( model=="nbinom")  nb.disp=data$NB.Dispersion
    if(is.vector(nb.disp)) nb.disp=matrix(rep(nb.disp,ncol(n)),nrow=nG)
    
    if(is.null(TMP)){
      method=method[1]
      if(method=="EM") TMP=rep(1,iter.max)
      if(method %in% c("DA","SA")) TMP=4*.9^(0:iter.max)
     }
    if(length(centers)==1) centers=logFC[sample(nG,centers),]
    C=centers
    P = matrix(1/nrow(C), nrow(n), nrow(C))
        M = matrix(0, nrow(n), nrow(C))
        for (k in 1:nrow(C))
            M[, k] = log(rowSums(n)/rowSums(exp(sweep(s, 2, C[k,t], "+"))))
    k = rep(1,nG)
    for(i in 1:length(TMP)){
        P = cl.mb.est.P(n, s, t, model, P, C, M, nb.disp, a=1/TMP[i],method=method)
        MC = cl.mb.est.MC(n, s, t, model, P, C, M, nb.disp,method=method)
        M = MC$M
        C.new = MC$C
        k.new=apply(P,1,which.max)
        dif =mean(k!=k.new)
        if (i>20 & dif < 0.001) break
        C = C.new
        k = k.new
   }
  print(paste("--->>>> iteration stops after",i,"steps"))
  return(list(probability=P,centers=C,cluster=k))
}

