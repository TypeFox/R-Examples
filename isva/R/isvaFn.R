`isvaFn` <-
function(data.m,pheno.v,ncomp=NULL){

  lm.o <- lm(t(data.m) ~ pheno.v);
  res.m <- t(lm.o$res);
  model <- model.matrix(~1+pheno.v);
  if(is.null(ncomp)){
    rmt.o <-  EstDimRMT(res.m)
    ncomp <- rmt.o$dim;
    print(paste("Number of candidate ISVs = ",ncomp,sep=""));
  }
  else {
    print("no need to estimate dimensionality");
  }

  ### perform ICA on residual matrix
  fICA.o <- fastICA(res.m,n.comp=ncomp);

  ### now construct ISV
  tmp.m <- t(fICA.o$A);
  isv.m <- tmp.m;
  sd <- 1/sqrt(ncol(data.m)-3);
  for(k in 1:ncol(tmp.m)){
   cor.v <- as.vector(cor(t(data.m),tmp.m[, k]))
   z.v <- 0.5*log((1+cor.v)/(1-cor.v));
   pv.v <- 2*pnorm(abs(z.v),0,sd,lower.tail=FALSE)
   tmp.s <- sort(pv.v,decreasing=FALSE,index.return=TRUE);
   qv.o <- qvalue(pv.v);
   nsig <- length(which(qv.o$qvalues<0.05));
   if( nsig < 500 ){
     nsig <- 500;
   }
   red.m <- data.m[tmp.s$ix[1:nsig],];
   fICA.o <- fastICA(red.m,n.comp=ncomp);
   cor.v <- abs(cor(tmp.m[,k],t(fICA.o$A)));
   kmax <- which.max(cor.v);
   isv.m[,k] <- t(fICA.o$A)[,kmax];
   print(paste("Built ISV ",k,sep=""));   
  }
  return(list(n.isv=ncomp,isv=isv.m));
}

