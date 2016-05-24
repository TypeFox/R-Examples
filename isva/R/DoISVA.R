`DoISVA` <-
function(data.m,pheno.v,cf.m=NULL,factor.log,pvthCF=0.01,th=0.05,ncomp=NULL){

 ### Main ISVA function
 isva.o <- isvaFn(data.m,pheno.v,ncomp);

 if(is.null(cf.m)==FALSE){

 ### study pattern of correlation of ISVA components to POI and CFs
 tmp.m <- cbind(pheno.v,cf.m);
 treatfactor <- c(FALSE,factor.log);
 pv.m <- matrix(nrow=ncol(isva.o$isv),ncol=1+ncol(cf.m));
 colnames(pv.m) <- c("POI",colnames(cf.m)); ## POI:phenotype of interest
 for(c in 1:ncol(tmp.m)){
  if(treatfactor[c]==FALSE){
   for(sv in 1:ncol(isva.o$isv)){
    lm.o <- lm(isva.o$isv[,sv] ~ as.numeric(tmp.m[,c]));
    pv.m[sv,c] <- summary(lm.o)$coeff[2,4];   
   }
  }
  else {
   for(sv in 1:ncol(isva.o$isv)){
    lm.o <- lm(isva.o$isv[,sv] ~ as.factor(tmp.m[,c]));
    pv.m[sv,c] <- pf(summary(lm.o)$fstat[1],summary(lm.o)$fstat[2],summary(lm.o)$fstat[3],lower.tail=FALSE);   
   }
  }
 }

 ### selection of ISVs
 print("Selecting ISVs");
 selisv.idx <- vector();
 for(sv in 1:nrow(pv.m)){

   ncf <- length(which(pv.m[sv,2:ncol(pv.m)]< pvthCF)) ## pvth=0.01
   minpv <- min(pv.m[sv,2:ncol(pv.m)]);
   phpv <- pv.m[sv,1];
   if(ncf > 0){
     if(minpv < phpv){
       selisv.idx <- c(selisv.idx,sv);
     }
   }
 }
 if (length(selisv.idx)==0 ){
   print("No ISVs selected because none correlated with the given confounders. Rerun ISVA with cf.m=NULL option"); stop;
 }
 
 }
 else { ### confounder matrix not given, so select all ISVs
  selisv.idx <- 1:ncol(isva.o$isv);
  pv.m <- NULL;
 }
 print("Running final multivariate regressions with selected ISVs");
 selisv.m <- matrix(isva.o$isv[,selisv.idx],ncol=length(selisv.idx));
 mod <- model.matrix( ~ pheno.v + selisv.m);
 modNULL <- model.matrix( ~ selisv.m);

 df1 <- dim(mod)[2]
 df0 <- dim(modNULL)[2]
 pv.v <- rep(0, nrow(data.m));
 Id <- diag(ncol(data.m))
 resid <- data.m %*% (Id - mod %*% solve(t(mod) %*% mod) %*%
        t(mod))
 rss1 <- rowSums(resid * resid)
 rm(resid)
 residNULL <- data.m %*% (Id - modNULL %*% solve(t(modNULL) %*% modNULL) %*%
        t(modNULL))
 rssNULL <- rowSums(residNULL * residNULL)
 rm(residNULL)
 fstats <- ((rssNULL - rss1)/(df1 - df0))/(rss1/(ncol(data.m) - df1))
 pv.v <- 1 - pf(fstats, df1 = (df1 - df0), df2 = (ncol(data.m) - df1))
 pv.s <- sort(pv.v,decreasing=FALSE,index.return=TRUE);
 qv.v <- qvalue(pv.s$x)$qvalue;
 ntop <- length(which(qv.v < th));
 print(paste("Number of DEGs after ISV adjustment = ",ntop,sep=""));
 if(ntop>0){
  pred.idx <- pv.s$ix[1:ntop];
  ### find t-stats of significant ones
  lm.o <- lm( t(data.m[pred.idx,]) ~ pheno.v + selisv.m );
  tstats.v <- unlist(lapply(summary(lm.o),function(x){ x$coeff[2,3];}));
  lm.m <- cbind(tstats.v,pv.s$x[1:ntop],qv.v[1:ntop]);
  colnames(lm.m) <- c("t-stat","P-value","q-value");
 }
 else {
  pred.idx <- NULL;
  lm.m <- NULL;
 }
 
 return(list(spv=pv.s$x,qv=qv.v,rk=pv.s$ix,ndeg=ntop,deg=pred.idx,lm=lm.m,isv=selisv.m,nsv=length(selisv.idx),pvCF=pv.m,selisv=selisv.idx));
 
} ### END OF FUNCTION

