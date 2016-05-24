prince <-
function(g,o,top=25,imputeknn=F,center=T,permute=F){

      if(is.matrix(g)!=T){stop("g is not a matrix")}
      if(is.data.frame(o)!=T){stop("o is not a data.frame")}
      classes<-unlist(lapply(unclass(o),class))
      if(any(classes=="character")){stop("o contains characters")}
      nrlevels<-unlist(lapply(unclass(o),function(x)length(levels(x))))
      if(any(nrlevels==1)){stop("o contains factors with only one level")}
      if(top>ncol(g)){stop("top is larger than ncol(g)")}
      if(top>nrow(g)){stop("top is larger than nrow(g)")}
      if(!identical(rownames(o),colnames(g))){warning("Colnames of g are not the same as rownames of o")}

      if (imputeknn==T){
         require(impute)
         gimpute<-impute.knn(g)
         g<-gimpute$data
         }
  if(center==T){pr<-prcomp(t(g))}
  if(center==F){pr<-prcomp(t(g),center=F)}
  linp<-matrix(ncol=top,nrow=ncol(o))  
  rownames(linp)<-colnames(o)
  rsquared<-matrix(ncol=top,nrow=ncol(o))
  rownames(rsquared)<-colnames(o)
  for (i in 1:ncol(o)){
  for (j in 1:top){
  fit<-lm(pr$x[,j]~o[,i])
  s<-summary(fit)
  linp[i,j]<-pf(s$fstatistic[1],s$fstatistic[2],s$fstatistic[3],lower.tail=FALSE)
  rsquared[i,j]<-s$r.squared[1]
  }}
  prop<-(pr$sdev[1:top]^2/sum(pr$sdev^2))*100

      if (permute==T){
      gperm<-g
      for (i in 1:nrow(g)){
      gperm[i,]<-sample(g[i,],ncol(g),replace=F)}
      if(center==T){prperm<-prcomp(t(gperm))}
      if(center==F){prperm<-prcomp(t(gperm),center=F)}
      linpperm<-matrix(ncol=top,nrow=ncol(o))
      rownames(linpperm)<-colnames(o)
      rsquaredperm<-matrix(ncol=top,nrow=ncol(o))
      rownames(rsquaredperm)<-colnames(o)
      for (i in 1:ncol(o)){
      for (j in 1:top){
      fitperm<-lm(prperm$x[,j]~o[,i])
      sperm<-summary(fitperm)
      linpperm[i,j]<-pf(sperm$fstatistic[1],sperm$fstatistic[2],sperm$fstatistic[3],lower.tail=FALSE)
      rsquaredperm[i,j]<-sperm$r.squared[1]
      }}
      propperm<-(prperm$sdev[1:top]^2/sum(prperm$sdev^2))*100
      }

  ret<-list(pr=pr,linp=linp,rsquared=rsquared,prop=prop,o=o,prperm=if(permute==T){prperm},linpperm=if(permute==T){linpperm},
            rsquaredperm=if(permute==T){rsquaredperm},propperm=if(permute==T){propperm},imputed=if(imputeknn==T){gimpute$data})

  class(ret)<-"prince"
  return(ret)
  }

