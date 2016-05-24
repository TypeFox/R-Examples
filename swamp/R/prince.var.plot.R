prince.var.plot <-
function(g,show.top=dim(g)[2],imputeknn=F,center=T,npermute=10){
      if(is.matrix(g)!=T){stop("g is not a matrix")}
            if(show.top>ncol(g)){stop("top is larger than ncol(g)")}
            if(show.top>nrow(g)){stop("top is larger than nrow(g)")}

      if (imputeknn==T){
         require(impute)             
         gimpute<-impute.knn(g)
         g<-gimpute$data
         }
      if(center==T){pr<-prcomp(t(g))}
      if(center==F){pr<-prcomp(t(g),center=F)}
      prop<-(pr$sdev[1:dim(g)[2]]^2/sum(pr$sdev^2))*100
      names(prop)<-paste("PC",1:ncol(g))
        permmat<-matrix(nrow=ncol(g),ncol=npermute,dimnames=list(paste("PC",1:ncol(g)),paste("Permutation",1:npermute)))
        for (j in 1:npermute){
        gperm<-g
        for (i in 1:nrow(g)){
        gperm[i,]<-sample(g[i,],ncol(g),replace=F)}
        if(center==T){prperm<-prcomp(t(gperm))}
        if(center==F){prperm<-prcomp(t(gperm),center=F)}
        permmat[,j]<-(prperm$sdev[1:dim(gperm)[2]]^2/sum(prperm$sdev^2))*100
        print(paste("Perm =",j))
        }
      medperm<-apply(permmat,1,median)
      plot(1:show.top,prop[1:show.top],xlab="Principal Components",ylab="Variation in %")
      points(1:show.top,medperm[1:show.top],col=2)
      legend(show.top*0.7,max(prop),legend=c("Observed Data","Permuted Data"),pch=1,col=c(1,2))
      
   return(list(real.variation=prop,permuted.variation=permmat))
   }

