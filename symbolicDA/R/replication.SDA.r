


.clust<-function(d,u,method){
      if(method=="pam")
      {
            cl<-pam(d,u,diss=TRUE)$clustering
      }
      else if(method=="diana"){
               cla<-cutree(as.hclust(diana(d)),k=u)
      }
      else if(method=="sclust" || method=="dclust"){
               cl<-DClust(d,u)
      }
      else
      {
         cl<-cutree(hclust(d,method=method),u)
      }
}

replication.SDA<-function(table.Symbolic, u=2, method="SClust", S=10,fixedAsample=NULL,...)
{
   x<-table.Symbolic
   method<-casefold(method)
   nrObjects<-nrow(x$individuals)
   if(!is.null(u))
   {
      if(u<2 || u>(nrObjects-1))
      stop("number of classes must be between 2 and ",(nrObjects-1))
   }
   half<-as.integer(nrObjects/2)
   a_A<-array(0,c(S,half))
   a_B<-array(0,c(S,nrObjects-half))
   a_centroA<-array(0,c(S,u))
   a_clA<-array(0,c(S,half))
   a_clB<-array(0,c(S,nrObjects-half))
   a_clBB<-array(0,c(S,nrObjects-half))
   a_cRand<-array(0,S)
   for(s in 1:S)
   {
      sampleA<-sample(1:nrObjects,half)
      if (!is.null(fixedAsample))
        {
        if(is.null(nrow(fixedAsample)) || nrow(fixedAsample)==1)
        {		
          SampleA<-fixedAsample
        }
        else
        {		
          SampleA<-fixedAsample[s,]
        }
      }
      sampleB<-(1:nrObjects)[-sampleA]
      #print(sampleA)
      #print(sampleB)
      #print("debug : 8")
      if (method=="sclust")
      {
          d<-.SDist(x,objectSelection=1:nrObjects,variableSelection=1:nrow(x$variables))
      #print("debug : 8a")
          dA<-.SDist(x,objectSelection=sampleA,variableSelection=1:nrow(x$variables))
      #print("debug : 8b")
          dB<-.SDist(x,objectSelection=sampleB,variableSelection=1:nrow(x$variables))
      #print("debug : 8c")
      }
      else{
          #t<-sample(1:33,33)
          #d<-dist(t)
          d<-dist.SDA(x,...)
          dA<-as.matrix(d)[sampleA,sampleA]
          dimnames(dA)[[1]]<-1:length(sampleA)
          dimnames(dA)[[2]]<-1:length(sampleA)
          dB<-as.matrix(d)[sampleB,sampleB]
          dimnames(dB)[[1]]<-1:length(sampleB)
          dimnames(dB)[[2]]<-1:length(sampleB)
          dA<-as.dist(dA)
          dB<-as.dist(dB)
          #print(dA)
          #print(dB)
      }
      a_A[s,]<-sampleA
      a_B[s,]<-sampleB
      #print("debug 8.5")
      clA<-.clust(dA,u,method)
      clB<-.clust(dB,u,method)
      #print(clA)
      #print(clB)
      for(i in 1:u){
        a_centroA[s,i]<-.medoid(as.matrix(dA),clA,i)
      }
      #print("debug : 9")
      d<-as.matrix(d)
      for(i in 1:length(sampleB))
      {
       t<-rep(0,u)
       for(j in 1:u){
        t[j]<-d[sampleA[a_centroA[s,j]],sampleB[i]]
       }
       a_clBB[s,i]<-which.min(t)
      }
      #print("debug : 10")
      a_clA[s,]<-clA
      a_clB[s,]<-clB
      ca<-classAgreement(table(clB,a_clBB[s,]),match.names=FALSE)
      #print(ca$crand)
      a_cRand[s]<-ca$crand
   }
   resulMedoids=a_centroA
   resulcRand<-mean(a_cRand)
   resul<-list(A=a_A,B=a_B,medoids=resulMedoids,clusteringA=a_clA,clusteringB=a_clB,clusteringBB=a_clBB,cRand=resulcRand)
   resul
}

