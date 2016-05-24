FPDC<- function(data=NULL,  k=2, nf=2, nu=2){
  
#   %%%%%%%%%%%%%%Performs Factor PD-Clustering%%%%%%%%%%%%%%%%%%
#     %%%%%%%%INPUT%%%%%%%%%
#   %data=data matrix  nxp
#   %k= number of clusters
#   %nf= number of factors
#   %%%%%%%OUTPUT%%%%%%%%%%
#   %la=cluster labels
#   %iter= number of iterations
#   %dist=distance matrix
#   %p=probability matrix,
#   %expl=explained variability
#   %c= center matrix
#   %JDF= Joint distance function

  data=as.matrix(data)
  
  n=nrow(data)
  J=ncol(data)
  
  #inizialization of convergence measure
  iter=0 #count the number of iterations
  nco=k #check if the number of clusters change
  JDF=(1/.Machine$double.eps)-1#  Join Distance Function
  JDFo=(1/.Machine$double.eps)
  JDFv=0 #to verify the convrgence
  expl=matrix(0,1,100)
  #Step0: inizialization of centers
 # km=kmeans(data, nc)
  #c=km$centers
  c=matrix(0,k,J)
  c[1,]=min(data)
  c[2,]=max(data)
  if(k>2){for(i in 3:k){
    c[i,]=runif(J)
  }}
  
  #Step1-5 iterative cicle, stop if JDF stop decreasing or at 100 iterations
  while((JDF-JDFo)<0.01 & iter<100){
    iter=iter+1
  JDFo=JDF
    
    
    #Step 1 Computation of the (n*nc)x p distance array 
    #distances are also computed as norm
    dis=matrix(0,n,k)
    ddt=matrix(0,n,(k*J))
    dd=matrix(0,n,J)
    for( j in 1:k){
    for( i in 1:n){
      v=t(as.vector(abs(data[i,]-c[j,])))
    dd[i,]=v
    }
    dis[,j]=sqrt(rowSums((data-matrix(c[j,],n,J,1))^2))#needed in the computation of probabilities
    ddt[,((1+J*(j-1)):(J*(j-1)+J))]=dd
    }
    
    #Computation of probabilities
    p=matrix(0,n,k)
    for( i in 1:k){
      t2=matrix(dis[,i],n,k)
      t=t2/dis
      p[,i]=rowSums(t)
    }
    p=1/p
    #Step 2 Tucker3
    
    nfc=k
    if( k>2){
      nfc=k-1# if there are more than two clusters the number of factor for clusters is nc-1;
    }
    tuk3=T3funcrep(ddt, n, J, k, nu, nf,nfc, start=0,conv=0.1)
    expl[iter]=tuk3$fp
    tuk2=tuk3$B
    
    #Step 3 projection of data in the new subspace
    tuk=data%*%tuk2 
    
    #Step 4 PDclustering
   pd=PDclust(tuk,k)
    c=pd$centers
    
#     %check for NaN
#     if sum(sum(isnan(c)))>0
#     idx= find(sum(isnan(c),2));
#     c(idx,:)=[]; %#ok<FNDSB>
#     nc=size(c,1);
#     
#     end
    
    
    #Step 5 update of centers
    c=c%*%t(tuk2)
    
    #labeling according to the higest probability
    l=max.col(p)#pd$p
    
#     %check: if two centers are equal K=K-1
#     if nc>2
#     for k=1:(nc-1)
#     app=sum(repmat(c(k,:),(nc-i),1)==c(i+1:nc,:),2)/ppp;
#     c(find(app),:)=[]; %#ok<FNDSB>
#     nc=size(c,1);
#     end
#     end
    #criteria update

    JDF=dis[,1]%*%p[,1]
    JDFv=cbind(JDFv,JDF)
    cat(paste("Iteration number", iter), fill = TRUE)
  }
  
  JDFv=JDFv[2:(iter+1)]
  expl=expl[1:iter]
  out=list(label=l, centers=c, probability=p, JDF=JDF, JDFIter=JDFv, iter=iter, explained=expl)
  out
}