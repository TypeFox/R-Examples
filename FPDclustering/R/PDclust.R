PDclust<- function(data=NULL,  k=2) {
  # Cluster the data whit pd-clustering Algoritmh
  #
  #
  #%%%%%%INPUT%%%%%%%%%
    #data=input data
  #k=number of cluster
  #
  #
  #%%%%%%OUTPUT%%%%%%%%
 #cnew=cluster's center  
#l=class label
#p=nxk matrix
#probability to belong to each class 
  #JDF  join distance function
  #cont=number of iterations until convergence
  if((!is.double(k))&(!is.integer(k))){stop("The number of clusters (k) must take an integer value.")}
  if(k<2){stop("The number of clusters (k) must be greater than one.");}
  if((k-round(k)!=0)){stop("The number of clusters (k) must be a whole number.");}
  data=as.matrix(data)
  if(!is.double(data)){stop("All elements of data must have type double.");}
  n=nrow(data)
  J=ncol(data)
  cnew=matrix(0,k,J)
  cnew[1,]=min(data)
  cnew[2,]=max(data)
  if(k>2){for(i in 3:k){
    cnew[i,]=runif(J)
  }}
#km=kmeans(data, k)
#cnew=km$centers

ver=100
cont=0

while(ver>0.001 & cont<1000){
  cont=cont+1
  c=cnew
  #STEP 1
  #computation of distance matrix
  dis=matrix(0,n,k)
  for( i in 1:k){
    dis[,i]=sqrt(rowSums((data-matrix(c[i,],n,J,1))^2))
    
  }
  
  #STEP 2
  #computation of centers and probabilities
  p=matrix(0,n,k)
  for( i in 1:k){
    t2=matrix(dis[,i],n,k)
    t=t2/dis
    p[,i]=rowSums(t)
  }
  p=1/p
  u=p^2/dis
  m=u/matrix(colSums(u),n,k,1)
  cnew=t(m)%*%data
  
  #check if centers move
  ver1=matrix(0,1,k)
  for( j in 1:k){
    ver1[j]=sqrt(sum((cnew[j,]-c[j,])^2))
  end
  ver=sum(ver1)}
  
}
  if( cont==1000){
  print('Convergence not reached')}

  #memebership definition according to the maximum probability
  l=max.col(p)
  # computation of JDF
  JDF=sum(mean(dis*p))
  out=list(label=l, centers=cnew, probability=p, JDF=JDF,iter=cont)
  out
}