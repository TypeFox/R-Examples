las_bycgd <-
function(Y,X=Y,lambda1=1,lambda2=1){
  eps<-0.0001; alpha<-1;
  lambda1<-max(apply(Y,1,FUN="mad"))*sqrt(dim(Y)[1])
  lambda2<-max(apply(Y,1,FUN="mad"))*sqrt(dim(Y)[1])*sqrt(log(dim(Y)[2]))
  D<-calculate_direction(Y,X,lambda1,lambda2)#计算x[k+1]=x[k]+alpha*d中的d
     
  while( alpha*max(abs(D))>eps ){#当alpha*max（|d|）>eps时，BCGD算法迭代继续
    alpha<-armijo_rule(Y,X,D,lambda1,lambda2,alpha)#计算x[k+1]=x[k]+alpha*d的alpha
    X<-X+alpha*D#x[k+1]=x[k]+alpha*d
    D<-calculate_direction(Y,X,lambda1,lambda2)#计算下一次迭代的方向
   }
  
  return (X)
}