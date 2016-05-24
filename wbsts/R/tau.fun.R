tau.fun <-
function(y){
  n=ifelse(length(y)<6000,length(y),6000);J=floor(log(n,2))
  scales<-J-1:floor(3*0.75*log(log(n)))
  sc=min(length(scales),6)
  coef=c(9.795793e-01, -1.084024e-04,  3.112115e+01,  1.131661e-08 )
  j1=sum(coef*c(1,n,1/n,n^2))
  coef=c(1.132108e+00, -1.110789e-04,  3.286444e+01,  1.093054e-08 )
  j2=sum(coef*c(1,n,1/n,n^2))
  coef=c(1.518703e+00, -1.585613e-04,  2.118412e+01,  1.773870e-08 )
  j3=sum(coef*c(1,n,1/n,n^2))
  coef=c(2.059255e+00, -1.798321e-04,  7.827924e+00,  1.632834e-08 )
  j4=sum(coef*c(1,n,1/n,n^2))
  coef=c(1.867909e+00 , 1.633426e-04,  4.867038e+02, -2.425376e-08 )
  j5=sum(coef*c(1,n,1/n,n^2))
  coef=c(5.782118e+00, -3.728447e-04, -5.025129e+03,  4.673449e-10 )
  j6=2.5  
  return(c(j1,j2,j3,j4,j5,j6))  
}
