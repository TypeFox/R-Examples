calculate_direction <-
function(Y,X,lambda1,lambda2){
  n<-dim(Y)[1];p<-dim(Y)[2];set.seed(123);D<-matrix(rnorm(n*p),n,p)#³õÊ¼»¯D
  tem<-.C("abcd",as.double(X),as.double(Y),as.double(D),as.double(D),
        as.double(lambda1),as.double(lambda2),
        as.integer(n),as.integer(p))
  temp<-D;
  D<-t(matrix(tem[[4]],p,n));

  while(max(abs(temp-D))>5*10^-3 &sum(is.nan(D)*1)==0 & sum((D==Inf)*1)==0 & sum((D==-Inf)*1)==0){
       temp<-D
       tem<-.C("abcd",as.double(X),as.double(Y),as.double(D),as.double(D),
             as.double(lambda1),as.double(lambda2),
            as.integer(n),as.integer(p))
       D<-t(matrix(tem[[4]],p,n));
    }
   return(temp)
}