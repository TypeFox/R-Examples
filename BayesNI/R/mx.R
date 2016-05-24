mx <-
function(x1,x2,n1,n2,w1,w2,m){
  
c<-function(i,j){
choose(n1,x1)*choose(n2,x2)*beta(x1+i+1,n1-x1+m-i+1)*beta(x2+j+1,n2-x2+m-j+1)/
(beta(i+1,m-i+1)*beta(j+1,m-j+1))}

c.vec<-Vectorize(c)

i_test<-0:m
j_test<-0:m

C<-outer(i_test,j_test,c.vec)

mx<-t(w1)%*%C%*%w2

return(mx)
}

