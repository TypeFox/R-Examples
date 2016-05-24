postNullFun <-
function(x1,x2,n1,n2,m=5,w1,w2,gfun){

wp1<-w1
wp2<-w2

wstar1<-function(q){w1[q+1]*(m+1)*choose(n1,x1)*choose(m,q)/((m+n1+1)*choose(m+n1,x1+q))}
wstar2<-function(q){w2[q+1]*(m+1)*choose(n2,x2)*choose(m,q)/((m+n2+1)*choose(m+n2,x2+q))}

q.all<-0:m
wp1<-sapply(q.all,wstar1)
wp2<-sapply(q.all,wstar2)

wp1<-wp1/sum(wp1)
wp2<-wp2/sum(wp2)

b<-function(i,j){
f1 <- function(y){return(pbeta(gfun(y),x2+j+1,m+n2-x2-j+1)*dbeta(y,x1+i+1,m+n1-x1-i+1))}
b<-integrate(f1,1e-4,1-1e-4)$value
}

b.vec<-Vectorize(b)

i_test<-0:m
j_test<-0:m

B<-outer(i_test,j_test,b.vec)


numBFtemp<-t(wp1)%*%B%*%wp2


if(numBFtemp>=1-1e-5){numBFtemp=1}
if(numBFtemp<=1e-5){numBFtemp=0}
postNull=numBFtemp

return(postNull)
}

