`weightedGramSchmidt` <-
function(x,w) {
  ss<-NULL;
  for (j in 1:dim(x)[2]) {
  if (j > 1) {xx<-x[,1:(j-1)]; x[,j]<-x[,j]-xx%*%(crossprod(xx,(w*x[,j])))}
  s<-sqrt(sum(w*x[,j]^2)); ss<-c(ss,s); x[,j]<-x[,j]/s;}
  list(pol=x,fac=ss)
}

