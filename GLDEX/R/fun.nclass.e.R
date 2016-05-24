"fun.nclass.e" <-
function(x){

len<-length(x)

# matrix of estimated mean and variance:
m.e<-t(sapply(1:length(x), function(i,x) fun.disc.estimation(x,i),x))
# matrix of actual mean and variance of data x.
m.a<-matrix(rep(c(mean(x),var(x)),len),nrow=len,byrow=TRUE)

r<-rowSums((m.e-m.a)^2)

    r<-which(r==min(r))


return(r)
}

