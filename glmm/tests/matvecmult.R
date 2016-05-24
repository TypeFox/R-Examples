library(glmm)
set.seed(1234)

 nrow <- 5
 ncol <- 3

 a <- matrix(rnorm(nrow * ncol), nrow = nrow)
 b <- rnorm(ncol)

 foo1 <- as.numeric(a %*% b)

 mout3 <- .C("matvecmult", a = as.double(a), b = as.double(b),
     nrow = as.integer(nrow), ncol = as.integer(ncol), result = double(nrow))
 identical(foo1, mout3$result)
 all.equal(foo1, mout3$result)

a<-matrix(1:8,nrow=2)
b<-matrix(1:4,nrow=2)
right<-t(a)%*%b

stuff<-.C("matTmatmult",as.double(a),as.double(b),as.integer(2),as.integer(4),as.integer(2),double(8))[[6]]
stuff2<-matrix(stuff,byrow=F,nrow=4) #at least
all.equal(right,stuff2)

#make sure the function to add a vector up works
a<-1:5
right<-sum(a)
stuff<-.C("sumup",as.double(a),as.integer(length(a)),double(1))[[3]]
all.equal(stuff,right)

#make sure the function to do a-b works
a<-6:10
b<-1:5
right<-a-b
stuff<-.C("subvec",as.double(a),as.double(b),length(a),result=double(length(a)))
all.equal(right,stuff$result)

# make sure matTvecmult works
a<-matrix(1:8,nrow=2)
b<-1:2
right<-t(a)%*%b
alsoright<-.C("matvecmult",as.double(t(a)),as.double(b),as.integer(nrow(t(a))),as.integer(ncol(t(a))),double(4))[[5]]
alsoright<-matrix(alsoright,ncol=1)
all.equal(right,alsoright)
stuff<-matrix(.C("matTvecmult",as.double(a),as.double(b),as.integer(nrow(a)),as.integer(ncol(a)),double(4))[[5]],ncol=1)
all.equal(stuff,right)


