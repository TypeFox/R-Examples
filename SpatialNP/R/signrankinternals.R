norm<-function(X) 
{
X<-as.matrix(X)
.C("norming", as.double(X), as.integer(dim(X)), res=double(dim(X)[1]),PACKAGE="SpatialNP")$res
}

pairdiff<-function(X)
{
X<-as.matrix(X)
d<-dim(X)
matrix(.C("pairdiff", as.double(X),as.integer(d), res=double(choose(d[1],2)*d[2]),PACKAGE="SpatialNP")$res,ncol=d[2],byrow=T)
}

pairsum<-function(X)
{
X<-as.matrix(X)
d<-dim(X)
matrix(.C("pairsum", as.double(X),as.integer(d), res=double(choose(d[1],2)*d[2]),PACKAGE="SpatialNP")$res,ncol=d[2],byrow=T)
}

sumsignout<-function(X)
{
#have to be careful with zero divisions:
ind<-numeric(0)
ind<-apply(X,1,setequal,y=0)
if(length(ind)>0) 
 X<-X[!ind,]

d<-dim(X)
tmp<-matrix(.C("sum_of_sign_outers", as.double(X),as.integer(d), res=double(d[2]^2),PACKAGE="SpatialNP")$res,ncol=d[2],byrow=T)
}

ranks<-function(X)
{
d<-dim(X)
if(any(d[2]==1,is.vector(X))) return(rank(X))
tmp<-matrix(.C("spatial_ranks", as.double(X),as.integer(d), res=double(d[1]*d[2]),PACKAGE="SpatialNP")$res,ncol=d[2],byrow=T)

#need to check whether multiple values caused a problem:
ind<-which(!is.finite(tmp[,1]))
if(length(ind)==0) return(tmp)
#else the problem ranks need to be recalculated:
for(i in ind) {
 s<-sweep(X,2,X[i,]) 
 r<-norm(s)
 r[r==0]<-1
 tmp[i,]<-apply(sweep(s,1,r,"/"),2,mean)
}
tmp
}

signranks<-function(X)
{
X<-as.matrix(X)
d<-dim(X)
tmp<-matrix(.C("signed_ranks", as.double(X),as.integer(d), res=double(d[1]*d[2]),PACKAGE="SpatialNP")$res,ncol=d[2],byrow=T)
#as in ranks, possibly need to recalculate stuff
ind<-which(!is.finite(tmp[,1]))
if(length(ind)==0) return(tmp)
#else 
for(i in ind) {
 sm<-sweep(X,2,X[i,]) 
 sp<-sweep(X,2,X[i,],"+")
 rm<-norm(sm)
 rp<-norm(sp)
 rm[rm==0]<-1
 rp[rp==0]<-1
 tmp[i,]<-(apply(sweep(sm,1,rm,"/"),2,mean)+apply(sweep(sp,1,rp,"/"),2,mean))/2
}
tmp
}

Q2internal<-function(X)
{
X<-as.matrix(X)
d<-dim(X)
.C("Q2internals",as.double(X),as.integer(d), res=double(d[2]^2+d[2]^4),PACKAGE="SpatialNP")$res
}

gen.inv<-function(M)
{
p<-sqrt(dim(M)[1])
eig<-eigen(M)
eig$values[1:((p+2)*(p-1)/2)]<-eig$values[1:((p+2)*(p-1)/2)]^-1
res<-eig$vectors%*%diag(eig$values)%*%t(eig$vectors)
Re(res) #because it really is real!
}

Cpp<-function(p)
{
I<-diag(p)
vecI<-I
dim(vecI)<-NULL
J<-vecI%*%t(vecI)
K<-matrix(0,ncol=p^2,nrow=p^2)
for (i in 1:p) for (j in 1:p) K<-K+kronecker(outer(I[,i],I[,j]),outer(I[,j],I[,i]))
(diag(p^2)+K)/2-J/p
}

covshape<-function(x,determ=NULL,trace=NULL,first=NULL) 
{
 to.shape(cov(x),determ=determ,trace=trace,first=first)
}

mat.sqrt<-function(A)
# Returns the square root matrix of the given matrix.
{
 eig<-eigen(A)
 eig$vectors%*%(diag(eig$values))^(1/2)%*%t(eig$vectors)
}

mat.norm<-function(A)
# Returns the matrix norm of the given matrix.
{
 sqrt(sum(diag(t(A)%*%A)))
}

center.step<-function(X,y)
# Returns AVE[U(xi-y)]/AVE[|xi-y|⁻¹]
{
X<-as.matrix(X)
d<-dim(X)
c(.C("center_step", as.double(X),as.integer(d),as.double(y), res=double(d[2]),PACKAGE="SpatialNP")$res)
}

spat.med.step<-function(X,y)
# Returns y+AVE[U(xi-y)]/AVE[|xi-y|⁻¹]
{
X<-as.matrix(X)
d<-dim(X)
c(.C("spat_med_step", as.double(X),as.integer(d),as.double(y), res=double(d[2]),PACKAGE="SpatialNP")$res)
}

SSCov.hub<-function(X,V,c.s,sig.s)
{
X<-as.matrix(X)
d<-dim(X)
matrix(.C("symm_huber", as.double(X), as.double(V), as.integer(d), as.double(c.s), as.double(sig.s),res=double(d[2]^2),PACKAGE="SpatialNP")$res, ncol=d[2],byrow=T)/(d[1]*(d[1]-1)/2)
}

hl.center.step <- function(X,y)
{
X<-as.matrix(X)
d<-dim(X)
c(.C("hl_center_step", as.double(X),as.integer(d),as.double(y),res=double(d[2]),PACKAGE="SpatialNP")$res)
}

hl.loc.step<-function(X,y)
{
X<-as.matrix(X)
d<-dim(X)
c(.C("hl_loc_step", as.double(X),as.integer(d),as.double(y), res=double(d[2]),PACKAGE="SpatialNP")$res)
}

spat.median <- function (X, init=NULL, steps=Inf, maxiter=500, eps=1e-06, na.action=na.fail)
{ 
    X <- na.action(X) 

    if (dim(X)[2] == 1) 
        return(median(X))
    if (is.null(init)) 
        y <- apply(X, 2, median)
    else y <- init
    if(is.finite(steps)) maxiter<-Inf
    iter <- 0
    differ <- Inf
    
   while (TRUE) {
        if(iter>=steps) return(y)
        if (iter == maxiter) 
            stop("maxiter reached without convergence")
        y.new <- spat.med.step(X,y)
        differ <- sqrt(sum((y -y.new)^2))
        iter <- iter + 1
        if(differ<eps) return(y.new)
        y <- y.new
    }
}

hl.location<-function(X, init=NULL, steps=Inf, maxiter=500, eps=1e-06, na.action=na.fail)
{ 
    X <- na.action(X)
    if (dim(X)[2] == 1) 
        return(median(X))
    if (is.null(init)) 
        y <- apply(X, 2, median)
    else y <- init
    if(is.finite(steps)) maxiter<-Inf
    iter <- 0
    differ <- Inf
   
   while (TRUE) {
        if(iter>=steps) return(y)
        if (iter == maxiter) 
            stop("maxiter reached without convergence")
        y.new <- y+hl.center.step(X,y)
        differ <- sqrt(sum((y -y.new)^2))
        iter <- iter + 1
        if(differ<eps) return(y.new)
        y <- y.new
    }
}


