
########################################################################
eucl.dist<-function(data,allloc=FALSE,distance="eucl"){
        #Calculates Euclidian distances from genetic data
nbloc<-dim(data)[2]-2
nbpop<-max(data[,1])
euclmat<-matrix(nrow=(nbpop*(nbpop-1)/2),ncol=nbloc)
for (i in 1:nbloc){
     eff<-apply(table(data[,(i+2)],data[,1]),2,sum,na.rm=TRUE)
     freq<-sweep(table(data[,(i+2)],data[,1]),2,eff,FUN="/")
     eucl<-t(freq)%*%freq
     dum<-0
     for (j in 2:nbpop){
     for (k in 1:(j-1)){
         dum<-dum+1
         euclmat[dum,i]<-eucl[j,j]+eucl[k,k]-2*eucl[j,k]
     }
     }
}
x<-apply(euclmat,1,sum,na.rm=TRUE)
euclmat<-cbind(euclmat,x)
if (allloc) {return(euclmat^0.5)} else {return(x^0.5)}
}
########################################################################
eucl.dist.trait<-function(data){
#calculate euclidian distance between populations for one trait
nbpop<-max(data[,1])
euclmat<-vector(length=(nbpop*(nbpop-1)/2))
eff<-tapply(data[,2],data[,1],mean)
     dum<-0
     for (j in 2:nbpop){
     for (k in 1:(j-1)){
         dum<-dum+1
         euclmat[dum]<-sqrt((eff[j]-eff[k])^2)
     }
     }
return(euclmat)
}
########################################################################
#eucl.dist.mtrait<-function(data,corr=FALSE){
#calculate euclidian distance between populations for multiple traits
#nbpop<-max(data[,1])
#nbtraits<-dim(data)[2]-1
#if (corr)  {data<-cbind(data[,1],scale(data[,-1]))}
#euclmat<-matrix(nrow=(nbpop*(nbpop-1)/2),ncol=nbtraits)
#for (i in 1:nbtraits){
#eff<-tapply(data[,i+1],data[,1],mean)
#     dum<-0
#     for (j in 2:nbpop){
#     for (k in 1:(j-1)){
#         dum<-dum+1
#         euclmat[dum,i]<-(eff[j]-eff[k])^2
#     }
#     }
#}
#x<-(apply(euclmat,1,sum))^0.5
#return(x)
#}
########################################################################
cfe.dist<-function(data,allloc=FALSE,distance="cfe"){
        #Calculates Cavalli-Sforza & Edwards Chord distance, according to Nei, 1987 (p 216, eq 9.15)
        #Note  I used d=(1-cos(theta))^0.5, for a distance between 0 & 1, see Nei
nbloc<-dim(data)[2]-2
nbpop<-max(data[,1])
cfemat<-matrix(nrow=(nbpop*(nbpop-1)/2),ncol=nbloc)
for (i in 1:nbloc){
     eff<-apply(table(data[,(i+2)],data[,1]),2,sum,na.rm=TRUE)
     freq<-sweep(table(data[,(i+2)],data[,1]),2,eff,FUN="/")
     cfe<-(1-t(freq^0.5) %*% freq^0.5)^0.5
     dum<-0
     for (j in 2:nbpop){
     for (k in 1:(j-1)){
         dum<-dum+1
         cfemat[dum,i]<-cfe[j,k]
     }
     }
}
x<-apply(cfemat,1,mean,na.rm=TRUE)
cfemat<-cbind(cfemat,x)
if (allloc) {return(cfemat)} else {return(x)}
}
########################################################################
da.dist<-function(data,allloc=FALSE,distance="da"){
        #Calculates Nei DA distance (Nei, 1986, p 216 eq 9.16)
nbloc<-dim(data)[2]-2
nbpop<-max(data[,1])
damat<-matrix(nrow=(nbpop*(nbpop-1)/2),ncol=nbloc)
for (i in 1:nbloc){
     eff<-apply(table(data[,(i+2)],data[,1]),2,sum,na.rm=TRUE)
     freq<-sweep(table(data[,(i+2)],data[,1]),2,eff,FUN="/")
     da<-1-t(freq^0.5) %*% freq^0.5
     dum<-0
     for (j in 2:nbpop){
     for (k in 1:(j-1)){
         dum<-dum+1
         damat[dum,i]<-da[j,k]
     }
     }
}
x<-apply(damat,1,mean,na.rm=TRUE)
damat<-cbind(damat,x)

if (allloc) {return(damat)} else {return(x)}
}

########################################################################
nei.dist<-function(data){
        #estimates Nei unbiased genetic distance (Nei, 1978)
nbloc<-dim(data)[2]-2
nbpop<-max(data[,1])
neivecn<-vector(length=(nbpop*(nbpop-1)/2))
neivecd1<-vector(length=(nbpop*(nbpop-1)/2))
neivecd2<-vector(length=(nbpop*(nbpop-1)/2))
for (i in 1:nbloc){
     ta<-table(data[,(i+2)],data[,1])
     eff<-apply(ta,2,sum,na.rm=TRUE)
     freq<-sweep(ta,2,eff,FUN="/")
     nei<-t(freq) %*% freq
     dum<-0
     for (j in 2:nbpop){
     for (k in 1:(j-1)){
         dum<-dum+1
         neivecn[dum]<-neivecn[dum]+nei[j,k]
         neivecd1[dum]<-neivecd1[dum]+(eff[j]*nei[j,j]-1)/(eff[j]-1)
         neivecd2[dum]<-neivecd2[dum]+(eff[k]*nei[k,k]-1)/(eff[k]-1)
     }
     }
}
return(-log(neivecn/(neivecd1*neivecd2)^0.5))
}
#################################################################
#sad.dist<-function(data){
#shared allele distance. Need to be completed
#x<-dim(data)
#if (max(data[,2],na.rm=TRUE)<1000000) modulo=1000
#if (max(data[,2],na.rm=TRUE)<10000) modulo<-100
#if (max(data[,2],na.rm=TRUE)<100) modulo<-10
#firstal<-pmin(data[,-1] %/% modulo,data[,-1] %% modulo)
#secal<-pmax(data[,-1]%/% modulo,data[,-1] %% modulo)
#nind<-x[1]
#das<-vector(length=nind*(nind-1)/2)
#cum<-0
#for (i in 2:nind){
#for (j in 1:(i-1)){
#cum<-cum+1
#}
#}
#}

#################################################################
vec2mat<-function(x){
	#fill a lower triangular matrix from a vector and copy it to upper triangle
nn<-length(x)
n<-(1+(1+8*nn)^0.5)/2 #dim of the matrix
mat<-matrix(rep(0,n*n),ncol=n,nrow=n)
cum<-0
for (i in 2:n) {
	for (j in 1:(i-1)){
		cum<-cum+1
		mat[i,j]<-x[cum]
		mat[j,i]<-mat[i,j]
	}
}
return(mat)
}
#################################################################
mat2vec<-function(mat){
#transform lower triangular matrix in vector 1.2,1.3,2.3,1.4,2.4,3.4 etc...
n<-dim(mat)[2]
nn<-n*(n-1)/2
x<-vector(length=nn)
cum<-0
for (i in 2:n){
for (j in 1:(i-1)){
cum<-cum+1
x[cum]<-mat[i,j]
}}
return(x)
}

#################################################################
pcoa<-function(mat,plotit=TRUE,...){
#principal coordinates analysis
#as described in Legendre & lengendre Numerical Ecology, p 426
n<-dim(mat)[1]
a<--0.5*mat^2
mr<-apply(a,1,mean)
mc<-apply(a,2,mean)
mm<-mean(a)
delta<-sweep(a,1,mr)
delta<-sweep(delta,2,mc)
delta<-delta+mm
eig.delta<-eigen(delta)
vec.delta<-eig.delta$vectors*matrix(rep(eig.delta$values^0.5,n),nrow=n,byrow=TRUE)
inertia<-eig.delta$values/sum(eig.delta$values)
if (plotit) {
par(mfrow=c(2,2))
plot(1:n,inertia);abline(h=0)
plot(vec.delta[,1],vec.delta[,2],type="n",xlab=paste("First axis (inertia=",round(inertia[1],2),")",sep=""),ylab=paste("Second axis (inertia=",round(inertia[2],2),")",sep=""))
text(vec.delta[,1],vec.delta[,2],...)
plot(vec.delta[,1],vec.delta[,3],type="n",xlab=paste("First axis (inertia=",round(inertia[1],2),")",sep=""),ylab=paste("Third axis (inertia=",round(inertia[3],2),")",sep=""))
text(vec.delta[,1],vec.delta[,3],...)

par(mfrow=c(1,1))
}

naxes<-min(100,sum(eig.delta$values>1.0e-7)) #to limit to max 100 axes
dum<-vec.delta[,1:(naxes-1)]
nn<-n*(n-1)/2
eucl<-matrix(rep(0,nn*(naxes-1)),ncol=naxes-1,byrow=TRUE)
cum<-0
for (i in 2:n){
for (j in 1:(i-1)){
cum<-cum+1
eucl[cum,]<-abs(dum[i,]-dum[j,])
}
}

list(valp=eig.delta$values[1:(naxes-1)],vecp=dum,eucl=eucl)
}

#################################################################
