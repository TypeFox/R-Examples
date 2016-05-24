setClass("fseries",representation(p="matrix",ap="matrix",m="matrix"),
validity=function(object){
if(any(ncol(object@ap)!=1,nrow(object@ap)!=nrow(object@p)))
stop("Wrong dimensions of matrix")
if(nrow(object@p)!=nrow(unique(object@p)))
stop("Rows of the coefficient matrix are not unique")
if(any(sum(matrix(as.integer(object@p),nrow(object@p),ncol(object@p))==object@p)!=nrow(object@p)*ncol(object@p),
sum(matrix(as.integer(object@p),nrow(object@p),ncol(object@p))>=0)!=nrow(object@p)*ncol(object@p)))
stop("Exponents are not natural")
return(TRUE)})

setMethod("initialize",c("fseries"),function(.Object,p,ap)
{
.Object@ap<-ap
.Object@p<-p
validObject(.Object)
w=which(.Object@ap[,1]!=0)
m<-as.matrix(rowSums(.Object@p))
ord=order(m[w])
.Object@p<-matrix(p[w[ord],],length(w),ncol(p))
.Object@ap<-matrix(ap[w[ord],],length(w),1)
.Object@m<-matrix(m[w[ord],],length(w),1)
colnames(.Object@p)<-paste("p",1:ncol(.Object@p),sep="")
colnames(.Object@ap)<-"ap"
colnames(.Object@m)<-"|p|"
rownames(.Object@p)<-NULL
rownames(.Object@ap)<-NULL
rownames(.Object@m)<-NULL
return(.Object)
})

setGeneric("rfseries",function(var,cf,k,m){standardGeneric("rfseries")})

setMethod("rfseries",c("numeric","numeric","numeric"),function(var,cf,k,m)
{
if((k+1)^var<cf)
stop("The parameter k is too low")
p=matrix(as.integer(runif(var*cf,0,k+1)),cf,var)
ap=matrix(2*rbinom(cf,1,1/2)-1,cf,1)*matrix(as.integer(runif(var*cf,1,m)),cf,1)
n=nrow(unique(p))
while(n!=nrow(p))
{
	p=unique(p)
	p=rbind(p,matrix(as.integer(runif(var*(cf-n),0,k+1)),cf-n,var))
	n=nrow(unique(p))	
}
c=new("fseries",p,ap)
return(c)
})

setMethod("print","fseries",function(x){
if(nrow(x@p)==0)return(print(0))
max=max(nchar(as.character(x@ap)))+1
for(i in 1:nrow(x@p)){
wartosc=paste(x@ap[i,]," ",sep="")
ma=nchar(wartosc)
if(max-ma>0)
for(k in 1:eval(max-ma))wartosc=paste(" ",wartosc,sep="")
for(j in 1:ncol(x@p)){
wartosc=paste(wartosc,"X",j,"^",x@p[i,j]," ",sep="")}
cat(wartosc,"\n")}})

setMethod("show","fseries",function(object){
if(nrow(object@p)==0)return(print(0))
max=max(nchar(as.character(object@ap)))+1
for(i in 1:nrow(object@p)){
wartosc=paste(object@ap[i,]," ",sep="")
ma=nchar(wartosc)
if(max-ma>0)
for(k in 1:eval(max-ma))wartosc=paste(" ",wartosc,sep="")
for(j in 1:ncol(object@p)){
wartosc=paste(wartosc,"X",j,"^",object@p[i,j]," ",sep="")}
cat(wartosc,"\n")}})

setMethod("+",c("fseries","fseries"),function(e1,e2)
{
if(ncol(e1@p)<ncol(e2@p)){
e1@p=cbind(e1@p,matrix(0,nrow(e1@p),ncol(e2@p)-ncol(e1@p)))
c=new("fseries",e1@p,e1@ap)}
if(ncol(e1@p)>ncol(e2@p)){
e2@p=cbind(e2@p,matrix(0,nrow(e2@p),ncol(e1@p)-ncol(e2@p)))
c=new("fseries",e2@p,e2@ap)}
n=nrow(e1@p)
colp=ncol(e1@p)
dlp=length(e1@p)
pcx=0
pcy=0
dodane=0
k=0
for(i in 1:n){
r=which(e1@m[i]==e2@m)
rowny=integer(0)
dlr=length(r)
if(dlr==1)
rowny=r[colSums(as.matrix(e1@p[i,]==e2@p[r,]))==colp]
if(dlr>1)
rowny=r[rowSums(matrix(rep(e1@p[i,]),dlr,colp,T)==e2@p[r,])==rep(colp,dlr)]
if(dlr!=0&length(rowny)!=0){
k=1
pcx=c(pcx,i)
pcy=c(pcy,rowny)
dodane=c(dodane,e1@ap[i]+e2@ap[rowny])
dodane=as.matrix(dodane)
}}
if(k==1)
c=new("fseries",p=rbind(matrix(e1@p[-pcx,],n-length(pcx)+1,colp),matrix(e2@p[-pcy,],nrow(e2@p)-length(pcy)+1,colp),matrix(e1@p[pcx,],length(pcx)-1,colp)),
ap=as.matrix(c(e1@ap[-pcx,],e2@ap[-pcy,],dodane[2:length(dodane),1])))
else
c=new("fseries",p=rbind(e1@p,e2@p),ap=rbind(e1@ap,e2@ap))
return(c)
})

setMethod("-",c("fseries"),function(e1)
{
ap=(-1)*e1@ap
c=new("fseries",e1@p,ap)
return(c)
})

setMethod("-",c("fseries","fseries"),function(e1,e2)
{
c=e1+(-e2)
return(c)
})

setMethod("+",c("fseries","numeric"),function(e1,e2)
{
e1[rep(0,ncol(e1@p))]=e1[rep(0,ncol(e1@p))]+e2
c=new("fseries",e1@p,e1@ap)
return(c)
})

setMethod("+",c("numeric","fseries"),function(e1,e2)
{
e2[rep(0,ncol(e2@p))]=e2[rep(0,ncol(e2@p))]+e1
c=new("fseries",e2@p,e2@ap)
return(c)
})

setMethod("-",c("fseries","numeric"),function(e1,e2)
{
e1[rep(0,ncol(e1@p))]=e1[rep(0,ncol(e1@p))]-e2
c=new("fseries",e1@p,e1@ap)
return(c)
})

setMethod("-",c("numeric","fseries"),function(e1,e2)
{
e2[rep(0,ncol(e2@p))]=e1-e2[rep(0,ncol(e2@p))]
c=new("fseries",e2@p,e2@ap)
return(c)
}) 

setMethod("*",c("fseries","fseries"),function(e1,e2)
{
if(ncol(e1@p)<ncol(e2@p)){
e1@p=cbind(e1@p,matrix(0,nrow(e1@p),ncol(e2@p)-ncol(e1@p)))
c=new("fseries",e1@p,e1@ap)}
if(ncol(e1@p)>ncol(e2@p)){
e2@p=cbind(e2@p,matrix(0,nrow(e2@p),ncol(e1@p)-ncol(e2@p)))
c=new("fseries",e2@p,e2@ap)}
p=matrix(e1@p[1,],nrow(e2@p),ncol(e1@p),T)+e2@p
ap=e1@ap[1,]*e2@ap
d=new("fseries",p,ap)
for(i in 2:nrow(e1@p))
{
p=matrix(e1@p[i,],nrow(e2@p),ncol(e1@p),T)+e2@p
ap=e1@ap[i,]*e2@ap
c=new("fseries",p,ap)
d=d+c
}
return(d)})

setMethod("*",c("numeric","fseries"),function(e1,e2)
{if(e1==0)return(0)
ap=e1*e2@ap
c=new("fseries",e2@p,ap)
return(c)})

setMethod("*",c("fseries","numeric"),function(e1,e2)
{if(e2==0)return(0)
ap=e2*e1@ap
c=new("fseries",e1@p,ap)
return(c)})

setMethod(
"[",
signature=c("fseries","numeric"),
function(x,i){
if(length(i)!=ncol(x@p))
stop("The vector length is incorrect")
k=which(rowSums(matrix(i,nrow(x@p),ncol(x@p),T)==x@p)==ncol(x@p))
if(length(k)!=0)return(as.numeric(x@ap[k,]))
else return(0)})

setReplaceMethod("[",signature=c("fseries","numeric","missing","numeric"),function(x,i,j,value){
if(length(i)!=ncol(x@p))
stop("The vector length is incorrect")
k=which(rowSums(matrix(i,nrow(x@p),ncol(x@p),T)==x@p)==ncol(x@p))
if(length(k)!=0)x@ap[k,]=value
else{
x@p=rbind(x@p,i)
x@ap=rbind(x@ap,value)}
x=new("fseries",x@p,x@ap)
return(x) 
})

setMethod(
"[",
signature=c("fseries","matrix"),
function(x,i){
if(ncol(i)!=ncol(x@p))
stop("The matrix dimensions is incorrect")
wynik=matrix(0,nrow(i),1)
for(j in 1:nrow(i)){
k=which(rowSums(matrix(i[j,],nrow(x@p),ncol(x@p),T)==x@p)==ncol(x@p))
if(length(k)!=0) wynik[j,]=as.numeric(x@ap[k,])
else wynik[j,]=0}
return(wynik)})

setReplaceMethod("[",signature=c("fseries","matrix","missing","matrix"),function(x,i,j,value){
if(nrow(i)!=nrow(value)||ncol(i)!=ncol(x@p))
stop("The matrix dimensions is incorrect")
for(j in 1:nrow(i)){
k=which(rowSums(matrix(i[j,],nrow(x@p),ncol(x@p),T)==x@p)==ncol(x@p))
if(length(k)!=0) x@ap[k,]=value[j,]
else {
x@p=rbind(x@p,i[j,])
x@ap=rbind(x@ap,value[j,])}}
x=new("fseries",x@p,x@ap)
return(x) 
})

aux11<-function(x){
i1=x[1]
w1=matrix(0,i1+1,1)
w2=matrix(0,i1+1,1)
k=0
for(j in i1:0){
k=k+1
w1[k,1]=j
w2[k,1]=i1-j
}
wynik=list(x=w1,y=w2)
return(wynik)
}

aux12<-function(x){
i1=x[1]
i2=x[2]
w1=matrix(0,(i1+1)*(i2+1),2)
w2=matrix(0,(i1+1)*(i2+1),2)
k=0
for(i in i2:0){
for(j in i1:0){
k=k+1
w1[k,1:2]=c(j,i)
w2[k,1:2]=c(i1-j,i2-i)
}}
wynik=list(x=w1,y=w2)
return(wynik)
}

aux13<-function(x){
i1=x[1]
i2=x[2]
i3=x[3]
w1=matrix(0,(i1+1)*(i2+1)*(i3+1),3)
w2=matrix(0,(i1+1)*(i2+1)*(i3+1),3)
w3=matrix(0,(i1+1)*(i2+1)*(i3+1),3)
k=0
for(z in i3:0){
for(i in i2:0){
for(j in i1:0){
k=k+1
w1[k,1:3]=c(j,i,z)
w2[k,1:3]=c(i1-j,i2-i,i3-z)
}}}
wynik=list(x=w1,y=w2)
return(wynik)
}

aux21<-function(i,p,ap){
wynik=matrix(0,length(i),1)
for(j in 1:length(i)){
k=which(rowSums(matrix(i[j],nrow(p),1,T)==p)==1)
if(length(k)!=0) wynik[j,]=as.numeric(ap[k,])
else wynik[j,]=0}
return(wynik)
}

aux2<-function(i,p,ap){
if(is.vector(i)){
k=which(rowSums(matrix(i,nrow(p),ncol(p),T)==p)==ncol(p))
if(length(k)!=0)return(as.numeric(ap[k,]))
else return(0)}
else{
wynik=matrix(0,nrow(i),1)
for(j in 1:nrow(i)){
k=which(rowSums(matrix(i[j,],nrow(p),ncol(p),T)==p)==ncol(p))
if(length(k)!=0) wynik[j,]=as.numeric(ap[k,])
else wynik[j,]=0}
return(wynik)
}}

setMethod("&",signature=c("fseries","numeric"),function(e1,e2)
{
if(ncol(e1@p)==1){
if(e1[0]==0)
stop("The formal series is not invertible")
for(i in 0:e2){
if(i==0){
ap=matrix(1/e1[0],1,1)
p=matrix(0,1,1)}
else{
l=aux11(i)
ap=rbind(ap,-(t(aux21(l$x[2:nrow(l$x),],p,ap))%*%e1[matrix(l$y[2:nrow(l$y),],length(l$y)-1,1)])/e1[0])
p=rbind(p,i)
}}
c=new("fseries",p,ap)
return(c)}

if(ncol(e1@p)==2){
if(e1[c(0,0)]==0)
stop("The formal series in invertible")
for(i in 0:e2){
for(j in 0:e2){
if(i==0&&j==0){
ap=matrix(1/e1[c(0,0)],1,1)
p=matrix(0,1,2)}
else{
l=aux12(c(j,i))
ap=rbind(ap,-(t(aux2(l$x[2:nrow(l$x),],p,ap))%*%e1[l$y[2:nrow(l$y),]])/e1[c(0,0)])
p=rbind(p,c(j,i))
}}}
c=new("fseries",p,ap)
return(c)}

if(ncol(e1@p)==3){
if(e1[c(0,0,0)]==0)
stop("The formal series in invertible")
for(z in 0:e2){
for(i in 0:e2){
for(j in 0:e2){
if(all(z==0,i==0,j==0)){
ap=matrix(1/e1[c(0,0,0)],1,1)
p=matrix(0,1,3)}
else{
l=aux13(c(j,i,z))
ap=rbind(ap,-(t(aux2(l$x[2:nrow(l$x),],p,ap))%*%e1[l$y[2:nrow(l$y),]])/e1[c(0,0,0)])
p=rbind(p,c(j,i,z))
}}}}
c=new("fseries",p,ap)
return(c)}
else stop("For more than 3 variables the method is not implemented")
})

setMethod("/",signature=c("fseries","fseries"),function(e1,e2)
{
c=e1*(e2&5)
return(c)
})

setMethod("/",signature=c("numeric","fseries"),function(e1,e2)
{
c=e1*(e2&5)
return(c)
})

setMethod("/",signature=c("fseries","numeric"),function(e1,e2)
{
ap=e1@ap/e2
c=new("fseries",e1@p,ap)
return(c)
})