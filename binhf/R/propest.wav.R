"propest.wav" <-
function (proportion=P2,binsize=1,length=256,times=100,meth="u",van=6,fam="DaubLeAsymm",min.level=3)
{

n<-length
binary<-NULL

mise<-matrix(0,times,n)
misea<-mise

x<-(1:n)/n
x<-x-1/(2*n)

y<-proportion(x)

if(length(y)==1){
cat("Proportion produces only one value...converting to length",n,"\n")  
y<-rep(y,length=n)
}


binmat<-matrix(0,times,n)
estmat<-matrix(0,times,n)
estmata<-estmat

for (k in 1:times){

for (j in 1:n){
binary[j]<-rbinom(1,binsize,y[j]/binsize)
}

binmat[k,]<-binary

fhat<-hfdenoise.wav(binary, binsize,transform="binhf",meth, van = van, fam = fam, min.level = min.level,coarse=FALSE)
fhata<-hfdenoise.wav(binary, binsize, transform="ansc",meth, van = van, fam = fam, min.level = min.level,coarse=FALSE)


estmat[k,]<-fhat
estmata[k,]<-fhata

mise[k,]<-(fhat[1:n]-y[1:n])^2
misea[k,]<-(fhata[1:n]-y[1:n])^2

}

amse<-sum(mise)/(times*n)
amsea<-sum(misea)/(times*n)

meanfhat<-apply(estmat,2,mean)
meanfhata<-apply(estmata,2,mean)

return(list(x=x,y=y,b=binmat,e=estmat,ea=estmata,meanfhat=meanfhat,meanfhata=meanfhata,amse=amse,amsea=amsea))

}

