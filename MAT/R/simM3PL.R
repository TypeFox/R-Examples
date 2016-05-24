simM3PL <-
function(ipar,cors,p,n.simulee=100,D=1.0,easiness=T,seed=NULL) {
call<-match.call()

ni<-nrow(ipar)

if (p==1) sigma<-1 
      else if (p>1) {
sigma<-as.matrix(cors)
sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(sigma)]
if (dim(sigma)[1]!=dim(sigma)[2] || p!= dim(sigma)[1]) stop("ERROR: p and cors non-conforming")
} 
      else stop("ERROR: p and cors non-conforming")

if (p==1) sigma.inv=1
else sigma.inv<-solve(sigma)

if (is.data.frame(ipar)) ipar<-cbind(as.matrix(ipar[paste("a",1:p,sep="")]),as.matrix(ipar["d"]),as.matrix(ipar["c"]))

aa<-matrix(ipar[,1:p],ncol=p)
dd<-ipar[,p+1]
cc<-ipar[,p+2]

if (!easiness) dd<- -dd

if (!is.null(seed)) set.seed(seed)
TH<-matrix(rnorm(n.simulee*p),n.simulee,p)
random<-matrix(runif(n.simulee*ni),n.simulee,ni)

if (all(sigma!=0) && p>1) TH<-TH%*%chol(sigma) #Cholesky decomposition

resp<-matrix(0,n.simulee,ni)
pp<-matrix(NA,n.simulee,ni)

for (i in 1:ni) {
Z<-matrix(0,nrow=n.simulee,ncol=1)
for (h in 1:p) {
Z<-Z+aa[i,h]*TH[,h]
}
pp[,i]<-cc[i] + (1-cc[i])/(1+exp(-D*(Z + dd[i])))

resp[,i]<-resp[,i]+ifelse(random[,i]<pp[,i],1,0)
}

resp<-as.data.frame(resp)
names(resp)<-paste("R",1:ni,sep="")

out<-list(call=call,theta=TH,resp=resp)
return(out)
}

