nbfactors <-
function(data,x=1,test=x[1],pvalues=NULL,maxnbfactors=8,diagnostic.plot=FALSE,min.err=1e-03){
dig = 2
if (class(data)[1]!="FAMTdata") stop("Class of data should be FAMTdata")
nbcovar = ncol(data$covariates)
if (any(!is.element(x,1:nbcovar))) stop(paste("x should be a subset of 1:",nbcovar,sep=""))
if (!is.element(test,1:nbcovar)) stop(paste("test should be one of ",x,sep=""))
m = nrow(data$expression)
n = ncol(data$expression)
falist = vector(length=maxnbfactors+1,"list")
falist[[1]] = list(B=matrix(0,ncol=1,nrow=m))
falist[-1] = lapply(1:maxnbfactors,emfa,data=data,x=x,test=test,pvalues=pvalues,min.err=min.err)
Blist = lapply(falist,function(fa,m) matrix(fa$B,nrow=m),m=m)

bivprob = function (rho,lower,upper=-lower,mean=0) {
nu = 0
low = rep(as.double((lower - mean)),2)
upp = rep(as.double((upper - mean)),2)
if (any(lower == upper)) return(0)
infin = c(2, 2)
infin = as.integer(infin)
low = replace(low, low == -Inf, 0)
upp = replace(upp, upp == Inf, 0)
rho = as.double(rho)
prob = as.double(0)
a = lapply(rho,function(r,low,upp) biv.nt.prob(df=Inf,lower=low,upper=upp,mean=rep(0,2),S=matrix(c(1,r,r,1),2,2)),
       low=low,upp=upp)
return(unlist(a))
}

Dt = function(rho) {
threshold=0.05
ut = qnorm(1 - threshold/2)
delta = unlist(lapply(rho,bivprob,lower=-ut)) - (1 - threshold)^2
dt = delta/(threshold * (1 - threshold))
return(dt)
}

VarInflation = function(data,Blist,x,test,maxnbfactors,pvalues,dig) {
n = ncol(data$expression)
m = nrow(data$expression)
vecrho = round(seq(10^(-dig),1,10^(-dig)),digits=dig)
vecdt = unlist(lapply(vecrho,Dt))
rdata = residualsFAMT(data,x,test,pvalues)
SelectH0 = rdata$SelectH0
m0 = length(SelectH0)
rdata = rdata$residuals[,SelectH0]
B0 = lapply(Blist,function(b,s,m) matrix(b[s,],nrow=m0),s=SelectH0,m=m0)
sampled = sample(1:m0,min(1000,m0))
sampsize = length(sampled)
cdata = scale(rdata[,sampled])
cordata = t(cdata)%*%cdata/(n-1)
sdt = rep(0,maxnbfactors+1)
names(sdt) = paste(0:maxnbfactors,"factors")
for (i in 1:(maxnbfactors+1)) {
   print(paste("Calculating criterion for the model with",i-1,"factors"))
   B = matrix(B0[[i]][sampled,],nrow=sampsize)
   sdb = sqrt(1-apply(B^2,1,sum))
   matrho = cordata - B%*%t(B)
   matrho = sweep(matrho,2,FUN="/",STATS=sdb)
   matrho = sweep(matrho,1,FUN="/",STATS=sdb)
   rho = matrho[col(matrho)>row(matrho)]
   rho[abs(rho)>=1] = 1
   veccor = sort(round(abs(rho),digits=dig))
   duplic = duplicated(veccor)
   vduplic = sort(unique(veccor[duplic]))
   vunic = setdiff(unique(veccor),vduplic)
   dtunic = vecdt[is.element(vecrho,vunic)]
   dtduplic = vecdt[is.element(vecrho,vduplic)]
   vmatch = match(vecrho,veccor,0)
   nboccur = diff(c(vmatch[vmatch>0],length(veccor)+1))
   nboccur = nboccur[nboccur>1]
   sdt[i] = 2*(m0-1)*(sum(dtunic)+crossprod(nboccur,dtduplic))/(sampsize*(sampsize-1))  }
return(sdt) }

sdt = VarInflation(data,Blist,x,test,maxnbfactors,pvalues,dig)
if (diagnostic.plot) {
   dev.new()
   plot(0:maxnbfactors,sdt,ylab="Variance Inflation Criterion",xlab="Number of factors",bty="l",
      lwd=1.25,type="b",pch=16,cex.lab=1.25,cex=1.25,cex.axis=1.25)
}
if (which.min(sdt)==1) opt = 0
if (which.min(sdt)>1) { 
   jumps = -diff(sdt)/sdt[-length(sdt)]
   opt = max((1:maxnbfactors)[jumps>0.05]) }
list(criterion=sdt,optimalnbfactors=opt)
}
