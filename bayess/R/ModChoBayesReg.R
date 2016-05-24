ModChoBayesReg=function(y,X,g=length(y),betatilde=rep(0,dim(X)[2]),niter=100000,prt=TRUE)
{

X=as.matrix(X)
n=length(y)
p=dim(X)[2]

for (i in 1:p)
{
X[,i]=X[,i]-mean(X[,i])
X[,i]=X[,i]/sqrt(mean(X[,i]^2))
}

if (det(t(X)%*%X)<=1e-7) stop("The design matrix has a rank lower than the number of explanatory variables!
Calculations cannot be done and the process should be stopped!",call.=FALSE)
alphaml=mean(y)
intlike0=(t(y-alphaml)%*%(y-alphaml))^(-(n-1)/2)

if (p<=15)
{
intlike=rep(0,2^p)
intlike[1]=intlike0
for (i in 2:2^p)
{
gam=as.integer(intToBits(i-1)[1:p]==1)
pgam=sum(gam)
Xgam=X[,which(gam==1)]
Ugam=solve(t(Xgam)%*%Xgam)%*%t(Xgam)
betatildegam=Ugam%*%X%*%betatilde
betamlgam=Ugam%*%y
s2gam=t(y-alphaml-Xgam%*%betamlgam)%*%(y-alphaml-Xgam%*%betamlgam)
kappagam=as.numeric(s2gam+t(betatildegam-betamlgam)%*%t(Xgam)%*%Xgam%*%(betatildegam-betamlgam)/(g+1))
intlike[i]=(g+1)^(-pgam/2)*kappagam^(-(n-1)/2)
}
intlike=intlike/sum(intlike)
modcho=order(intlike)[2^p:(2^p-9)]
probtop10=intlike[modcho]
modtop10=rep("",10)
for (i in 1:10)
{
modtop10[i]=paste(which(intToBits(modcho[i]-1)==1),collapse=" ")
}

if (prt==TRUE)
{
cat("\n")
cat("Number of variables less than 15")
cat("\n")
cat("Model posterior probabilities are calculated exactly")
cat("\n")
cat("\n")
print(data.frame(Top10Models=modtop10,PostProb=round(probtop10,4)))
cat("\n")
cat("\n")
}
list(top10models=modtop10,postprobtop10=probtop10)
}

else
{
gamma=rep(0,niter)
mcur=sample(c(0,1),p,replace=TRUE)
gamma[1]=sum(2^(0:(p-1))*mcur)+1
pcur=sum(mcur)
if (pcur==0) intlikecur=intlike0 else
{
Xcur=X[,which(mcur==1)]
Ucur=solve(t(Xcur)%*%Xcur)%*%t(Xcur)
betatildecur=Ucur%*%X%*%betatilde
betamlcur=Ucur%*%y
s2cur=t(y-alphaml-Xcur%*%betamlcur)%*%(y-alphaml-Xcur%*%betamlcur)
kappacur=as.numeric(s2cur+t(betatildecur-betamlcur)%*%t(Xcur)%*%Xcur%*%(betatildecur-betamlcur)/(g+1))
intlikecur=(g+1)^(-pcur/2)*kappacur^(-(n-1)/2)
}
for (i in 1:(niter-1))
{
mprop=mcur
j=sample(1:p,1)
mprop[j]=abs(mcur[j]-1)
pprop=sum(mprop)
if (pprop==0) intlikeprop=intlike0 else
{
Xprop=X[,which(mprop==1)]
Uprop=solve(t(Xprop)%*%Xprop)%*%t(Xprop)
betatildeprop=Uprop%*%X%*%betatilde
betamlprop=Uprop%*%y
s2prop=t(y-alphaml-Xprop%*%betamlprop)%*%(y-alphaml-Xprop%*%betamlprop)
kappaprop=as.numeric(s2prop+t(betatildeprop-betamlprop)%*%t(Xprop)%*%Xprop%*%(betatildeprop-betamlprop)/(g+1))
intlikeprop=(g+1)^(-pprop/2)*kappaprop^(-(n-1)/2)
}
if (runif(1)<=(intlikeprop/intlikecur))
{
mcur=mprop
intlikecur=intlikeprop
}
gamma[i+1]=sum(2^(0:(p-1))*mcur)+1
}
gamma=gamma[20001:niter]
res=as.data.frame(table(as.factor(gamma)))
odo=order(res$Freq)[length(res$Freq):(length(res$Freq)-9)]
modcho=res$Var1[odo]
probtop10=res$Freq[odo]/(niter-20000)
modtop10=rep("",10)
for (i in 1:10)
{
modtop10[i]=paste(which(intToBits(as.integer(paste(modcho[i]))-1)==1),collapse=" ")
}
if (prt==TRUE)
{
cat("\n")
cat("Number of variables greather than 15")
cat("\n")
cat("Model posterior probabilities are estimated by using an MCMC algorithm")
cat("\n")
cat("\n")
print(data.frame(Top10Models=modtop10,PostProb=round(probtop10,4)))
cat("\n")
cat("\n")
}
list(top10models=modtop10,postprobtop10=probtop10)
}
}
