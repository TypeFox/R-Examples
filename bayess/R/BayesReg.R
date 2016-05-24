BayesReg=function(y,X,g=length(y),betatilde=rep(0,dim(X)[2]),prt=TRUE)
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
U=solve(t(X)%*%X)%*%t(X)
alphaml=mean(y)
betaml=U%*%y
s2=t(y-alphaml-X%*%betaml)%*%(y-alphaml-X%*%betaml)
kappa=as.numeric(s2+t(betatilde-betaml)%*%t(X)%*%X%*%(betatilde-betaml)/(g+1))
malphabayes=alphaml
mbetabayes=g/(g+1)*(betaml+betatilde/g)
msigma2bayes=kappa/(n-3)
valphabayes=kappa/(n*(n-3))
vbetabayes=diag(kappa*g/((g+1)*(n-3))*solve(t(X)%*%X))
vsigma2bayes=2*kappa^2/((n-3)*(n-4))
postmean=c(malphabayes,mbetabayes)
postsqrt=sqrt(c(valphabayes,vbetabayes))
intlike=(g+1)^(-p/2)*kappa^(-(n-1)/2)

bayesfactor=rep(0,p)
if (p>=2)
{
for (i in 1:p)
{
p0=p-1
X0=X[,-i]
U0=solve(t(X0)%*%X0)%*%t(X0)
betatilde0=U0%*%X%*%betatilde
betaml0=U0%*%y
s20=t(y-alphaml-X0%*%betaml0)%*%(y-alphaml-X0%*%betaml0)
kappa0=as.numeric(s20+t(betatilde0-betaml0)%*%t(X0)%*%X0%*%(betatilde0-betaml0)/(g+1))
intlike0=(g+1)^(-p0/2)*kappa0^(-(n-1)/2)
bayesfactor[i]=intlike/intlike0
}
}
if (p==1)
{
intlike0=(t(y-alphaml)%*%(y-alphaml))^(-(n-1)/2)
bayesfactor=intlike/intlike0
}

bayesfactor=log10(bayesfactor)
evid=rep("    ",p+1)
for (i in 1:p)
{
if (bayesfactor[i]<0) evid[i+1]="      "
if (0<=bayesfactor[i] & bayesfactor[i]<=0.5) evid[i+1]="   (*)"
if (0.5<bayesfactor[i] & bayesfactor[i]<=1) evid[i+1]="  (**)"
if (1<bayesfactor[i] & bayesfactor[i]<=2) evid[i+1]=" (***)"
if (bayesfactor[i]>2) evid[i+1]="(****)"
}

if (prt==TRUE)
{
vnames="Intercept"
for (i in 1:p) vnames=c(vnames,paste("x",i,sep=""))
cat("\n")
print(data.frame(PostMean=round(postmean,4),PostStError=round(postsqrt,4),
Log10bf=c("",round(bayesfactor,4)),EvidAgaH0=evid,row.names=vnames))
cat("\n")
cat("\n")
cat(paste("Posterior Mean of ","Sigma2",":"," ",round(msigma2bayes,4),sep=""))
cat("\n")
cat(paste("Posterior StError of ","Sigma2",":"," ",round(sqrt(vsigma2bayes),4),sep=""))
cat("\n")
cat("\n")
}
list(postmeancoeff=postmean,postsqrtcoeff=postsqrt,log10bf=bayesfactor,postmeansigma2=msigma2bayes,
postvarsigma2=vsigma2bayes)
}
