`Shannon` <-
function(p)
{
-sum(p*log(p))
}

`estShannon` <-
function(x, Nspec=NULL)
{

if(is.null(Nspec))
 {Nspec<-length(x)}

Nind<-sum(x)

estp<-as.numeric(x)/Nind

estpo0 <- estp[estp>0]

pvec <- matrix(estpo0, ncol=1)

m1 <- diag(estpo0) - (pvec %*% t(pvec))

sigma2 <- (t(log(pvec)) %*% m1) %*% log(pvec)

Nspecpres <- length(estpo0)

estH <- (-1)*sum(estpo0*log(estpo0))

estHBC <- estH + (Nspec-1)/(2*Nind)

estHBCo0 <- estH + (Nspecpres-1)/(2*Nind)

 varest<-as.numeric(sigma2)/Nind

out<-list( estimate=estHBC,
 estraw=estH,
 estHBCo0=estHBCo0,
 varest=varest)

return(out)

}

`estShannonf` <-
function(X, f)
{

X<-as.data.frame(X)

J<-length(f)

if(nrow(X)!=J)
 {stop("The number of columns in X must be equal to the length of f!")}

if(any(floor(X)!=X))
 {warning("The elements of X should be integers!")}

X<-as.data.frame(X)

names<-colnames(X)

ff<-as.factor(f)
Xs<-split(X, f=ff)

Xcs<-lapply(X=Xs, FUN=function(x){apply(X=x, MARGIN=2, FUN=sum)})

tab<-matrix(ncol=ncol(X), nrow=length(Xcs))
for(i in seq(along.with=Xcs))
{
tab[i,]<-Xcs[[i]]
}

colnames(tab)<-names
rownames(tab)<-names(Xcs)


ngroups<-nrow(tab)

estimate<-numeric(length=ngroups)
estraw<-numeric(length=ngroups)
estHBCo0<-numeric(length=ngroups)
varest<-numeric(length=ngroups)

for(i in 1:ngroups)
{
temp<-estShannon(tab[i,])

estimate[i]<-temp$estimate
estraw[i]<-temp$estraw
estHBCo0[i]<-temp$estHBCo0
varest[i]<-temp$varest

}

gnames<-rownames(tab)
names(estimate)<-names(estraw)<-names(estHBCo0)<-names(varest)<-gnames

out<-list(estimate=estimate,
estraw=estraw,
estHBCo0=estHBCo0,
varest=varest,
table=tab
)

return(out)
}

`Shannonci` <-
function(X, f, cmat=NULL, type="Dunnett", alternative = "two.sided", conf.level = 0.95, dist = "MVN", ...)
{

aargs<-list(...)

type<-match.arg(type, choices=c("Dunnett","Tukey","Sequen"))
alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

if(length(conf.level)!=1)
 {stop("argument conf.level should be a single numeric value")}

if(conf.level<=0.5|conf.level>=1) 
 {stop("argument conf.level should be a single numeric value between 0.5 and 1")}

if(!is.numeric(conf.level)) 
 {stop("argument conf.level should be a single numeric value")}

if(!is.data.frame(X))
 {stop("X must be of an object of class 'data.frame'")}

if(!is.factor(f))
 {f<-as.factor(f)}

k<-length(levels(f))

if(k<=1)
 {stop("The factor variable f should have at least 2 levels to be compared")}

esti<-estShannonf(X=X, f=as.factor(f))

n<-apply(X=esti$table, MARGIN=1, FUN=sum)

if(is.null(cmat))
 {
  if(is.null(aargs$base)){base<-1}
   else{base<-aargs$base}
 cmat<-contrMat(n=n, type=type, base=base)}
else
 {
 if(ncol(cmat)!=k)
  {stop("Number of columns in cmat should be the same as the number of levels in f")}
 }


out <- Waldci(cmat=cmat, estp=esti$estimate, varp=esti$varest,
 varcor=esti$varest, alternative = alternative, 
    conf.level = conf.level, dist = dist) 

TAB<-esti$table

estimate <- cmat %*% matrix(esti$estimate,ncol=1)

out$estimate <- estimate
out$cmat <-cmat
out$methodname <- "Confidence intervals for differences of Shannon Indices"
out$sample.estimates <- esti

class(out)<-c("Shannonci", "sci")
return(out)

}


`print.Shannonci` <-
function(x,...)
{

# A table of confidence intervals

aargs<-list(...)

if(is.null(aargs$digits))
 {digits<-4}
else
 {digits<-aargs$digits}

dist<-attr(x$quantile, which="dist")

if(dist=="MVN"){
cat("Simultaneous", round(x$conf.level*100,3),"percent-confidence intervals","\n",
"for differences of Shannon indices", "\n")
}
else{
cat("Local", round(x$conf.level*100,3),"percent-confidence intervals","\n",
"for differences of Shannon indices", "\n")
}

conf.int <- cbind(x$estimate, x$conf.int)

colnames(conf.int)[1] <- "estimate"

print(round(conf.int, digits=digits))

invisible(x)
}


`summary.Shannonci` <-
function(object,...)
{

aargs<-list(...)
if(is.null(aargs$digits))
 {digits<-4}
else
 {digits<-aargs$digits}

cat("\n Data: \n")

print(object$sample.estimates$table)

cat("\n Summary statistics: \n")

N<-apply(X=object$sample.estimates$table, MARGIN=1, FUN=function(x){sum(x)})

summarystat<-rbind(
"Total number of individuals"=round(N,0),
"Shannon index, bias corrected estimate"=round(object$sample.estimates$estimate, digits),
"Shannon index, raw estimate"=round(object$sample.estimates$estraw, digits),
"Variance estimate"=round(object$sample.estimates$varest,digits)
)

colnames(summarystat)<-names(object$sample.estimate$estimate)

print(summarystat, digits=digits)

cat("\n Contrast matrix: \n")

print(object$cmat, digits=digits)

cat("\n")

print(object, digits=digits)

invisible(object)
}


