`Simpson` <-
function(p)
{
1-sum(p^2)
}

`estSimpson` <-
function(x)
{
n<-sum(x)
estp<-x/n
Simpson(estp)*n/(n-1)
}

`estSimpsonf` <-
function(X, f)
{

if(any(floor(X)!=X))
 {warning("The elements of X should be integers!")}

J<-length(f)

if(nrow(X)!=J)
 {stop("The number of columns in X must be equal to the length of f!")}

varestSimpson<-function(x)
{
n<-sum(x)
estp<-x/n
S2<-sum(estp^2)
S3<-sum(estp^3)

2*(S2 + 2*(n-2)*S3 + (3-2*n)*S2^2)/(n*(n-1))
}

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

estlist<-lapply(X=Xcs, FUN=function(x){estSimpson(x)})
estv<-unlist(estlist)
 
varestlist<-lapply(X=Xcs, FUN=function(x){varestSimpson(x)})
varestv<-unlist(varestlist)

return(list(estimate=estv, varest=varestv, table=tab))

}

`Simpsonci` <-
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

esti<-estSimpsonf(X=X, f=as.factor(f))

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
out$methodname <- "Confidence intervals for differences of Simpsons Indices"
out$sample.estimates <- esti


class(out)<-c("Simpsonci", "sci")
return(out)

}



`print.Simpsonci` <-
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
"for differences of Simpson indices", "\n")
}
else{
cat("Local", round(x$conf.level*100,3),"percent-confidence intervals","\n",
"for differences of Simpson indices", "\n")
}

conf.int <- cbind(x$estimate, x$conf.int)

colnames(conf.int)[1] <- "estimate"

print(round(conf.int, digits=digits))

invisible(x)
}


`summary.Simpsonci` <-
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
"Estimated Simpson index"=round(object$sample.estimates$estimate, digits),
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



