
###########


checkargssim<-function(coef, vcov, cmat)
{

if(!is.matrix(vcov))
 {stop("'vcov' must be a matrix")}

if(!is.numeric(vcov))
 {stop("'vcov' must be a numeric matrix")}

if(!is.matrix(cmat))
 {stop("'cmat' must be a matrix")}

if(!is.numeric(cmat))
 {stop("'cmat' must be a numeric matrix")}

if(!is.numeric(coef))
 {stop("'coef' must be a numeric vector")}


k<-length(coef)

if(k<=2)
 {warning("Applying this function makes hardly sense for only two parameters!")}

if(any(dim(vcov)!=k))
 {stop("'vcov' must be a k-times-k matrix if 'coef' has length k")}

if(ncol(cmat)!=k)
 {stop("'cmat' must be a matrix with k columns, if 'coef' has length k")}

}



#################

simplesimint<-function(coef, vcov, cmat, df=NULL, conf.level = 0.95,
 alternative = c("two.sided", "less", "greater") )
{

checkargssim(coef=coef, vcov=vcov, cmat=cmat)

alternative<-match.arg(alternative)

if(length(conf.level)!=1)
 {stop("'conf.level' must be a single numeric value")}

if(!is.numeric(conf.level))
 {stop("'conf.level' must be a single numeric value")}

if(conf.level<=0.5 | conf.level>=1)
 {stop("'conf.level' must be a single numeric value between 0.5 and 1")}

estC <- cmat %*% coef

vcovC <- cmat %*% vcov %*% t(cmat)

corrC <- cov2cor(vcovC)

stderr <- sqrt(diag(vcovC))

M<-nrow(cmat)

if(is.null(df))
{

switch(alternative,

two.sided={ quanti <- qmvnorm(p=conf.level, corr=corrC, tail="both.tails")$quantile
   lCI <- estC-quanti*stderr
   uCI <- estC+quanti*stderr
 },
less={ quanti <- qmvnorm(p=conf.level, corr=corrC, tail="lower.tail")$quantile
     lCI <- rep(-Inf, M)
     uCI <- estC+quanti*stderr
   },
greater={quanti <- qmvnorm(p=conf.level, corr=corrC, tail="upper.tail")$quantile
      lCI <- estC+quanti*stderr
      uCI <- rep(Inf, M)
     })

attr(quanti, which="dist") <- paste(M, "-variate normal distribution", collapse="")

}
else{

if(length(df)!=1)
 {stop("'df' must be a single number")}

if(!is.numeric(df))
 {stop("'df' must be numeric!")}

if(df<2)
 {warning("degree of freedom less than 2 leads to computational problems!")}


switch(alternative,

two.sided={ quanti <- qmvt(p=conf.level, corr=corrC, tail="both.tails", df=df)$quantile
   lCI <- estC-quanti*stderr
   uCI <- estC+quanti*stderr
 },
less={ quanti <- qmvt(p=conf.level, corr=corrC, tail="lower.tail", df=df)$quantile
     lCI <- rep(-Inf, M)
     uCI <- estC+quanti*stderr
   },
greater={quanti <- qmvt(p=conf.level, corr=corrC, tail="upper.tail", df=df)$quantile
      lCI <- estC+quanti*stderr
      uCI <- rep(Inf, M)
     })

attr(quanti, which="dist")<-paste(M, "-variate t-distribution", collapse="")

}

out<-list(
estimate=estC,
lower=lCI,
upper=uCI,
cmat=cmat,
alternative=alternative,
conf.level=conf.level,
quantile=quanti,
df=df,
stderr=stderr,
vcovC=vcovC)

class(out)<-"simplesimint"

return(out)

}


print.simplesimint <- function(x,...)
{
pargs<-list(...)

LEVEL<-round(x$conf.level*100,2)

CI<-data.frame(x$estimate, x$lower, x$upper)
colnames(CI)<-c("Estimate","Lower","Upper")

cat("\n Simultaneous", LEVEL, "% confidence intervals: \n")
pargs$x<-CI

do.call("print", pargs)
invisible(x)
}

##########

summary.simplesimint <- function(object,...)
{

pargs<-list(...)

LEVEL<-round(object$conf.level*100,2)

CI<-data.frame(object$estimate, object$lower, object$upper)
colnames(CI)<-c("Estimate","Lower","Upper")

CMAT<-object$cmat

QUANT<-round(object$quantile, 4)
DIST<-attr(object$quantile, which="dist")
DF<-object$df


VCOV<-object$vcovC

cat("\n Simultaneous ", LEVEL, "% confidence intervals: \n", sep="")
pargs$x<-CI

do.call("print", pargs)


cat("\n Used quantile: ", QUANT, ",\n", sep="")

cat("obtained from a ",DIST, "\n", sep="")

if(!is.null(DF)){cat( "with ", DF, " degrees of freedom.\n", sep="")}


#

cat("\n Used contrast matrix: \n")
pargs$x<-CMAT

do.call("print", pargs)

#

cat("\n Resulting variance covariance matrix of the contrasts: \n")
pargs$x<-VCOV

do.call("print", pargs)

# 

cat("\n Corresponding correlation matrix of the contrasts: \n")
pargs$x<-cov2cor(VCOV)

do.call("print", pargs)

invisible(object)

}


##########

plotCI.simplesimint<-function(x, ...)
{
pargs<-list(...)

estimate<-as.numeric(x$estimate)
lower<-as.numeric(x$lower)
upper<-as.numeric(x$upper)

names(estimate)<-rownames(x$cmat)

pargs$estimate <- estimate
pargs$lower<-lower
pargs$upper<-upper
pargs$alternative<-x$alternative

do.call("plotCII", pargs)

}





####




