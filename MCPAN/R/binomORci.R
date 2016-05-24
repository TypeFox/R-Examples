`binomORci` <-
function(x, ...){UseMethod("binomORci")}


`binomORci.default` <-
function(x, n, names=NULL, type="Dunnett", method="GLM",
 cmat=NULL, alternative="two.sided", conf.level=0.95, dist="MVN", ...)
{

aargs<-list(...)

method<-match.arg(method, choices=c("GLM", "Woolf"))

# check x, n

    k <- length(x)
    if(k<2)
     {stop("x, n must contain at least two elements")}

    if(any(n<2)) 
     {stop("sample size must be at least 2 in each group")}   

    if(length(n)!=k)
     {stop("x, n must be of equal length")}

    if(!is.numeric(x) | !is.numeric(n))
     {stop("x, n must be numeric vectors")}



# check the names
# if no names are given, available names of x are taken
# if names are specified, names of x/n are overwritten

    gnames<-names

    if(is.null(gnames))
     {
      if(is.null(names(x)))
       {gnames <- names(x) <- names(n) <- as.character(1:k)}
       else
        {gnames <- names(n) <- names(x)}
     }
     else
      {
       if(length(gnames)!=k)
        {stop("length of names does not fit to length of x,n")}
        else
         {names(x) <- names(n) <- as.character(gnames)}
      }

# check the contrast matrix

    if (is.null(cmat)) {
      if(type=="Dunnett") {
        if(is.null(aargs$base)){base<-1}
        else{base<-aargs$base}
        cmat <- contrMat(n=n, type=type, base=base)
       }
       else{cmat <- contrMat(n = n, type = type)}
    }
    else {
        if (!is.matrix(cmat) || ncol(cmat) != k)
         {stop("cmat must be a matrix with number of columns = number of groups")}
	attr(cmat, which="type")<-"User defined"
    }


# get the point and variance estimates

switch(method,
GLM={
    est <- binomest.default(x=x, n=n, names=gnames, success=aargs$success)
    grp<-factor(est$names, levels=est$names)
    logisticfit <- glm(cbind(x,n-x) ~ 0 + grp, family=binomial(link="logit"))
    eta <- coefficients(logisticfit)
    sigma <- vcov(logisticfit)

    out<-Waldci( cmat=cmat,
     estp=eta,
     varp=diag(sigma),
     varcor=diag(sigma),
     alternative = alternative,
     conf.level=conf.level,
     dist=dist)


    conf.int <- exp(out$conf.int)
    estimate <- exp(cmat %*% eta)

    if(!is.null(dimnames(cmat)[[1]]))
    {
    cnamesD<-dimnames(cmat)[[1]]
    cnamesSlist<-strsplit(cnamesD, "-")
    cnamesOR<-unlist(lapply(cnamesSlist, FUN=function(x){paste(x, collapse="/")}))
    }
    else{cnamesOR<-paste("C",1:nrow(cmat),sep="")}

    rownames(estimate)<-cnamesOR
    rownames(conf.int)<-cnamesOR    

    out$conf.int<-conf.int
    out$estimate<-estimate
    colnames(cmat) <- gnames
    out$x <- x
    out$n <- n
    out$p <- est$estimate
    out$success <- est$success
    out$names <- gnames
    out$method <- "glm"
    out$cmat <- cmat
},
Woolf={

# Restrict Bounds to certain values, when x=0 or x=n occurs

restrictboundsOR<-function(x, n, cmat, conf.int)
{

if(all(x!=0) & all(x!=n) )
 {return(conf.int)}
else{

M<-nrow(cmat)

warning("0 occured in the data and the risk ratio might not be defined")

for(i in 1:M)
{

# 'Uninformative events'
# all x=0 in the ith contrast:

cNOT0<-sign(cmat[i,])!=0



if( all( x[cNOT0]==0 ) )
 {
  conf.int[i,1]<-0
  conf.int[i,2]<-Inf
#  cat("Contrast",i,"All x=0 \n")
 }

# all x=n in the ith contrast:

if( all( x[cNOT0]==n[cNOT0] ) )
 {
  conf.int[i,1]<-0
  conf.int[i,2]<-Inf
# cat("Contrast",i,"All x=n \n")
 }


# all numerator x=0 in contrast i

cINDN<-sign(cmat[i,])==1

if(all ( x[cINDN]==0 ) )
 {
 conf.int[i,1]<-0
# cat("Numerator Contrast",i,"All x=0 \n")
 }


# all numerator x=n in contrast i

if(all ( x[cINDN]==n[cINDN] ) )
 {
  conf.int[i,2] <- Inf
 # cat("Numerator Contrast",i,"All x=n \n")
 }


# all denominator x=0

cINDD<-sign(cmat[i,])==(-1)

if(all (x[cINDD]==0 ) )
 {
  conf.int[i,2]<-Inf
 # cat("Denominator Contrast",i,": all x=0 \n")
 }

# all denominator x=n

if(all (x[cINDD]==n[cINDD] ) )
 {
  conf.int[i,1] <- 0 
 # cat("Denominator Contrast",i,": all x=n \n")
 }


}

return(conf.int)
}

}


    est <- binomest.default(x=x, n=n, names=gnames, success=aargs$success)
    
    XPI <- est$Y+0.5 
    XQI <- est$n-est$Y+0.5

    estI <- log(XPI/XQI)

    varI <- 1/XPI + 1/XQI

    out<-Waldci( cmat=cmat,
     estp=estI,
     varp=varI,
     varcor=varI,
     alternative = alternative,
     conf.level=conf.level,
     dist=dist)


    conf.int <- exp(out$conf.int)
    estimate <- exp(cmat %*% estI)

    if(!is.null(dimnames(cmat)[[1]]))
    {
    cnamesD<-dimnames(cmat)[[1]]
    cnamesSlist<-strsplit(cnamesD, "-")
    cnamesOR<-unlist(lapply(cnamesSlist, FUN=function(x){paste(x, collapse="/")}))
    }
    else{cnamesOR<-paste("C",1:nrow(cmat),sep="")}

    rownames(estimate)<-cnamesOR
    rownames(conf.int)<-cnamesOR    
 
    conf.int<-restrictboundsOR(x=x, n=n, cmat=cmat, conf.int=conf.int)

    out$conf.int<-conf.int
    out$estimate<-estimate
    colnames(cmat) <- gnames
    out$x <- x
    out$n <- n
    out$p <- est$estimate
    out$success <- est$success
    out$names <- gnames
    out$method <- "Adjusted Woolf"
    out$cmat <- cmat
})




    class(out) <- c("binomORci", "sci")
    return(out)
}


`binomORci.formula` <-
function(formula, data, type="Dunnett", method="GLM",
 cmat=NULL, alternative="two.sided", conf.level=0.95, dist="MVN", ...)
{
aargs<-list(...)

est<-binomest.formula(formula=formula, data=data, success=aargs$success)

aargs$x<-est$Y
aargs$n<-est$n
aargs$names<-est$names
aargs$type<-type
aargs$method<-method
aargs$cmat<-cmat
aargs$alternative<-alternative
aargs$conf.level<-conf.level
aargs$success<-est$success
aargs$dist<-dist

out<-do.call("binomORci.default", args=aargs)

return(out)

}




`binomORci.table` <-
function(x, type="Dunnett", method="GLM",
 cmat=NULL, alternative="two.sided", conf.level=0.95, dist="MVN", ...)
{
aargs<-list(...)

est<-binomest.table(x=x, success=aargs$success)

aargs$x<-est$Y
aargs$n<-est$n
aargs$names<-est$names
aargs$type<-type
aargs$method<-method
aargs$cmat<-cmat
aargs$alternative<-alternative
aargs$conf.level<-conf.level
aargs$success<-est$success
aargs$dist<-dist

out<-do.call("binomORci.default", args=aargs)

return(out)

}

`binomORci.matrix` <-
function(x, type="Dunnett", method="GLM",
 cmat=NULL, alternative="two.sided", conf.level=0.95, dist="MVN", ...)
{
aargs<-list(...)

tab<-as.table(x)

est<-binomest.table(x=tab, success=aargs$success)

aargs$x<-est$Y
aargs$n<-est$n
aargs$names<-est$names
aargs$type<-type
aargs$method<-method
aargs$cmat<-cmat
aargs$alternative<-alternative
aargs$conf.level<-conf.level
aargs$success<-est$success
aargs$dist<-dist

out<-do.call("binomORci.default", args=aargs)

return(out)

}

