`binomRRci` <-
function(x, ...){UseMethod("binomRRci")}


# # # THis is the extension of Gart-Nams (1988) "Crude" interval to the situation
# of multiple contrasts

# With the following changes:

# 1)
# If all x going to the denominator AND numerator of the contrast are ==0
# --> [lower, upper]<-[0, Inf]

# 2)
# If all x going to the numerator of the contrast are == 0
# --> [lower,]<-0

# 3)
# If all x going to the denominator of the contrast are == 0
# --> [,upper]<-Inf





`binomRRci.default` <-
function(x, n, names=NULL, type="Dunnett",
 cmat=NULL, alternative="two.sided", conf.level=0.95, dist="MVN", ...)
{

aargs<-list(...)

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


#### RESTRICT BOUNDS function

restrictbounds<-function(x, cmat, conf.int)
{

if(all(x!=0)) {return(conf.int)}
else{

M<-nrow(cmat)

warning("0 occured in the data and the risk ratio might not be defined")

for(i in 1:M)
{

# all 0
if( all( x[sign(cmat[i,])!=0]==0 ) )
 {conf.int[i,1]<-0; conf.int[i,2]<-Inf}

# all numerator 0

if(all ( x[sign(cmat[i,])==1]==0 ) )
 {conf.int[i,1]<-0}


# all denominator 0

if(all (x[sign(cmat[i,])==(-1)]==0 ) )
 {conf.int[i,2]<-Inf}

}

return(conf.int)
}

}




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
    }


# get the point and variance estimates

    est <- binomest.default(x=x, n=n, names=gnames, success=aargs$success)

  #  print(est)

if(any(x<4) & any(n-x<4))
 {warning("Normal approximation might be inappropriate")}


    xI<-est$Y+0.5
    nI<-est$n+0.5

    estI<-log(xI/nI)
    varI<- (1/xI)-(1/nI)

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
    cnamesRR<-unlist(lapply(cnamesSlist, FUN=function(x){paste(x, collapse="/")}))
    }
    else{cnamesRR<-paste("C",1:nrow(cmat),sep="")}

    rownames(estimate)<-cnamesRR
    rownames(conf.int)<-cnamesRR    

    # expand the CI to the estimate:
    conf.int<-restrictbounds(x=x, cmat=cmat, conf.int=conf.int)

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
    out$method <-"Crude normal approximation for the risk ratio"

    class(out) <- c("binomRRci", "sci")
    return(out)
}


`binomRRci.formula` <-
function(formula, data, type="Dunnett",
 cmat=NULL, alternative="two.sided", conf.level=0.95, dist="MVN", ...)
{
aargs<-list(...)

est<-binomest.formula(formula=formula, data=data, success=aargs$success)

aargs$x<-est$Y
aargs$n<-est$n
aargs$names<-est$names
aargs$type<-type
aargs$cmat<-cmat
aargs$alternative<-alternative
aargs$conf.level<-conf.level
aargs$success<-est$success
aargs$dist<-dist

out<-do.call("binomRRci.default", aargs)

return(out)

}




`binomRRci.table` <-
function(x, type="Dunnett",
 cmat=NULL, alternative="two.sided", conf.level=0.95, dist="MVN", ...)
{
aargs<-list(...)

est<-binomest.table(x=x, success=aargs$success)

aargs$x<-est$Y
aargs$n<-est$n
aargs$names<-est$names
aargs$type<-type
aargs$cmat<-cmat
aargs$alternative<-alternative
aargs$conf.level<-conf.level
aargs$success<-est$success
aargs$dist<-dist

out<-do.call("binomRRci.default", aargs)

return(out)

}

`binomRRci.matrix` <-
function(x, type="Dunnett",
 cmat=NULL, alternative="two.sided", conf.level=0.95, dist="MVN", ...)
{
aargs<-list(...)

tab<-as.table(x)

est<-binomest.table(x=tab, success=aargs$success)

aargs$x<-est$Y
aargs$n<-est$n
aargs$names<-est$names
aargs$type<-type
aargs$cmat<-cmat
aargs$alternative<-alternative
aargs$conf.level<-conf.level
aargs$success<-est$success
aargs$dist<-dist

out<-do.call("binomRRci.default", aargs)

return(out)

}






