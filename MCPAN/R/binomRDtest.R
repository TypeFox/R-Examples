'binomRDtest' <- function(x,...){UseMethod("binomRDtest")}


'binomRDtest.default' <- function(x, n, names=NULL, type="Dunnett",
 cmat=NULL, method="Wald", alternative="two.sided", dist="MVN", ...)
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


# include checks for validity of methods


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


    est <- binomest.default(x=x, n=n, names=gnames, method = method, success=aargs$success)

    out <- Waldtest(estp = est$estp,
     varp = est$varcor,
     cmat = cmat, 
     alternative = alternative,
     dist=dist)

    out$estimate <- cmat %*% est$estimate
    colnames(cmat) <- gnames
    out$x <- x
    out$n <- n
    out$p <- est$estimate
    out$success <- est$success
    out$names <- gnames
    out$method <- method
    out$cmat <- cmat
    class(out) <- "binomRDtest"
    return(out)
}


'binomRDtest.table' <- function(x, type="Dunnett",
 cmat=NULL, method="Wald", alternative="two.sided", dist="MVN", ...)
{

aargs<-list(...)


 est<-binomest.table(x, method=method, success=aargs$success)
 n<-est$n
 k<-length(n)

    if(any(n<2)) 
     {stop("Sample size must be at least 2 in each group")}   

    gnames<-est$names

# check the contrast matrix

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

    out <- Waldtest(estp = est$estp,
     varp = est$varcor,
     cmat = cmat, 
     alternative = alternative,
     dist=dist)

    out$estimate <- cmat %*% est$estimate
    colnames(cmat) <- gnames
    out$x <- est$Y
    out$n <- est$n
    out$p <- est$estimate
    out$success <- est$success
    out$names <- gnames
    out$method <- method
    out$cmat <- cmat

    class(out) <- "binomRDtest"
    return(out)
}




'binomRDtest.matrix' <- function(x, type="Dunnett",
 cmat=NULL, method="Wald", alternative="two.sided", dist="MVN",...)
{

 mat<-as.table(x)

 aargs<-list(...)

 aargs$x<-mat
 aargs$type<-type
 aargs$cmat<-cmat
 aargs$method<-method
 aargs$alternative<-alternative
 aargs$dist<-dist 

 out<-do.call("binomRDtest.table", args=aargs)
 
 return(out)
}




'binomRDtest.formula' <- function(formula, data, type="Dunnett",
 cmat=NULL, method="Wald", alternative="two.sided", dist="MVN", ...)
{

 testaargs<-aargs<-list(...)
 aargs$formula<-formula
 aargs$data<-data
 aargs$method<-method

 est<-do.call("binomest.formula", args=aargs)
 testaargs$x<-est$Y
 testaargs$n<-est$n
 testaargs$names<-est$names
 testaargs$type<-type
 testaargs$cmat<-cmat
 testaargs$method<-method
 testaargs$alternative<-alternative
 testaargs$dist<-dist 
 testaargs$success<-est$success

 out<-do.call("binomRDtest.default", args=testaargs)
 return(out)

}




