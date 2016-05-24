
'binomRDci' <- function(x,...){UseMethod("binomRDci")}

'binomRDci.default' <- function(x, n, names=NULL, type="Dunnett",
 cmat=NULL, method="Wald", alternative="two.sided", conf.level=0.95, dist="MVN", ...)
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

    est <- binomest.default(x=x, n=n, names=gnames, method = method, success=aargs$success)

    out <- Waldci(estp = est$estp,
     varcor=est$varcor,
     varp = est$varp,
     cmat = cmat, 
     alternative = alternative,
     conf.level=conf.level,
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


    class(out) <- c("binomRDci", "sci")
    return(out)
}




'binomRDci.table' <- function(x, type="Dunnett",
 cmat=NULL, method="Wald", alternative="two.sided", conf.level=0.95, dist="MVN",...)
{

 aargs<-list(...)

 est<-binomest.table(x, method=method, success=aargs$success)

 n<-est$n
 k<-length(n)

   if(any(n<2)) 
     {stop("Sample size must be at least 2 in each group")}   

    gnames<-est$names

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

    out <- Waldci(estp = est$estp,
     varcor=est$varp,
     varp = est$varcor,
     cmat = cmat, 
     alternative = alternative,
     conf.level=conf.level,
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

    class(out) <- "binomRDci"
    return(out)
}


'binomRDci.matrix' <- function(x, type="Dunnett",
 cmat=NULL, method="Wald", alternative="two.sided", conf.level=0.95, dist="MVN",...)
{

 mat<-as.table(x)

 out<-binomRDci.table(mat,
 type=type, cmat=cmat,
 method=method, alternative=alternative, conf.level=conf.level, dist=dist, ...)
    return(out)
}



'binomRDci.formula' <- function(formula, data, type="Dunnett",
 cmat=NULL, method="Wald", alternative="two.sided", conf.level=0.95, dist="MVN",...)
{

 aargs<-list(...)
 aargs$formula<-formula
 aargs$data<-data
 aargs$method<-method

 est<-do.call("binomest.formula", args=aargs)

 out<-binomRDci.default(x=est$Y, n=est$n, names=est$names,
 type=type, cmat=cmat,
 method=method, alternative=alternative, conf.level=conf.level, success=est$success, dist=dist)
    return(out)

}




