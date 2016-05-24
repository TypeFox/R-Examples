##########################################################################
# math functions
##########################################################################


#binomial coefficient
binom.coef=function(n,m) { prod((n-m+1):n)/prod(1:m) }

expit=function (x) {1/(1+exp(-x))}
logit=function (x) {log(x/(1-x))}


as.binary <- function(n, base=2 , r=FALSE) {
   out <- NULL
   while(n > 0) {
     if(r) {
       out <- c(out , n%%base)
     } else {
       out <- c(n%%base , out)
     }   
     n <- n %/% base
   }
   return(out)
}



##-- Stirling numbers of the 2nd kind
##-- (Abramowitz/Stegun: 24,1,4 (p. 824-5 ; Table 24.4, p.835)

##> S^{(m)}_n = number of ways of partitioning a set of $n$ elements into $m$
##> non-empty subsets

Stirling2 <- function(n,m)
{
    ## Purpose:  Stirling Numbers of the 2-nd kind
    ##      S^{(m)}_n = number of ways of partitioning a set of
    ##                      $n$ elements into $m$ non-empty subsets
    ## Author: Martin Maechler, Date:  May 28 1992, 23:42
    ## ----------------------------------------------------------------
    ## Abramowitz/Stegun: 24,1,4 (p. 824-5 ; Table 24.4, p.835)
    ## Closed Form : p.824 "C."
    ## ----------------------------------------------------------------
    ## maechler@_stat.math.ethz.ch
    
    if (0 > m || m > n) stop("'m' must be in 0..n !")     
    k <- 0:m
    sig <- rep(c(1,-1)*(-1)^m, length= m+1)
    # 1 for m=0; -1 1 (m=1)     
    ## The following gives rounding errors for (25,5) :     
    ## r <- sum( sig * k^n /(gamma(k+1)*gamma(m+1-k)) )     
    ga <- gamma(k+1)
    round(sum( sig * k^n /(ga * rev(ga)))) 
}

logSumExp=function (logx){
    logMeanExp(logx, 1)
}
logSumExpFor2=function (logx, logy){
    c=max(logx, logy)
    dif=abs(logx-logy)
    if (dif>300) return (c)
    else {
        log(sum(exp(logx-c), exp(logy-c)))+c
    }
}
# log( exp(logx1)-exp(logx2) )
logDiffExp=function (logx1,logx2){
    c=logx1
    c+ log(1-exp(logx2-logx1))
}

logMeanExp=function (logx,B=NULL){
# mean function for small numbers
# logx is a vector of large negative values
# return log (sum(exp(logx))/B)
# calculate log of the mean of a vector which contains 0 and very small real numbers (logged)
# return log of the mean
    if (is.null(B)) B=length(logx)
    c=max(logx)
    log(sum(exp(logx-c))/B)+c
}
# logMeanExp(log(1:5), 5) # test, should return log(3)

logDiffExp=function (logx1, logx2){
# diff function for small numbers
# logx1 and logx2 are typically large negative values, logx1>logx2
# return log (exp(logx1)-exp(logx2))
    if (logx1<logx2) {cat("\nlWarning [logDiffExp]: first argument smaller than second, return NaN.\n\n"); return (NaN);}
    c=max(logx1, logx2)
    log (exp(logx1-c)-exp(logx2-c))+c
}
# logDiffExp(log(2), log(1)) # test, should return 0

# from combinat package
permn=function (x, fun = NULL, ...) 
{
    if (is.numeric(x) && length(x) == 1 && x > 0 && trunc(x) == 
        x) 
        x <- seq(x)
    n <- length(x)
    nofun <- is.null(fun)
    out <- vector("list", gamma(n + 1))
    p <- ip <- seqn <- 1:n
    d <- rep(-1, n)
    d[1] <- 0
    m <- n + 1
    p <- c(m, p, m)
    i <- 1
    use <- -c(1, n + 2)
    while (m != 1) {
        out[[i]] <- if (nofun) 
            x[p[use]]
        else fun(x[p[use]], ...)
        i <- i + 1
        m <- n
        chk <- (p[ip + d + 1] > seqn)
        m <- max(seqn[!chk])
        if (m < n) 
            d[(m + 1):n] <- -d[(m + 1):n]
        index1 <- ip[m] + 1
        index2 <- p[index1] <- p[index1 + d[m]]
        p[index1 + d[m]] <- m
        tmp <- ip[index2]
        ip[index2] <- ip[m]
        ip[m] <- tmp
    }
    out
}



##########################################################################
# array and matrix functions
##########################################################################


getMidPoints=function(x){
    ((c(0,x)+c(x,0))/2)[2:length(x)] 
}
#getMidPoints(1:10)


# the input is a list, typically the output from a sapply call that should be matrix, but have different length
fill.jagged.array=function(a) {
    # don't check is.matrix, because for some reason, it will return true
    max.len=max(sapply(a, length))
    sapply(a, function (e) {
        c(e, rep(NA, max.len-length(e)))
    })    
}


last = function (x, n=1, ...) {
    if (is.vector(x)) x[length(x)]
    else if (is.array(x)) x[length(x)]
    else if (is.list(x)) x[[length(x)]]
    else if (is.character(x)) tail (readLines(x), n=n, ...) # read file
    else stop ("last(): x not supported")
}

# return a subset of data that is 1 row every thin.factor rows
ThinRows = function (dat, thin.factor=10) {
    NumRows = nrow(dat)
    dat[1:(NumRows/thin.factor)*thin.factor,]
}
thin.rows=ThinRows

#mix two arrays in an interlacing way
mix = function (a, b) {
    if (length(a)!=length(b)) print ("Length of two arguments to mix function not equal.")
    out = rep (a, each=2)
    for (i in 1:length(a)) {
        out[2*i]=b[i]
    }
    out
}


# like lag, move vector to the right/left by given number of steps
# x is a vector
shift.right = function (x, k=1) {
    p=length(x)
    x[(1+k):p]=x[1:(p-k)]
    x[1:k]=NA
    x
}

shift.left = function (x, k=1) {
    p=length(x)
    x[1:(p-k)] = x[(1+k):p]
    x[(1+p-k):p]=NA
    x
}



# trace
tr=function(m) sum(diag(m))


# serial covariance matrix
AR1 = function (p, w) {
    m = matrix(1, p, p)
    for (i in 1:p) {
        for (j in 1:p) {
            m [i,j]=w**abs(i-j)
        }
    }
    m
}

# exchangeable covariance matrix
EXCH = function (p, rho) {
    m = matrix(1, p, p)
    for (i in 1:p) {
        for (j in 1:p) {
            if (i!=j) m [i,j]=rho
        }
    }
    m
}

getUpperRight = function (matri, func=NULL) {
    n=nrow (matri)
    out= numeric ( (n-1)*n/2 )
    index=0
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            index=index+1
            out[index]=matri[i,j]
        }
    }
    if (is.null(func)) {
        out
    } else {
        func(out)
    }
}


#repeat a matrix in a block diagonal fashion
rep.matrix.block = function (x, times=2, ...) {
    orig.dim = nrow (x)     
    m = matrix (0, orig.dim * times, orig.dim * times)
    for (i in 1: times) {
        m[(1+(i-1)*orig.dim):(i*orig.dim), (1+(i-1)*orig.dim):(i*orig.dim)]=x
    }
    m
}


#it does not work on data.frame
rep.matrix = function (x, times=1, each=1, by.row=TRUE, ...) {
    if (by.row) {
        colnames.=colnames(x)
        new.matrix=matrix(0, nrow(x)*each*times, ncol(x) )
        for (i in 1:nrow(x)) {
            for (j in 1:each) {
                new.matrix[(i-1)*each + j,] = x[i,]
            }
        }
        
        if(times>1) {
            for (i in 2:times) {
                new.matrix[ ((i-1)*nrow(x)*each+1) : (i*nrow(x)*each), ] = new.matrix[1:(nrow(x)*each), ]
            }
        }
    
        dimnames(new.matrix)[[2]]=colnames.
        new.matrix
    }
    else {
        t ( rep.matrix(t(x), times, each, by.row=TRUE) )
    }    
}

# rep.data.frame(chi[21,], 2)
rep.data.frame = function (x, times=1, ...){
    out = x
    if (time==1) return (out)
    for (i in 2:times) {
        out = rbind (out, x)
    }
    out
}



##########################################################################
# string functions
##########################################################################


#### file name manipulation

getFileStem=function(file.name){
    substr(file.name, 1, lastIndex(file.name, ".")-1 )
}
getStem=getFileStem
fileStem=getFileStem

getExt = function(file.name){
    substr(file.name, lastIndex(file.name, ".")+1, nchar(file.name) )
}

# s="prefix_name" is changed to "name"
remove.prefix=function(s,sep="_"){
    tmp=strsplit(s,sep)    
    sapply(tmp, function (x) {
        if (length(x)==1) return (x)
        concatList(x[-1],sep)
    })
}



#### misc

escapeUnderline=function (name) {
    gsub("_", "\\_", name)
}



#### string search

firstIndex =function (s1, s2) {
    k=nchar (s2)
    ret=-1;
    for (i in 1:(nchar(s1)-k+1) ) {
        if (substr(s1, i, i+k-1)==s2) {
            ret=i;
            break;
        }
    }
    ret
}

lastIndex =function (s1, s2) {
    k=nchar (s2)
    ret=-1;
    for (i in 1:(nchar(s1)-k+1) ) {
        if (substr(s1, i, i+k-1)==s2) {
            ret=i;
        }
    }
    ret
}

# return TRUE if s1 starts with s2? 
startsWith=function(s1, s2){
    sapply (s1, function (s) {
        if ( substring (s, 1, nchar(s2)) == s2 ) {
            return (TRUE);
        } else {
            return (FALSE);
        }
    })
}

# return TRUE if s1 ends with s2, s1 can be a vector
endsWith=function(s1, s2){
    sapply (s1, function (s) {
        if ( substring (s, nchar(s)-nchar(s2)+1, nchar(s)) == s2 ) {
            return (TRUE);
        } else {
            return (FALSE);
        }
    })
}



# return TRUE if s1 contains s2
contain =function (s1, s2) {
    sapply (s1, function (s) {
        k=nchar (s2)
        matched=FALSE
        for (i in 1:(nchar(s)-k+1) ) {
            if (substr(s, i, i+k-1)==s2) {
                matched=TRUE
            }
        }
        matched
    })
}



# paste two strings together
# e.g. "a" %+% "b"

"%+%" <- function (a, b) {
    out=paste(a,b,sep="")
    out
}

concatList = function (lis, sep=""){
    out=lis[[1]]
    i=2
    while (i<=length(lis)){
        out=out%+%sep%+%lis[[i]]
        i=i+1
    }
    out
}



myprint <- function(object, ...) UseMethod("myprint") 

# this function is placed at the bottom of the file because it contains "\""), which makes all the following line be miss-interpreted as being in quotes
myprint.default = function (..., newline=TRUE, digits=3) {   
    digits.save=getOption("digits")
    options(digits=digits)
    object <- as.list(substitute(list(...)))[-1]
    x=list(...)
    for (i in 1:length(x)) {
        if (is(x[[i]],"formula")) {cat(as.character(x[[i]])); next}
        tmpname <- deparse(object[[i]])[1]
        #str(tmpname)
        #str(gsub("\\\\","\\",gsub("\"", "", tmpname)))
        #str(x[[i]])
        #if (gsub("\\\\","\\",gsub("\"", "", tmpname))!=x[[i]]) {
        if (contain(tmpname, "\"") | contain(tmpname, "\\")) {
            for (a in x[[i]]) cat(a," ")
        } else {
            cat (tmpname %+% " = ")
            for (a in x[[i]]) cat(a," ")
            if (i!=length(x)) cat ("; ")
        }
    }
    if (newline)  cat("\n")
    options(digits=digits.save)
}
#a="hello"; b="you"; myprint (a, b); myprint ("test"); myprint.default ("\t")
