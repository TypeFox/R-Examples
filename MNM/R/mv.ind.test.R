######################
###
#
# Test for independence based in different score function
#
###
######################


mv.ind.test <- function(X, Y, score="identity", method = "approximation", n.simu = 1000 ,na.action=na.fail)
    {
    dname.x <- deparse(substitute(X))
    dname.y <- deparse(substitute(Y))

    DNAME <- paste(dname.x, " and ", dname.y, sep = "")

    score <- match.arg(score,c("identity","sign","rank","symm"))
    method <- match.arg(method,c("approximation","permutation"))
    
    n <- dim(X)[1]
    if(n != dim(Y)[1]) stop("'X' and 'Y' must have the same number of rows")
    
    X <- na.action(X)
        if (!is.null(attr(X, "na.action"))) 
            Y <- Y[-(attr(X, "na.action")), ]
    Y <- na.action(Y)
        if (!is.null(attr(Y, "na.action"))) 
            X <- X[-(attr(Y, "na.action")), ]
    
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    X<-as.matrix(X)
    
    if(!all(sapply(Y, is.numeric))) stop("'Y' must be numeric")
    Y<-as.matrix(Y)
    
    p.X <- dim(X)[2]
    p.Y <- dim(Y)[2]
    p <- p.X + p.Y
    
    res1<-switch(score,
        "identity"={
               iind(X = X, Y = Y, method = method, n.simu = n.simu, n = n, p.X = p.X, p.Y = p.Y)
               }
        ,
        "sign"={
                ssind.inner(X, Y, method =  method, n.simu = n.simu, n = n, p.X = p.X, p.Y = p.Y, p = p)
                }
        ,
        "rank"={
                srind.inner(X, Y, method =  method, n.simu = n.simu, n = n, p.X = p.X, p.Y = p.Y, p = p)
                }
        ,
        "symm"={
                symmsind.inner(X, Y, method =  method, n.simu = n.simu, n = n, p.X = p.X, p.Y = p.Y, p = p)
                }
        )
        
     ALTERNATIVE <- "two.sided"
     NVAL <- 0
     names(NVAL) <- "measure of dependence"

    
    res<-c(res1,list(data.name=DNAME,alternative=ALTERNATIVE,null.value=NVAL))
    class(res) <- "htest"    
    return(res)
    
    }
