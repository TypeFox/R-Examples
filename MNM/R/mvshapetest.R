mv.shape.test <- function(X,score="identity", location="est", na.action=na.fail, ...)
    {
    DNAME <- deparse(substitute(X))
    
    score <- match.arg(score,c("identity", "sign", "rank", "symmsign"))
    #stand <- match.arg(estimate,c("inner", "outer"))
    location <- match.arg(location,c("est", "origin"))
    
    X<-na.action(X)
    
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    X<-as.matrix(X)
    n<-dim(X)[1]
    p<-dim(X)[2]
    
    #if(stand=="inner") stop('Sorry, only estimate="outer" is possible')
    if(score=="rank") stop('Sorry, score="rank" is not implemented')
    
    res1<-switch(score,
        "identity"={
               mauchly.test(X,location,n,p)
               }
        ,
        "sign"={
               SCov.test(X,location, n, p)
               }
        ,
        "symmsign"={
               SSCov.test(X, n, p)
               }
        )
    res<-c(res1,list(data.name=DNAME))
    class(res) <- "htest" 
    return(res)
    }
    
 
