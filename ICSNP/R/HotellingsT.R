`HotellingsT2`<-function(X,...)
    {
    UseMethod("HotellingsT2")
    }

`HotellingsT2.default` <-
function(X,Y=NULL,mu=NULL,test="f",na.action=na.fail,...)
    {
    if (is.null(Y)) 
       {
       DNAME<-deparse(substitute(X))
       }
    else
       {
       DNAME=paste(deparse(substitute(X)),"and",deparse(substitute(Y)))
       }
       
    X<-na.action(X)
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    X<-as.matrix(X)
    
    p<-dim(X)[2]
    
    if (!is.null(Y))
        {
        Y<-na.action(Y)
       if(!all(sapply(Y, is.numeric))) stop("'Y' must be numeric")
       if (p!=dim(Y)[2]) stop("'X' and 'Y' must have the same number of columns")
       Y<-as.matrix(Y)
        }

    if (is.null(mu)) mu<-rep(0,p) 
    else if (length(mu)!=p) stop("length of 'mu' must equal the number of columns of 'X'")

    test<-match.arg(test,c("f","chi"))
    
    if (is.null(Y) & test=="f") version<-"one.sample.f"
    if (is.null(Y) & test=="chi") version<-"one.sample.chi"
    if (!is.null(Y) & test=="f") version<-"two.sample.f"
    if (!is.null(Y) & test=="chi") version<-"two.sample.chi"
    
    res1<-switch(version,
            "one.sample.f"={
                result<-HotellingsT.internal(X,mu=mu,test=test)
                STATISTIC<-result$test.statistic
                names(STATISTIC)<-"T.2"
                PVAL<-result$p.value
                METHOD<-"Hotelling's one sample T2-test"
                PARAMETER<-c(result$df.1,result$df.2)
                names(PARAMETER)<-c("df1","df2")
                RVAL<-list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)
               
               RVAL}
            ,
            "one.sample.chi"={
                result<-HotellingsT.internal(X,mu=mu,test=test)
                STATISTIC<-result$test.statistic
                names(STATISTIC)<-"T.2"
                PVAL<-result$p.value
                METHOD<-"Hotelling's one sample T2-test"
                PARAMETER<-c(result$df.1)
                names(PARAMETER)<-c("df")
                RVAL<-list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)
               
               RVAL}
            ,
            "two.sample.f"={
                result<-HotellingsT.internal(X,Y,mu,test)
                STATISTIC<-result$test.statistic
                names(STATISTIC)<-"T.2"
                PVAL<-result$p.value
                METHOD<-"Hotelling's two sample T2-test"
                PARAMETER<-c(result$df.1,result$df.2)
                names(PARAMETER)<-c("df1","df2")
                RVAL<-list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)
               
               RVAL}
            ,
            "two.sample.chi"={
                result<-HotellingsT.internal(X,Y,mu,test)
                STATISTIC<-result$test.statistic
                names(STATISTIC)<-"T.2"
                PVAL<-result$p.value
                METHOD<-"Hotelling's two sample T2-test"
                PARAMETER<-c(result$df.1)
                names(PARAMETER)<-c("df")
                RVAL<-list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)
               
               RVAL}
            )
    ALTERNATIVE="two.sided"
    NVAL<-paste("c(",paste(mu,collapse=","),")",sep="")
    if (is.null(Y)) names(NVAL)<-"location" else names(NVAL)<-"location difference"
    res<-c(res1,list(data.name=DNAME,alternative=ALTERNATIVE,null.value=NVAL))
    class(res)<-"htest"
    return(res)
    }

`HotellingsT2.formula` <-
function(formula, na.action = na.fail, ...)
   {
    if(missing(formula)
       || (length(formula) != 3)
       || (length(attr(terms(formula[-2]), "term.labels")) != 1))
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if(nlevels(g) != 2)
        stop("grouping factor must have exactly 2 levels")
    DATA <- split(as.data.frame(mf[[response]]), g)
    names(DATA) <- c("X", "Y")
    y <- do.call("HotellingsT2", c(DATA, list(...)))
    y$data.name <- DNAME
    return(y)
  }
