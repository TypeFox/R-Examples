`ojaCsampleTest`<-function(X,...)
    {
    UseMethod("ojaCsampleTest")
    }

`ojaCsampleTest.default` <-
function(X,Y,mu=NULL, scores="sign", p = 1, method = "approximation", n.simu = 1000, center = "ojaMedian", na.action=na.fail, ...)
    {
    DNAME<-paste(deparse(substitute(X)),"and",deparse(substitute(Y)))
    
    X<-na.action(X)
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    X<-as.matrix(X)   
    k<-dim(X)[2]
   
    if (is.null(mu)) mu<-rep(0,k)
    else if (length(mu)!=k) stop("length of 'mu' must equal the number of columns of 'X'")

    Y<-na.action(Y)
    if(!all(sapply(Y, is.numeric))) stop("'Y' must be numeric")
    Y<-as.matrix(Y)
    if (k!=dim(Y)[2]) stop("'X' and 'Y' must have the same number of columns")
    if (dim(X)[1]<2 | dim(Y)[1]<2) stop("both 'X' and 'Y' must have at least two observations")

    
    X<-sweep(X,2,mu)
    
    scores<-match.arg(scores,c("sign","rank"))
    method<-match.arg(method,c("approximation","permutation"))
    
    res1<-switch(scores,
        "sign"={
               oja2sampleSignTest(X=X,Y=Y,p=p,method=method,n.simu=n.simu, center=center,...)
               }
        ,
        "rank"={
               oja2sampleRankTest(X=X,Y=Y,p=p,method=method,n.simu=n.simu,...)
               }
        )
    
    ALTERNATIVE="two.sided"
    NVAL<-paste("c(",paste(mu,collapse=","),")",sep="")
    names(NVAL)<-"location difference"
    res<-c(res1,list(data.name=DNAME,alternative=ALTERNATIVE,null.value=NVAL))
    class(res)<-"htest"
    return(res)
    }

`ojaCsampleTest.formula` <-
function(formula, scores="sign", p = 1, method = "approximation", n.simu = 1000, center = "ojaMedian", data, subset, na.action,...)
   {
   if(missing(formula)
       || (length(formula) != 3)
       || (length(attr(terms(formula[-2]), "term.labels")) != 1))
        stop("'formula' missing or incorrect")
    
    mf <- match.call(expand.dots = FALSE)
    
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]                # Retain only the named arguments
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    
    DNAME <- paste(names(mf), collapse = " by ")
    
    ## Get the data matrices
    Y <- model.response(mf, "numeric")
    if (is.vector(Y)) stop("response must be at least bivariate")
    
    g <- factor(mf[[2]])
    n.g<-nlevels(g)
    if(n.g < 2)
        stop("grouping factor must have two or more levels")
    {    
    
    scores<-match.arg(scores,c("sign","rank"))
    method<-match.arg(method,c("approximation","permutation"))
    
    if(n.g==2)
        {
        DATA <- split(as.data.frame(Y), g)
        names(DATA) <- c("X", "Y")
        RVAL <- do.call("ojaCsampleTest", c(DATA, scores=scores, p = p, method = method, n.simu = n.simu, center = center, list(...)))
        RVAL$data.name <- DNAME
        }   
    else
        { 
        if ("mu" %in% names(list(...)) ) {if(any( list(...)$mu!=0 )) {stop("if there are more than two groups 'mu' should not be specified or 0")}}
        DATA <- as.matrix(Y)
        names(DATA)<-"X"
        names(g)<-g
        y <- do.call("ojaCsampleTestFormula",list(DATA,g, scores = scores, p = p, method = method, n.simu = n.simu, center = center,...))
        
        
        METHOD<- switch(scores,"sign"="OJA C SAMPLE SIGN TEST",
                               "rank"="OJA C SAMPLE RANK TEST"
                               )
        PARAMETER<-y$df
        RVAL<-c(y,list(method=METHOD, data.name=DNAME))
        class(RVAL)<-"htest"           
        
        }
        return(RVAL)
    }
  }
