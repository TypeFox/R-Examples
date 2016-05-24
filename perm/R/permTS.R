`permTS` <-
function (x, ...){
    UseMethod("permTS")
}

`permTS.formula` <-
function(formula, data, subset, na.action, ...){
    ## mostly copied from wilcox.test.formula
    if (missing(formula) || (length(formula) != 3) || (length(attr(terms(formula[-2]), 
        "term.labels")) != 1)) 
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    groupname<-names(mf)[2]

    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if (nlevels(g) != 2) 
        stop("grouping factor must have exactly 2 levels")
    DATA <- split(mf[[response]], g)
    names(DATA) <- c("x", "y")
    out <- do.call("permTS", c(DATA, list(...)))
    out$data.name <- DNAME
    glevels<-levels(g)
    
    names(out$null.value)<-names(out$estimate)<-paste("mean ",groupname,"=",glevels[1]," - mean ",groupname,"=",glevels[2],sep="")
    out
}
`permTS.default` <-
function(x, y, alternative = c("two.sided", "less", "greater"), 
    exact = NULL, method=NULL, methodRule=methodRuleTS1, control=permControl(),...){
    
    cm<-control$cm
    nmc<-control$nmc
    seed<-control$seed
    digits<-control$digits
    p.conf.level<-control$p.conf.level
    setSEED<-control$setSEED
    tsmethod<-control$tsmethod


    ## program originally written so that tsmethod="abs" was 
    ## alternative="two.sidedAbs", changed to avoid confusion with coin package
    if (alternative[1]=="two.sidedAbs"){
        warning("alternative='two.sidedAbs' may be deprecated in the future,
            use alternative='two.sided' and control=permControl(tsmethod='abs'))")
        alternative<-"two.sided"
        tsmethod<-"abs"
    }
    alternative <- match.arg(alternative)

    if (!(tsmethod=="central" | tsmethod=="abs") & alternative=="two.sided"){
        stop("only tsmethod='central' and tsmethod='abs' allowed")
    }
    if (tsmethod=="abs" & alternative=="two.sided"){
        alternative<-"two.sidedAbs"
    } else if (tsmethod=="central" & alternative=="two.sided"){
        alternative<-"two.sided"
    }		
																																																																																																																																																																																																																																																																																																																																																																																																																										
    if (!is.numeric(x) | !is.numeric(y) | !is.vector(x) | !is.vector(y) ) stop("x and y must be numeric vectors")

    W<-c(x,y)
    Z<-c(rep(1,length(x)),rep(0,length(y)))
    
    if (is.null(method))    method<-methodRule(W,Z,exact)
    method.OK<-(method=="pclt" | method=="exact.mc" | method=="exact.network" | method=="exact.ce")
    if (!method.OK) stop("method not one of: 'pclt', 'exact.mc'. 'exact.network', 'exact.ce'")
    mout<-switch(method,
        pclt=twosample.pclt(W,Z),
        exact.network=twosample.exact.network(W,Z,digits),
        exact.ce=twosample.exact.ce(W,Z,cm,digits),
        exact.mc=twosample.exact.mc(W,Z,alternative,nmc,seed,digits,p.conf.level,setSEED))
    p.values<-mout$p.values
    PVAL <- switch(alternative,two.sided=p.values["p.twosided"],greater=p.values["p.gte"],
        less=p.values["p.lte"],two.sidedAbs=p.values["p.twosidedAbs"]) 
    if (method=="exact.network") METHOD<-"Exact Permutation Test (network algorithm)"
    else if (method=="pclt") METHOD<-"Permutation Test using Asymptotic Approximation"
    else if (method=="exact.mc") METHOD<-"Exact Permutation Test Estimated by Monte Carlo"
    else if (method=="exact.ce") METHOD<-"Exact Permutation Test (complete enumeration)"
    m<-match.call()
    xname<-deparse(substitute(x))
    yname<-deparse(substitute(y))
    #xname<-as.character(m$x)
    #yname<-as.character(m$y)
    if (length(xname)>1 || nchar(xname)>10) xname<-c("GROUP 1")
    if (length(yname)>1 || nchar(yname)>10) yname<-c("GROUP 2")
    DNAME <- paste(xname, "and", yname)
    Z<-mout$Z
    if (!is.null(Z)) names(Z)<-"Z"
   
    null.value<-0
    estimate<-mean(x)-mean(y)
    names(estimate)<-names(null.value)<-paste("mean",xname,"- mean",yname)
    p.conf.int<-mout$p.conf.int
    if (method!="exact.mc") nmc<-NULL 
    OUT <- list(statistic = Z, estimate=estimate, parameter = NULL, p.value = as.numeric(PVAL), 
        null.value = null.value, alternative = alternative, method = METHOD, 
        data.name = DNAME, p.values=p.values,p.conf.int=p.conf.int, nmc=nmc)
    if (method=="exact.mc") class(OUT) <- "mchtest"
    else class(OUT) <- "htest"
    return(OUT)
}



