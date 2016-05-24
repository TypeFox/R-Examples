`permKS` <-
function (x, ...){
    UseMethod("permKS")
}

`permKS.formula` <-
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
    resp <- mf[[response]]
    out <- do.call("permKS", c(list(x=resp,g=g), list(...)))
    out$data.name <- DNAME
    out
}

`permKS.default` <-
function(x, g, exact = NULL, method=NULL, methodRule=methodRuleKS1,
    control=permControl(),...){
   
    cm<-control$cm
    nmc<-control$nmc
    seed<-control$seed
    digits<-control$digits
    p.conf.level<-control$p.conf.level
    setSEED<-control$setSEED
 
    if (!is.numeric(x) | (!is.character(g) & !is.factor(g))) stop("x must be numeric and g must be character or factor vectors")

    if (is.null(method))    method<-methodRule(x,g,exact)

    method.OK<-(method=="pclt" | method=="exact.mc")
    if (!method.OK) stop("method not one of: 'pclt', 'exact.mc'")
    mout<-switch(method,
        pclt=ksample.pclt(x,g),
        exact.mc=ksample.exact.mc(x,g,nmc,seed,digits,p.conf.level,setSEED))
    p.values<-mout$p.values
    PVAL<-p.values["p.twosided"]
    if (method=="pclt") METHOD<-"K-Sample Asymptotic Permutation Test"
    else if (method=="exact.mc") METHOD<-"K-Sample Exact Permutation Test Estimated by Monte Carlo"
    xname<-deparse(substitute(x))
    gname<-deparse(substitute(g))
    if (length(xname)>1) xname<-c("x")
    if (length(gname)>1) gname<-c("g")
   
     DNAME <- paste(xname, "and", gname)
   
    chisq<-mout$chisq.value
    if (!is.null(chisq)) names(chisq)<-"Chi Square"
    df<-mout$df
    if (!is.null(df)) names(df)<-"df"
    if (method!="exact.mc") nmc<-NULL 
    OUT <- list(statistic = chisq, parameter=df, estimate=NULL, 
        p.value = as.numeric(PVAL), method = METHOD, 
        data.name = DNAME, p.conf.int=mout$p.conf.int, nmc=nmc)
    if (method=="exact.mc"){ class(OUT) <- "mchtest"
    } else class(OUT) <- "htest"
    return(OUT)
}


