`permTREND` <-
function (x, ...){
    UseMethod("permTREND")
}

`permTREND.formula` <-
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
    DNAME <- paste(names(mf), collapse = " and ")
    if (nchar(DNAME)>15) DNAME<-"x and y"
    groupname<-names(mf)[2]

    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- mf[[-response]]
    resp<-mf[[response]]
    out <- do.call("permTREND", c(list(x=resp,y=g), list(...)))
    out$data.name <- DNAME
    names(out$estimate)<-names(out$null.value)<-paste("correlation of",DNAME)

    out
}

`permTREND.default` <-
function(x, y, alternative = c("two.sided", "less", "greater"), 
    exact = NULL, method=NULL, methodRule=methodRuleTREND1, 
    control=permControl(),...){

    cm<-control$cm
    nmc<-control$nmc
    seed<-control$seed
    digits<-control$digits
    p.conf.level<-control$p.conf.level
    setSEED<-control$setSEED
    tsmethod<-control$tsmethod    

    if (!is.numeric(x) | !is.numeric(y) | !is.vector(x) | !is.vector(y) ) stop("x and y must be numeric vectors")

    if (length(x)!=length(y)) stop("x and y must be the same length")
 
    if (is.null(method))    method<-methodRule(x,y,exact)

    method.OK<-(method=="pclt" | method=="exact.mc")
    if (!method.OK) stop("method not one of: 'pclt', 'exact.mc'")



    ## program originally written so that tsmethod="abs" was 
    ## alternative="two.sidedAbs", changed to avoid confusion with coin package
    ## change it back internally
    if (alternative[1]=="two.sidedAbs"){
        warning("alternative='two.sidedAbs' may be deprecated in future versions,
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

    mout<-switch(method,
        pclt=trend.pclt(x,y),
        exact.mc=trend.exact.mc(x,y,alternative,nmc,seed,digits,p.conf.level,setSEED))
    p.values<-mout$p.values
    PVAL <- switch(alternative,
        two.sided=p.values["p.twosided"],
        greater=p.values["p.gte"],
        less=p.values["p.lte"],two.sidedAbs=p.values["p.twosidedAbs"]) 
    if (method=="pclt") METHOD<-"Permutation Test using Asymptotic Approximation"
    else if (method=="exact.mc") METHOD<-"Exact Permutation Test Estimated by Monte Carlo"
    m<-match.call()
    #xname<-deparse(substitute(x))
    #yname<-deparse(substitute(y))
    xname<-as.character(m$x)
    yname<-as.character(m$y)

    if (length(xname)>1 || nchar(xname)>10) xname<-c("RESPONSE")
    if (length(yname)>1 || nchar(yname)>10) yname<-c("COVARIATE")
    DNAME <- paste(xname, "and", yname)
    Z<-mout$Z
    if (!is.null(Z)) names(Z)<-"Z"
 
    null.value<-0
    estimate<-cor(x,y)
    names(estimate)<-names(null.value)<-paste("correlation of ",xname," and ",yname,sep="")
    p.conf.int<-mout$p.conf.int
    if (method!="exact.mc") nmc<-NULL 
    OUT <- list(statistic = Z, estimate=estimate, parameter = NULL, p.value = as.numeric(PVAL), 
        null.value = null.value, alternative = alternative, method = METHOD, 
        data.name = DNAME, p.values=p.values, p.conf.int=p.conf.int, nmc=nmc)
    if (method=="exact.mc") class(OUT) <- "mchtest"
    else class(OUT) <- "htest"
     return(OUT)
}


