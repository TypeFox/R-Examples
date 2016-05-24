##' @export
categorical <- function(x,formula,K,beta,p,liability=FALSE,regr.only=FALSE,...) {

    if (is.character(formula)) {
        regr <- FALSE
        X <- formula
    } else {
        y <- getoutcome(formula)
        X <- attributes(y)$x
        regr <- TRUE
        if (length(y)==0) regr <- FALSE
        if (length(attributes(y)$x)==0) {
            X <- y; regr <- FALSE
        }
    }
    if (!missing(p)) {
        if (!missing(K) && K!=length(p)+1) stop("Wrong dimension of 'p'")
        if (is.numeric(p) && sum(p)>1) warning("'p' sum > 1")
        if (is.logical(all.equal(1.0,sum(p)))) p <- p[-length(p)]
    }
    if (missing(K)) {
        if (!is.null(list(...)$labels)) K <- length(list(...)$labels)
        if (!missing(beta)) K <- length(beta)
        if (!missing(p)) K <- length(p)+1
    }
    if (!regr.only) {
        if (missing(p)) p <- rep(1/K,K-1)
        pname <- names(p)
        if (is.null(pname)) pname <- rep(NA,K-1)
        ordinal(x,K=K,liability=liability,p=p,constrain=pname,...) <- X
        if (!regr) return(x)
    }

    if (missing(beta)) beta <- rep(0,K)
    fname <- paste(gsub(" ","",deparse(formula)),seq(K)-1,sep=":")
    fpar <- names(beta)
    if (is.null(fpar)) fpar <- fname

    parameter(x,fpar,start=beta) <- fname
    val <- paste0("function(x,p,...) p[\"",fpar[1],"\"]*(x==0)")
    for (i in seq(K-1)) {
        val <- paste0(val,"+p[\"",fpar[i+1],"\"]*(x==",i,")")
    }
    functional(x,formula) <- eval(parse(text=val))
    return(x)
}

##' @export
'categorical<-' <- function(x,...,value) categorical(x,value,...)
