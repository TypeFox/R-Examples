`exact2x2` <-
function(x, y = NULL, or = 1, alternative = "two.sided", 
    tsmethod=NULL,
    conf.int = TRUE, conf.level = 0.95, 
    tol=0.00001,conditional=TRUE, paired=FALSE, plot=FALSE) 
{
    if (!conditional) stop("unconditional exact tests not supported")
    if (tol<.Machine$double.eps^.5) warning("tol set very small, make sure pnhyper in exact2x2CI is this accurate")
    # copied setup code from fisher.test, july 1,2009
    DNAME <- deparse(substitute(x))
    METHOD <- "Fisher's Exact Test for Count Data"
    if (is.data.frame(x)) 
        x <- as.matrix(x)
    if (is.matrix(x)) {
        if (any(dim(x) < 2)) 
            stop("'x' must have at least 2 rows and columns")
        if (!is.numeric(x) || any(x < 0) || any(is.na(x))) 
            stop("all entries of 'x' must be nonnegative and finite")
        if (!is.integer(x)) {
            xo <- x
            x <- round(x)
            if (any(x > .Machine$integer.max)) 
                stop("'x' has entries too large to be integer")
            if (!identical(TRUE, (ax <- all.equal(xo, x)))) 
                warning("'x' has been rounded to integer: ", 
                  ax)
            storage.mode(x) <- "integer"
        }
    }
    else {
        if (is.null(y)) 
            stop("if 'x' is not a matrix, 'y' must be given")
        if (length(x) != length(y)) 
            stop("'x' and 'y' must have the same length")
        DNAME <- paste(DNAME, "and", deparse(substitute(y)))
        OK <- complete.cases(x, y)
        ## 2011-4-20: fix so that if x and y are input as factors, they remain factors as input
        ##   previously had
        #
        #    x <- factor(x[OK])
        #    y <- factor(y[OK])
        ##
        ##   so if x<-factor(c(0,1,0)) and y<-factor(c(0,0,0),levels=c(0,1)) then you would get an 
        ##   error
        x <- as.factor(x[OK])
        y <- as.factor(y[OK])

        if ((nlevels(x) < 2) || (nlevels(y) < 2)) 
            stop("'x' and 'y' must have at least 2 levels")
        x <- table(x, y)
    }
    alternative <- char.expand(alternative, c("two.sided", 
        "less", "greater"))
    if (length(alternative) > 1 || is.na(alternative)) 
        stop("alternative must be \"two.sided\", \"less\" or \"greater\"")
    # end of setup code from fisher.test

    ## default to usual two-sided Fisher's exact test for paired=FALSE
    ## and usual exact McNemar confidence intervals for paired=TRUE
    if (is.null(tsmethod) & paired){ tsmethod<-"central"
    } else if (is.null(tsmethod)) tsmethod<-"minlike"

    nr <- nrow(x)
    nc <- ncol(x)
    if (nr!=2 || nc!=2) stop("table must be 2 by 2")

    if (paired){
        OUT<-mcnemar.exact.calc(x[1,2],x[2,1],or=or,alternative=alternative,tsmethod=tsmethod,conf.level=conf.level)
    } else {

        if (alternative=="less" | alternative=="greater"){
            OUT<-fisher.test(x,or=or,alternative=alternative,conf.int=conf.int,conf.level=conf.level)
            OUT$method<-"One-sided Fisher's Exact Test"
        } else {
            #xmat<-x
            #m <- sum(x[, 1])
            #n <- sum(x[, 2])
            #k <- sum(x[1, ])
            #x <- x[1, 1]
            tsmethod<-char.expand(tsmethod,c("minlike","blaker","central"))
            if (tsmethod!="blaker" & tsmethod!="minlike" & tsmethod!="central") stop("tsmethod must be 'blaker', 'central' or 'minlike'.")
            if (tsmethod=="minlike"){
                OUT<-fisher.test(x,or=or,alternative="two.sided",conf.int=FALSE) 
                if (conf.int){
                    OUT$conf.int<-exact2x2CI(x,tsmethod="minlike",conf.level=conf.level,
                        tol=tol)
                } else {
                    OUT$conf.int<-NULL
                }
                OUT$method<-"Two-sided Fisher's Exact Test (usual method using minimum likelihood)"
            } else if (tsmethod=="blaker"){
                OUT<-fisher.test(x,or=or,alternative="two.sided",conf.int=FALSE) 
                if (conf.int){
                    OUT$conf.int<-exact2x2CI(x,tsmethod="blaker",conf.level=conf.level,
                        tol=tol)
                } else {
                    OUT$conf.int<-NULL
                }
                OUT$p.value<-exact2x2Pvals(x,or,tsmethod="blaker")$pvals
                OUT$method<-"Blaker's Exact Test"
            } else if (tsmethod=="central"){
                OUT<-fisher.test(x,or=or,alternative=alternative,conf.int=conf.int,
                    conf.level=conf.level)
                pless<-fisher.test(x,or=or,alternative="less",conf.int=FALSE)$p.value
                pgreater<-fisher.test(x,or=or,alternative="greater",conf.int=FALSE)$p.value
                OUT$p.value<-min(1,2*min(pless,pgreater))
                OUT$method<-"Central Fisher's Exact Test"
            } 
        }
    }
    if (plot){
        ci<-OUT$conf.int
        if (is.null(ci)){
             ci<-fisher.test(x,or=or,alternative=alternative,
                    conf.level=conf.level)$conf.int
        }
        lowerRange<-ifelse(ci[1]==0, .8*min(or,ci[2]), .8*min(ci[1],or) )
        upperRange<-ifelse(ci[2]==Inf, 1.2*max(or,ci[1]), 1.2*max(ci[2],or) )
        exact2x2Plot(x, tsmethod = tsmethod, alternative=alternative, orRange=c(lowerRange,upperRange),paired=paired, doci=TRUE, alphaline=TRUE,col="gray",pch=16,cex=.5)
        points(or,OUT$p.value,col="black")
    }


    OUT$data.name<-DNAME
    return(OUT)
}