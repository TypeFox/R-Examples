exact2x2Plot<-function(x, y=NULL,
    OR=NULL,
    ndiv=1000,
    tsmethod=NULL,
    method=NULL,
    paired=FALSE,
    orRange=NULL,
    dolog=TRUE,
    dolines=FALSE,
    dopoints=TRUE,
    doci=TRUE, 
    alternative = c("two.sided", "less", "greater"), 
    conf.level=.95,
    alphaline=TRUE,
    newplot=TRUE,...){

    ## method parameter is only for backwards compatability
    if (!is.null(method) & !is.null(tsmethod)){
        stop("only one of method or tsmethod can be non-NULL")
    } else if (!is.null(method)){
        tsmethod<-method
    }


    # copied setup code from fisher.test, july 1,2009, deleted some things not needed
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
        OK <- complete.cases(x, y)
        x <- factor(x[OK])
        y <- factor(y[OK])
        if ((nlevels(x) < 2) || (nlevels(y) < 2)) 
            stop("'x' and 'y' must have at least 2 levels")
        x <- table(x, y)
    }
    alternative<-match.arg(alternative)
    alternative <- char.expand(alternative, c("two.sided", 
        "less", "greater"))
    if (length(alternative) > 1 || is.na(alternative)) 
        stop("alternative must be \"two.sided\", \"less\" or \"greater\"")
    # end of setup code from fisher.test

    ## default to usual two-sided Fisher's exact test for paired=FALSE
    ## and usual exact McNemar confidence intervals for paired=TRUE
    if (is.null(tsmethod) & paired){ tsmethod<-"central"
    } else if (is.null(tsmethod)) tsmethod<-"minlike"


    if (is.null(orRange)){
        ## use ci to get a default range for the plot
        eout<-exact2x2(x,tsmethod=tsmethod,alternative=alternative,paired=paired,conf.level=conf.level)
        ci<-eout$conf.int
        est<-eout$estimate
        lowerRange<-ifelse(ci[1]==0, .8*min(est,ci[2]), .8*min(ci[1],est) )
        upperRange<-ifelse(ci[2]==Inf, 1.2*max(est,ci[1]), 1.2*max(ci[2],est) )
        
        orRange<-c(lowerRange,upperRange)
        if (ci[1]==0 & !dolog){
            orRange[1]<-0
        }
    }
    if (dolog & is.null(OR)){
        minor<-log10(orRange[1])
        maxor<-log10(orRange[2])
        OR<-10^(minor + (0:ndiv)*(maxor-minor)/ndiv )
    } else if (is.null(OR)){
        minor<-(orRange[1])
        maxor<-(orRange[2])
        OR<-(minor + (0:ndiv)*(maxor-minor)/ndiv )
    }

    alternative <- match.arg(alternative)
    if (paired){
        #bb<-x[1,2]
        #cc<-x[2,1]
        #x<-bb
        #n<-bb+cc
        p<-OR/(1+OR)
        ## earlier versions used method not tsmethod
        if (packageDescription("exactci",fields="Version")<="1.0-0"){
            stop("update exactci package to a Version>1.0-0")
            # for old version use method=tsmethod
            # see following lines...
            #p2sided<-exactbinomPvals(x[1,2],
            #    x[1,2]+x[2,1],p,
            #    method=tsmethod)$pvals
        } else {
            p2sided<-exactbinomPvals(x[1,2],x[1,2]+x[2,1],p,tsmethod=tsmethod)$pvals
        }
        pval<-switch(alternative,
          less=pbinom(x[1,2],x[1,2]+x[2,1],p),
          greater=pbinom(x[1,2]-1,x[1,2]+x[2,1],p,lower.tail=FALSE),
          two.sided=p2sided)
    } else {
        pval<-exact2x2Pvals(x,OR,tsmethod=tsmethod,alternative=alternative)$pvals
    }

    if (newplot){
        if (dolog){ LOG<-"x"
        } else LOG<-""
        plot(OR,pval,log=LOG,xlab="Null Hypothesis Odds Ratio",ylab="p-value",...)
    } 
    if (dopoints){
        points(OR,pval,...)
    } 
    if (dolines){
        lines(OR,pval,...)
    }
    if (doci){
        ci<-exact2x2(x,tsmethod=tsmethod,alternative=alternative,paired=paired,
            conf.level=conf.level)$conf.int
        alpha<-1-conf.level
        if (alphaline) lines(c(ci[1]*.1,ci[2]*10),c(alpha,alpha),lty=2)
        lines(c(ci[1],ci[1]),c(0,alpha),...)
        lines(c(ci[2],ci[2]),c(0,alpha),...)
    }

}
