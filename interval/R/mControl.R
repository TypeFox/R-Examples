`mControl`<- function(cm=NULL,nmc=10^3-1,seed=1234321,digits=12,
    p.conf.level=.99,
    setSEED=TRUE,
    tol.svd=10^-8,nwsr=10^3-1,np=10^3-1,
    tsmethod="central"){
    if (!is.numeric(nmc) || nmc <= 0) 
        stop("value of 'nmc' must be > 0")
    if (!is.numeric(digits) || digits <= 0) 
        stop("value of 'digits' must be > 0")
    if (!is.numeric(p.conf.level) || p.conf.level < 0 || p.conf.level>1) 
        stop("must have 0< p.conf.level < 1")
    if (!is.null(cm)){
        if (!is.matrix(cm) || !all(cm==1 | cm==0) )
        stop(" 'cm' must be either NULL or a matrix of zeros and ones")
    }
    if (!is.numeric(tol.svd) || tol.svd <= 0) 
        stop("value of 'tol.svd' must be > 0")
    if (!is.numeric(nwsr) || nwsr <= 0) 
        stop("value of 'nwsr' must be > 0")
    if (!is.numeric(np) || np <= 0) 
        stop("value of 'np' must be > 0")
    if (!is.logical(setSEED)){
        stop(" setSEED must be logical")
    }
    if (!(class(tsmethod)=="character" & (tsmethod=="central" | tsmethod=="abs"))){
        stop(" tsmethod must be 'central' or 'abs' ")
    }

    list(cm=cm,nmc=nmc,seed=seed,digits=digits,p.conf.level=p.conf.level,
        setSEED=setSEED,
        tol.svd=tol.svd,nwsr=nwsr,np=np,tsmethod=tsmethod)    
}
