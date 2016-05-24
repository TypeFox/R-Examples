combnBtRc=function (x, n, m, R, FUN, simplify, ...) 
{   ## this function should be called from combn2R ONLY,
    ## which guarantees x,n,m are all sensible
    if(m>n-m && m!=n){
        complement=TRUE
        m=n-m
    }else
        complement=FALSE
    
    a <- 1:m

    r=FUN(x[a*(if(complement) -1 else 1)],...)
    len.r <- length(r)
    
    if (simplify) {
        out <- matrix(r, nrow = len.r, ncol = R)
        d <- dim(r)
        dim.use <- if (length(d) > 1) c(d, R)
                   else if (len.r > 1) c(len.r, R)
                   else c(d, R)
    }else {
        out <- vector("list", R)
    }

    if (simplify & complement)
    {	evalFUN=function() out[,nexttoselect]<<-FUN(x[-a],...)
    }else if (simplify)
    {	evalFUN=function() out[,nexttoselect]<<-FUN(x[a],...)
    }else  if (complement)
    {   evalFUN=function() out[[nexttoselect]]<<-FUN(x[-a],... )
    }else 
        evalFUN=function() out[[nexttoselect]]<<-FUN(x[a],... )

    
    
    count <- round(choose(n, m))

    if(R<count){
        oldwarn=options('warn')$warn
        options(warn=2) ## treat warnings as errors
        trythis=try({selected.idx=as.integer(sort(sample(count,R)))},silent=TRUE)
        options(warn=oldwarn)
        if (class(trythis)=='try-error') return(trythis)
        allCombns=as.integer(0)
    }else{
        allCombns=as.integer(1)
        selected.idx=as.integer(0)
    }

    env=environment()
    nexttoselect=as.integer(0)
    a=as.integer(rep(0,m))
    R=as.integer(R)
    .Call("combnBtRc", n, m, a, selected.idx, R, nexttoselect, quote(evalFUN()), env,allCombns, PACKAGE='pi0');


    if (simplify) dim(out)=dim.use
    attr(out,'sample.method')='noReplace'
    out
}

