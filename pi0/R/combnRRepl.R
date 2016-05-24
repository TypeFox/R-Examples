combnRRepl=function (x, n, m, R, FUN, simplify, ...) 
{   ## this function should be called from combn2R ONLY,
    ## which guarantees x,n,m,x2,n2,m2, are all sensible
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



    for(nexttoselect in 1:R){
        a=sample(n,m)
        evalFUN()
    }

    if (simplify) dim(out)=dim.use
    attr(out,'sample.method')="replace"
    out
}

