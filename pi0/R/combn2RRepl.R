combn2RRepl=function (x, n, m, x2, n2, m2, R, FUN, simplify, ...) 
{   ## this function should be called from combn2R ONLY,
    ## which guarantees x,n,m,x2,n2,m2, are all sensible
    if(m>n-m && m!=n){
        complement=TRUE
        m=n-m
    }else
        complement=FALSE
    if(m2>n2-m2 && m2!=n2){
        complement2=TRUE
        m2=n2-m2
    }else
        complement2=FALSE    

    a <- 1:m
    a2 <- 1:m2
    r=FUN(x[a*(if(complement) -1 else 1)],x2[a2*(if(complement2) -1 else 1)],...)
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

    evalFUN=function(){
        r=FUN(x[a*(if(complement) -1 else 1)],x2[a2*(if(complement2) -1 else 1)],...)
        if(simplify){
            out[,L]<<-r
        }else{
            out[[L]]<<-r
        }
    }



    for(L in 1:R){
        a=sample(n,m)
        a2=sample(n2,m2)
        evalFUN()
    }


#    if (simplify) 
#        array(out, dim.use)
#    else out
    if (simplify) dim(out)=dim.use
    attr(out,'sample.method')="replace"
    out
}

