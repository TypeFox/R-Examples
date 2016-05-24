combnBt2Rc=function (x, n, m, x2, n2, m2, R, FUN, simplify, ...) 
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

    
    count <- round(choose(n, m))
    count2 <- round(choose(n2, m2))
    count.total=count*count2

    if(R<count.total){
        oldwarn=options('warn')$warn
        options(warn=2) ## treat warnings as errors
        trythis=try(selected.idxtotal<-sort(round(sample(count.total,R))),
                    silent=TRUE)
        options(warn=oldwarn)
        if(class(trythis)=='try-error') return(trythis)
        selected.idx.all=ceiling(selected.idxtotal/count2)
        selected.idx=as.integer(unique(selected.idx.all))
        R1=length(selected.idx)
        selected.idx2.all=as.integer(selected.idxtotal-(selected.idx.all-1)*count2)
        R2perGrp1=as.integer(table(selected.idx.all))
        allCombns=as.integer(0)
    }else{
        allCombns=as.integer(1)
        selected.idx=as.integer(1)
        selected.idx2.all=as.integer(1)
        R1=as.integer(count)
        R2perGrp1=as.integer(rep(count2,count))
    }


    a=as.integer(rep(0,m))
    a2=as.integer(rep(0,m2))
    env=environment()
    L=as.integer(0)

    .Call("combnBt2Rc",n,  m,  n2,  m2, a, a2, 
            selected.idx, R1,selected.idx2.all, R2perGrp1,
             L, quote(evalFUN()), env, allCombns,
             PACKAGE='pi0')

#    if (simplify) 
#        array(out, dim.use)
#    else out
    if (simplify) dim(out)=dim.use
    attr(out,'sample.method')="all"
    out
}

