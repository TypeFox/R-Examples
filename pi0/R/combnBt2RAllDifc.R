combnBt2RAllDifc=function (x, n, m, x2, n2, m2, R, FUN, simplify, ...) 
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
    count2 <- round(choose(n2, m2))  ########

    oldwarn=options('warn')$warn
    options(warn=2)
    trythis=try({
        selected.idx=as.integer(sort(sample(count,R)))
        selected.idx2.all=as.integer(sample(count2,R))
        },silent=TRUE)
    options(warn=oldwarn)
    if(class(trythis)=='try-error') return(trythis)
    R1=length(selected.idx)
    R2perGrp1=as.integer(rep(1,R1))
    allCombns=as.integer(0)

    a=as.integer(rep(0,m))
    a2=as.integer(rep(0,m2))
    env=environment()
    L=as.integer(0)

    .Call("combnBt2Rc",n,  m,  n2,  m2, a, a2, 
            selected.idx, R1,selected.idx2.all, R2perGrp1,
             L, quote(evalFUN()), env, allCombns, PACKAGE='pi0')

#    out=if (simplify) 
#        array(out, dim.use)
#    else out
    if (simplify) dim(out)=dim.use
    attr(out,'sample.method')='diff2'
    out
}

