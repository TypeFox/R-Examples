combn2R=function (x, m, x2, m2, R, FUN = c, simplify = TRUE,
            sample.method="all", try.rest=TRUE, ...) 
{
    if(!is.function(FUN))stop("FUN needs to be a function")
    missingR=as.integer(missing(R))
    oneGrp=FALSE
    if(missing(x2) || missing(m2)) {m2=0;oneGrp=TRUE}
    if(length(m) != 1 || length(m2)!=1 || m<0 || m2 < 0) stop("m and m2 needs to be non-negative scalar")
    if (m == 0){ oneGrp=TRUE;x=x2;m=m2 }
    if (m2==0) oneGrp=TRUE

    if (is.numeric(x) && length(x) == 1 && x > 0 && trunc(x) == x)    x <- seq.int(x)
    n <- length(x)
    m <- as.integer(m)
    if(n<m) stop("n needs to be as large as m")
    count <- round(choose(n, m))
    sample.method=match.arg(sample.method,c("diff2","all","replace","noReplace"))
    if (oneGrp) {
         if(missingR){
            R=switch(sample.method,all=count,diff2=count,replace=NA)
            if(is.na(R))stop("R cannot be missing when sample.method='replace'")
         }
         if(length(R)!=1 || R<0) stop("R needs to be non-negative scalar")
         if(R==0) return(if (simplify) vector(mode(x), 0) else list())
         if (sample.method=='all' || sample.method=='diff2'|| sample.method=='noReplace'){
            if (R>count) R=count
            if (try.rest){
                trythis=try(combnBtRc(x=x,n=n, m=m,R=R,FUN=FUN,simplify=simplify,...=...))
                if (class(trythis)=='try-error') sample.method='replace'
                else return(trythis)
            }
         }
         if (sample.method=='replace')
            return(combnRRepl(x=x,n=n,m=m,R=R,FUN=FUN,simplify=simplify,...=...))  
    }

    if (is.numeric(x2) && length(x2) == 1 && x2 > 0 && trunc(x2) == x2)   x2 <- seq.int(x2) 
    n2 <- length(x2)
    m2 <- as.integer(m2)  
    if (n2 < m2) stop("n2 needs to be as large as m2")

    count2 <- round(choose(n2, m2))
    count.total=count*count2 
    mincount=min(c(count,count2))

    if(missingR){
        R=switch(sample.method,all=count.total,diff2=mincount,replace=NA)
        if(is.na(R))stop("R cannot be missing when sample.method='replace'")
    }
    if(length(R)!=1 || R<0) stop("R needs to be non-negative scalar")
    if(R==0) return(if (simplify) vector(mode(x), 0) else list())


    if (sample.method=="noReplace") stop("sample.method=noReplace is only for one group situation.")
    if (sample.method=='diff2'){
        if(R>mincount){
            if(try.rest) sample.method='all'
            else R=mincount
        }
        if(R<=mincount){
            if(try.rest){
                trythis=try(combnBt2RAllDifc(x=x,n=n,m=m,x2=x2,n2=n2,m2=m2,R=as.integer(R),
                            FUN=FUN,simplify=simplify,...=...),silent=TRUE)
                if(class(trythis)=='try-error')sample.method='replace'
                else return(trythis)
            }else{
                return(combnBt2RAllDifc(x=x,n=n,m=m,x2=x2,n2=n2,m2=m2,R=as.integer(R),
                        FUN=FUN,simplify=simplify,...=...))
            }
        }
    }
    if(sample.method=="all"){
        if (R>=count.total) R=count.total
        if(try.rest){
            trythis=try(combnBt2Rc(x=x,n=n,m=m,x2=x2,n2=n2,m2=m2,R=as.integer(R),
                        FUN=FUN,simplify=simplify,...=...),silent=TRUE)
            if(class(trythis)=='try-error') sample.method='replace'
            else return(trythis)
        }else return(combnBt2Rc(x=x,n=n,m=m,x2=x2,n2=n2,m2=m2,R=as.integer(R),
                    FUN=FUN,simplify=simplify,...=...))
    }
    if (sample.method=="replace")
        return(combn2RRepl(x=x,n=n,m=m,x2=x2,n2=n2,m2=m2,R=R,FUN=FUN,simplify=simplify,...=...))
}
