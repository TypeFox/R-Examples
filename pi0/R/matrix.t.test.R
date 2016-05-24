`matrix.t.test` <-
function(x,       MARGIN=1,
        n1=    if(MARGIN==1)floor(ncol(x)/2) else floor(nrow(x)/2),
        n2=    if(MARGIN==1)ncol(x)-n1 else nrow(x)-n1,
        pool=TRUE, pOnly=TRUE, tOnly=FALSE
        ) 
{
#        stopifnot(inherits(x,c('matrix','data.frame')))
        if(!inherits(x,c('matrix','data.frame'))) stop("x needs to be a matrix or data.frame")
        MARGIN=as.integer(MARGIN)[1]
        if(!(MARGIN==1 || MARGIN==2)) stop("MARGIN can only be 1 or 2")
        if(n1<1 || n2<1 || n1+n2<3) stop("sample sizes error")
        n1=as.integer(n1)[1]
        n2=as.integer(n2)[1]
        if({if(MARGIN==1)ncol(x) else nrow(x)}!=n1+n2) stop("sample sizes error")
        ntests=as.integer(if(MARGIN==1)nrow(x) else ncol(x))
        pool=as.integer(if(pool) 1 else 0)
        x=as.double(as.matrix(x))
       
        
#        x;MARGIN;n1;n2;ntests;pool
        tmp=.C("tstatistic",dat=x,n1=n1,n2=n2,ntests=ntests,MARGIN=MARGIN,pool=pool,tstat=rep(0.0,ntests),
                df=rep(0.0,ntests), PACKAGE='pi0')
                
        if(tOnly) return(tmp$tstat)

        pvals=pt(abs(tmp$tstat),if(pool==1)n1+n2-2 else tmp$df, lower.tail=FALSE)*2
        if(pOnly) return(pvals)
        
        return(list(stat=tmp$tstat,df=if(pool==1)n1+n2-2 else tmp$df, p.value=pvals))
}

