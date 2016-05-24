subt=function(dat,
              n1=round(ncol(dat)/2),
              n2=ncol(dat)-n1,
              f1method=c("lastbin","qvalue"),
              max.reps=if (balanced) 20 else 5,
              balanced=FALSE,
              ...
             )
{
    if(!inherits(dat,c('matrix','data.frame'))) stop("dat needs to be a matrix or data.frame")
    N=ncol(dat)
    G=nrow(dat)
    if(N!=n1+n2) stop("ncol(dat) does match n1+n2")
    f1method=f1method[1]
    if(f1method=='qvalue') {
        #library("qvalue")
        q.pi0=function(p,...)qvalue(p,...)$pi0
    }


    if(f1method=="lastbin"){
       histf1.stat=function(tstat,bw=.2,trunc=TRUE)
       {
        tstat=tstat[!is.na(tstat)]
        cutoff=qt(.5+bw/2, sub.n1+sub.n2-2)
        if(trunc)max(min(mean(abs(tstat)<=cutoff)/bw,1),0) else mean(abs(tstat)<=cutoff)/bw
       }
       get.sub.pi0=function(idx,idx2,sub.n1,sub.n2){
            tstat=matrix.t.test(dat[,c(idx,idx2)],1,sub.n1,sub.n2,tOnly=TRUE)
            histf1.stat(tstat)
       }
    }else if (f1method=='qvalue'){
        get.sub.pi0=function(idx,idx2,sub.n1,sub.n2){
            pvals=matrix.t.test(dat[,c(idx,idx2)],1,sub.n1,sub.n2)
            pvals=pvals[!is.na(pvals)]
            do.call("q.pi0",list(p=pvals,...=...))
        }
    }else{
        get.sub.pi0=function(idx,idx2,sub.n1,sub.n2){
            pvals=matrix.t.test(dat[,c(idx,idx2)],1,sub.n1,sub.n2)
            pvals=pvals[!is.na(pvals)]
            do.call(f1method,list(p=pvals,...=...))
        }
    }

    total.rslt=sum(sapply(seq.int(3,N),function(sub.N)
        {sub.n1=seq_len(n1);
         sub.n2=sub.N-sub.n1;
         idx=(sub.n2<=n2 & sub.n2>=1 & (!balanced | sub.n1==sub.n2) );
         if(any(idx)) sum(
            pmin(max.reps,choose(n1,sub.n1[idx])*choose(n2,sub.n2[idx]))
         )else 0
         }))
    sub.pi0=matrix(,total.rslt,3)
    cur.row=0
    for(sub.N in seq.int(3,N)){
        for(sub.n1 in seq_len(n1)){
            sub.n2=sub.N-sub.n1
            if((balanced && sub.n1!=sub.n2)  || sub.n2<=0 || sub.n2>n2)next
            R=min(max.reps,choose(n1,sub.n1)*choose(n2,sub.n2))
            combn2R.rslt=combn2R(n1,sub.n1,1:n2+n1,sub.n2,R,get.sub.pi0,
                    sample.method='diff2',try.rest=TRUE,sub.n1=sub.n1,sub.n2=sub.n2)
            n.rslt=seq_along(combn2R.rslt)
            sub.pi0[cur.row+n.rslt,]=cbind(combn2R.rslt,sub.n1,sub.n2)
            cur.row=cur.row+length(n.rslt)
            if(cur.row>=total.rslt) break
        }
        if(cur.row>=total.rslt) break
    }
    colnames(sub.pi0)=c('f1','n1','n2')
    attr(sub.pi0,'n1')=n1
    attr(sub.pi0,'n2')=n2
    attr(sub.pi0,'f1method')=f1method
    attr(sub.pi0,'max.reps')=max.reps
    attr(sub.pi0,'balanced')=balanced
    class(sub.pi0)=c('subt','matrix')
    sub.pi0
}
