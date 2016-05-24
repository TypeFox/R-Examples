ec <-
function(mg1,mg2,sdg1,sdg2, df){
    l1=length(mg1);l2=length(mg2)
    y=sum(l2*mg1)-sum(l1*mg2)
    d1=rep(1,l1);d2=rep(1,l2)
    d1=l2*d1;d2=l1*d2;d=c(d1,d2);d=d^2
    sdg1=sdg1^2;sdg2=sdg2^2
    sd=c(sdg1,sdg2); sd=sd*d;sd=sum(sd); sd=sd^0.5
    tcal=y/sd; tcal=(tcal^2)^0.5
    p=(1-pt(tcal,df))*2; grupos="group.1 vs group.2"
    res=data.frame(grupos, contrast=round(y,2),standard.error=round(sd,4), tcal=round(tcal,2), p.value=round(p,4))
    return(res)
}
