`extrp.pi0` <-
function(dat,slope.constraint=TRUE,
         gamma2.range=2^c(-4,3),
         rate.margin=c(0.5,0.5),
         plotit=TRUE
         ) {
#library("limSolve")

    gamma2.range=range(gamma2.range)
    stopifnot(gamma2.range[1]>0)
    rate.margin=range(rate.margin)
    stopifnot(rate.margin[1]<=0,rate.margin[1]>-.5,rate.margin[2]>=0)
    
    known.gamma2=(gamma2.range[1]==gamma2.range[2])

    if (slope.constraint && rate.margin[1]==rate.margin[2]) freeparm='none'
    else if (slope.constraint) freeparm='rate'
    else if (rate.margin[1]==rate.margin[2]) freeparm='slope'
    else freeparm='both'

    
    y=dat[,'f1']
#    idx=y!=1
#    if(sum(idx)<=1)  ### need fix
#    y=y[idx]
    idx=!is.na(y)
    y=y[idx]
    n1=dat[,'n1'][idx]
    n2=dat[,'n2'][idx]


    if(known.gamma2){
        if (freeparm=='none'){
            ## direct solution available
            anslist=extrp.pi0.only(n1,n2,y,gamma2.range[1])
        }else if (freeparm=='slope'){
            ## direct QP solution available
            anslist=extrp.pi0.slope(n1,n2,y,gamma2.range[1])
        }else if (freeparm=='rate') {
            ## 1-dim search on c
            anslist=extrp.pi0.rate(n1,n2,y,gamma2.range[1],rate.margin+.5)
        }else if (freeparm=='both') {
            ## 1-dim search on c, with QP nested in it
            anslist=extrp.pi0.both(n1,n2,y,gamma2.range[1],rate.margin+.5)
        }
    }else{## unknown.gamma2
        if (freeparm=='none'){
            ## 1-dim search on gamma2
            anslist=extrp.pi0.gam2(n1,n2,y,gamma2.range)
        }else if (freeparm=='slope'){
            ## 1-dim search on gamma2 with QP
            anslist=extrp.pi0.slope.gam2(n1,n2,y,gamma2.range)
        }else if (freeparm=='rate'){
            ## 2-dim search on gamma2 and c
            anslist=extrp.pi0.rate.gam2(n1,n2,y,gamma2.range,rate.margin+.5)
        }else if (freeparm=='both'){
            ## 2-dim search on gamma2 and c with QP
            anslist=extrp.pi0.both.gam2(n1,n2,y,gamma2.range,rate.margin+.5)
        }
    }


    ans=anslist$pi0
    names(ans)='pi0'
    attr(ans,'par')=anslist$par
    attr(ans,'sse')=anslist$sse
    attr(ans,'fitted.obj')=anslist$fitted.obj
    attr(ans,'freeparm')=freeparm
    attr(ans,'start.val')=start
    attr(ans,'subt.data')=dat
    class(ans)='extrpi0'
    if(plotit)plot(ans)
    ans
}

