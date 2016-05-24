subex=function(dat,
              n1=round(ncol(dat)/2),
              n2=ncol(dat)-n1,
              f1method=c("lastbin","qvalue"),
              max.reps=20,
              balanced=FALSE,
              nparm=c(2,4),
              extrpFUN=c('constrOptim','genoud'),
              starts=c(pi0=.75,gam2=1,a=.5,c=.5), 
              plotit=TRUE
#              ,f1.args
#              ,extrp.args
){
    sub.pi0s=subt(dat,n1,n2,f1method,max.reps,balanced)     #,...=f1.args)
    ext.pi0=extrp.pi0(sub.pi0s,nparm,extrpFUN,starts,plotit)    #,...=extrp.args)
    pvals=matrix.t.test(dat,1,n1,n2)
    qvals=fdr(pvals,ext.pi0[1])
    ans=list(pi0=ext.pi0[1],
             extrp.fit=ext.pi0,
             pvalues=pvals,
             qvalues=qvals
            )
    class(ans)='subex'
    ans
}
