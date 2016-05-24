dtn.mix=function(t,df,mu.ncp, sd.ncp, log=FALSE, approximation=c('int2','saddlepoint','laplace','none'),...)
{
    approximation=match.arg(approximation)
    if(all(is.infinite(df))) return( dnorm(t, mu.ncp, sqrt(1+sd.ncp*sd.ncp), log=log) )
    if(approximation!='none' && approximation!='int2' ) df[is.infinite(df)]=500 
	
	if(approximation=='none' && all(df==round(df))) approximation='int2'
    
	scale.fact=sqrt(1+sd.ncp*sd.ncp)
	ncp=mu.ncp/scale.fact
	
	fname=paste('dt', switch(approximation, int2='.int2',saddlepoint='.sad',laplace='.lap', none=''), sep='')
	
	lst=list(x=t/scale.fact, df=df, ncp=ncp, log=log, ...)
	
	ans=do.call(fname, lst)
	if(isTRUE(log)) ans-log(scale.fact) else ans/scale.fact
}
