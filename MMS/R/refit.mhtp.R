refit.mhtp=function(object,Ynew,z,grp,D,fix,rand,alpha,step,num,ordre,m,show,IT,maxq,speed,...)
{
if(missing(D))
{
	if("mhtp_nodiag"%in%class(object))
	{			out=refit.mhtp_nodiag(object,Ynew,z,grp,fix,rand,alpha,step,num,ordre,m,show,IT,maxq,speed,...)
	}
	if("mhtp_diag"%in%class(object))
	{					out=refit.mhtp_diag(object,Ynew,z,grp,fix,rand,alpha,step,num,ordre,m,show,IT,maxq,speed,...)
	}
}else{
	if(D==0){
out=refit.mhtp_nodiag(object,Ynew,z,grp,fix,rand,alpha,step,num,ordre,m,show,IT,maxq,speed,...)
	}else{
out=refit.mhtp_diag(object,Ynew,z,grp,fix,rand,alpha,step,num,ordre,m,show,IT,maxq,speed,...)
	}
	
	
}
	
out$call=match.call()
		
out	
}