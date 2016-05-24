mhtp=function(data,Y,z,grp,D,fix,rand,alpha,step,num,ordre,m,show,IT,maxq,speed)
{
q=ncol(z)
if(missing(D)&(q==1)){D=1}
if(missing(D)&(q>1)){D=0}

if(D==0){
out=mhtp_nodiag(data,Y,z,grp,fix,rand,alpha,step,num,ordre,m,show,IT,maxq,speed)
	}else{
out=mhtp_diag(data,Y,z,grp,fix,rand,alpha,step,num,ordre,m,show,IT,maxq,speed)
	}
out$call=match.call()

out	
}