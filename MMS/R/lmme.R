lmme=function(data,Y,z,grp,D,step,showit)
{
	q=ncol(z)
if(missing(D)&(q==1)){D=1}
if(missing(D)&(q>1)){D=0}

if(D==0){
	out=lmme_nodiag(data,Y,z,grp,step,showit)
	}else{
	out=lmme_diag(data,Y,z,grp,step,showit)
	}
out$call=match.call()
	
out
}


