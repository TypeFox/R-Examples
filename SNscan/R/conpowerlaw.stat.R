conpowerlaw.stat <-
function (obs, pop=1, zloc,xmin=1,zetatable=NULL) 
{
	if(sum(length(zloc)==c(1,length(obs),length(obs)-1))>0) {t.stat=0;a0=0;az=0} else{
		if (xmin==1) {
			conalpha=function(x)
			{
				n=length(x)
				hatalpha=1+n/(sum(log(x)))
				return(hatalpha)
			}
			n0 = length(obs)
			nc=length(obs[-zloc])
			nz=length(zloc)
			m0 = sum(log(obs))
			mc = sum(log(obs[-zloc]))
			mz = sum(log(obs[zloc]))
			az=conalpha(x=obs[zloc])
			ac=conalpha(x=obs[-zloc])
			a0=conalpha(x=obs)
			t.stat = ((nz*log(az-1)-az*mz)+(nc*log(ac-1)-ac*mc)-(n0*log(a0-1)-a0*m0))*(az<=ac)+
				0*(az>ac)
		} else{
			Oz = conpl$new(obs[zloc])
			Oc = conpl$new(obs[-zloc])
			O0 = conpl$new(obs)
			Oz$setXmin(xmin)
			Oc$setXmin(xmin)
			O0$setXmin(xmin)
			az=estimate_pars(Oz)$pars
			ac=estimate_pars(Oc)$pars
			a0=estimate_pars(O0)$pars
			n0 = length(obs)
			nc=length(obs[-zloc])
			nz=length(zloc)
			m0 = sum(log(obs))
			mc = sum(log(obs[-zloc]))
			mz = sum(log(obs[zloc]))
			t.stat = ((nz*log(az-1)-az*mz)+(nc*log(ac-1)-ac*mc)-(n0*log(a0-1)-a0*m0))*(az<=ac)+
				0*(az>ac)
		}
	}
    return(c(t.stat, a0, az))
}
