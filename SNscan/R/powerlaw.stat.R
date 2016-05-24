powerlaw.stat <-
function (obs, pop=1, zloc,xmin=1,zetatable=NULL) 
{
	if(sum(length(zloc)==c(1,length(obs),length(obs)-1))>0) {t.stat=0;a0=0;az=0} else{
		if (xmin==1) {
			n0 = length(obs)
			nc=length(obs[-zloc])
			nz=length(zloc)
			m0 = sum(log(obs))
			mc = sum(log(obs[-zloc]))
			mz = sum(log(obs[zloc]))
			az=zetatable[which.min(abs(-(mz/nz)-zetatable[,2])),1]
			ac=zetatable[which.min(abs(-(mc/nc)-zetatable[,2])),1]
			a0=zetatable[which.min(abs(-(m0/n0)-zetatable[,2])),1]
			zeta0=as.numeric(zeta(a0))
			zetaz=as.numeric(zeta(az))
			zetac=as.numeric(zeta(ac))
			t.stat = ((a0*m0+n0*zeta0)-(az*mz+nz*zetaz)-(ac*mc+nc*zetac))*(az<=ac)+
				0*(az>ac)
		} else{
			Oz = displ$new(obs[zloc])
			Oc = displ$new(obs[-zloc])
			O0 = displ$new(obs)
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
			zeta0=as.numeric(zeta(a0))
			zetaz=as.numeric(zeta(az))
			zetac=as.numeric(zeta(ac))
			t.stat = ((a0*m0+n0*zeta0)-(az*mz+nz*zetaz)-(ac*mc+nc*zetac))*(az<=ac)+
				0*(az>ac)
		}
	}
    return(c(t.stat, a0, az))
}
