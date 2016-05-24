fdr_od <-
function(obsp,permp,pnm,ntests,thres,cl=.95,c1=NA){
      z_ = qnorm(1 - (1 - cl)/2) # two-tailed test
	pcount = rep(NA,length(permp))
	for(p_ in 1:length(permp)){
	   permp[[p_]][,pnm] = ifelse(permp[[p_]][,pnm]<thres,1,0)
	   pcount[p_] = sum(permp[[p_]][,pnm],na.rm=TRUE)
	}
	# over-dispersion parameter is observed variance of p in permuted data / expected
	p = mean(pcount,na.rm=TRUE)/ntests # estimate p
	e_vr = ntests*p*(1 - p)
	o_vr = var(pcount,na.rm=TRUE)
	if( is.na( c1 ) ) {
	   c1 = o_vr / e_vr
	   if( !is.na( c1 ) ) if( c1 < 1) c1 = 1
	} 

	nperm = length(permp)
	mo = ntests
	ro = sum(obsp < thres)
	vp = sum(pcount)
	vp1 = vp
  rslt = rep(NA,4)
  if(ro > 0){
		if(vp == 0) vp = 1
		mp = nperm * mo
		prod1 = vp / (mp-vp)
		prod2 = mo / ro - 1
		fdr = prod1 * prod2
	
		t1var = mp / (vp * (mp - vp))
		t2var = mo / (ro * (mo - ro))
	
		evar_jm = (t1var + t2var) * c1
	
		ul = exp(log(fdr) + z_ * sqrt(evar_jm))
		ll = exp(log(fdr) - z_ * sqrt(evar_jm))
	
		pi0 = (mo - ro) / (mo - (vp / nperm))
	
		rslt = c(fdr,ll,ul,pi0)
		rslt = ifelse(rslt > 1, 1, rslt) # FDR > 1 does not make sense, thus set to 1 in this case
		rslt = c(rslt,c1,ro,vp1)
	} 
	return(rslt)
}

