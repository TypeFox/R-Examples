chisq.sim.factory <-
function(R, gam_ca, gam_co, ppv, npv, homRR, N_co, maf, prev, model){
	out = do.call(chisq.setup, as.list(environment()))
	chisq.sim.func = function(){
		p_CA = rmultnom(out$numP_CA, out$popGtFreq[1,])
		p_CO = rmultnom(out$numP_CO, out$popGtFreq[2,])
		m_CA = rmultnom(out$numM_CA, out$perGtFreq[1,])
		m_CO = rmultnom(out$numM_CO, out$perGtFreq[2,])
		xtable = rbind(cases=p_CA+m_CA, controls=p_CO+m_CO)
		OR = allelicOR(xtable)
		ch.test = chisq.test(xtable)
		stat = ch.test$statistic
		return(c(stat=stat, OR=OR))
	}
	return(chisq.sim.func)

}

