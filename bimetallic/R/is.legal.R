is.legal <-
function(R, gam_ca, gam_co, ppv, npv, homRR, N_co, maf, prev, model){
	#is the set of parameters legal?
	out = do.call(chisq.setup, as.list(environment()))
	with(out, all(popGtFreq>0))
	}

