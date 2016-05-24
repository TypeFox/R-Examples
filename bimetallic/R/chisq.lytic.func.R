chisq.lytic.func <-
function(R, gam_ca, gam_co, ppv, npv, homRR, N_co, maf, prev, model){
	out = do.call(chisq.setup, as.list(environment()))
	misclass.lambda = with(out, 
		(numP_CA+numM_CA)*(numP_CO+numM_CO)*sum((combGtFreq[1,]-combGtFreq[2,])^2/((numP_CA+numM_CA)*combGtFreq[1,]+(numP_CO+numM_CO)*combGtFreq[2,]))
	)
	}

