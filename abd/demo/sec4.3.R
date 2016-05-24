xYplot(Cbind(estimate,lower,upper) ~ sample,
         data=CIsim(n=100,samples=20, rdist=rnorm, args=list(mean=2622.0, sd=2036.9), 
					estimand=2622),
		 par.settings=col.mosaic(), 
		 groups=cover)
