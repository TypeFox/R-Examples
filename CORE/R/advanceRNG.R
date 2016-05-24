advanceRNG <-
function(randopt=c("SIMPLE","RESCALE"),nrand,nevents){
	randopt<-match.arg(randopt)
	switch(randopt,
		SIMPLE=replicate(nrand,{a<-runif(2*nevents);rm(a)}),
		RESCALE=replicate(nrand,{a<-sample(nevents);a<-runif(nevents);rm(a)})
	)
	return()
}
