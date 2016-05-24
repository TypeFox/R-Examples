fitTirmPartition <-
function(data, max.pop, n.boots.part=100){

	if (nrow(data) == 1){
		return(print("All capture counts are the same. Partitioning not possilbe."))
	}
	
	part.dat <- partitionCountData(data=data, n.boots.part=n.boots.part, max.pop=max.pop)
	
	tirm.part <- fitTirm(part.dat$low.2.thirds, max.pop)
	
	excluded <- part.dat$up.1.third
	
	model <- "Two.innate.rates.partitioned"
	
	likelihood <- tirm.part$likelihood
	
	ml.pop.size <- tirm.part$ml.pop.size + sum(excluded[,2])
	
	ml.na <- tirm.part$ml.na
	
	ml.nb <- tirm.part$ml.nb
	
	alpha <- tirm.part$alpha
	
	cap.ind <- tirm.part$cap.ind
	
	sampled.ind <- sum(data[,2]) # prior to partitioning
	
	sample.size <-  tirm.part$sample.size # after partitioning
	
	max.pop <- max.pop
	
	p.value.partition <- part.dat$p.2.3
	
				
	return(list(model=model, likelihood=likelihood, ml.pop.size=ml.pop.size, ml.na=ml.na, ml.nb=ml.nb, alpha=alpha, cap.ind=cap.ind, sampled.ind=sampled.ind, sample.size=sample.size, max.pop=max.pop, excluded=excluded, p.value.partition=p.value.partition))

}
