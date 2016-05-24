`Stem.Bootstrap` <-
function(StemModel, B) {

	seed.list = list()
	for (i in 1:B) {
		seed.list[[i]] = as.integer(runif(1,min=-1,max=1)*(10^8))
	}

	output = lapply(seed.list, Stem.Bootstrap.fn, StemModel = StemModel, seed.list.out = seed.list)
	return(list(boot.output=output))
}

