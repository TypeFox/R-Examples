load_posterior_predictive_samples <- function(posterior.predictive.sample.file){
	tmpenv <- environment()
	tmp <- load(posterior.predictive.sample.file,envir=tmpenv)
	posterior.predictive.samples <- lapply(tmp,get,envir=tmpenv)
	names(posterior.predictive.samples) <- tmp
	stopifnot(setequal(names(posterior.predictive.samples),
								c("observed.Fst","posterior.sample.Fst",
								"D","E","posterior.predictive.sample.size")))
	return(posterior.predictive.samples)
}