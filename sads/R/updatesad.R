## updates a fitsad/fitrad object by running the optimizer again starting
## on better fit returned by profile
updatesad <- function(object, ...) {
	dots <- list(...)
	prof <- profile(object)
	if(class(prof) == "profile.mle2") stop("Cannot update, profile did not find a better fit!")
	newcall <- as.list(object@call)
	newcall[[1]] <- NULL # removes the "mle2" function name from the call
	newcall[["control"]] <- NULL # removes the "control" slot
	# merges the original call with arguments supplied by ...
	for (v in names(dots)) newcall[[v]] <- dots[[v]]
	# bugfix? profile returns coefficient "prob" name as "prob.prob", sometimes changes parameters
	# from start to fixed...
	names(newcall$start) -> name 
	for (v in name) {
		if (v %in% names(prof@coef))
			newcall$start[[v]] <- as.numeric(prof@coef[[v]])
		else # then it must be in fixed
			newcall$start[[v]] <- as.numeric(prof@call.orig$fixed[v])
	}
	newobj <- do.call("mle2", newcall)
	if(class(object) == "fitsad") 
		return (new("fitsad", newobj, sad=object@sad, distr=object@distr, trunc=object@trunc))
	else # fitrad
		return (new("fitrad", newobj, rad=object@rad, distr=object@distr, trunc=object@trunc, rad.tab=object@rad.tab))
}
updaterad <- function(object, ...) updatesad(object, ...)
