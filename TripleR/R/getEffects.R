merge_recurse <- function (dfs, ...) {
    if (length(dfs) == 2) {
        merge(dfs[[1]], dfs[[2]], all = TRUE, sort = FALSE, ...)
    }
    else {
        merge(dfs[[1]], Recall(dfs[-1]), all = TRUE, sort = FALSE, 
            ...)
    }
}


#' @export
getEffects <- function(formule, data, varlist, by=NA, na.rm=TRUE, minVar=localOptions$minVar, gm=FALSE, ...) {

	# run a RR analysis for each variable and store results in a list
	res_list <- list()
	for (v in 1:length(varlist)) {
		print(paste("Calculate:",varlist[v]))
		f1 <- formula(paste(varlist[v], paste(as.character(formule), collapse="")))
		RR1 <- RR(f1, data=data, na.rm=na.rm, verbose=FALSE, minVar=minVar, ...)
		
		if (gm==FALSE) {
			eff <- RR1$effects
		} else {
			eff <- RR1$effects.gm
		}
		res_list <- c(res_list, list(eff))
	}

	# now combine all effects in a single data frame; merge by id
	if (is.na(by)) {
		if (length(RR1$groups) > 1) {by <- c("id", "group.id")} 
		else {by <- "id"}
	}
	
	res <- merge_recurse(res_list, by=by)
}