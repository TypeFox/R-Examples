#' Get the found variable subset(s)
#'
#' Get a list of variable indices/names of the found variable subsets.
#'
#' This method is used to get the names or indices of the variables used in specified variable subsets.
#'
#' @param object The GenAlg object returned by \code{\link{genAlg}}.
#' @param indices The indices of the subsets or empty if all subsets should be returned.
#' @param names Should the names or the column numbers of the variables be returned.
#' @return A logical matrix where each column represents a variable subset
#' @export
#' @example examples/subsets.R
subsets <- function(object, indices, names = TRUE) {
    if (!is(object, "GenAlg")) {
        stop("`object` must be of type `GenAlg`. See ?subsets for details.");
    }

    if (missing(indices)) {
        indices <- seq_len(ncol(object@subsets));
    }

	subs <- object@subsets[ , indices, drop = FALSE];
	if(names == TRUE && !is.null(colnames(object@covariates))) {
		cnames <- colnames(object@covariates);
		lapply(split(subs, rep(seq_len(length(indices)), each = nrow(subs))), function(s) { cnames[s]; });
	} else {
		lapply(split(subs, rep(seq_len(length(indices)), each = nrow(subs))), function(s) { which(s); });
	}
}
