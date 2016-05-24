#########################################################################
# GET CORES
#########################################################################

get_cores <- function(runs) {

	x <- max(detectCores(all.tests = FALSE, logical = TRUE) - 1, 1)

	for (i in x:1) {
		if (runs %% i == 0) {
			break
		}
	}
	x <- i
	if (is.na(x)) x <- 2

	if (Sys.info()[['sysname']] == 'Windows') x <- 1

	return(x)

}
