`getErrorReports` <- 
function(tmp, 
	tmp.tf, 
	par.sgp.config) {

	tmp.errors <- tmp[tmp.tf]
	tmp.config <- par.sgp.config[tmp.tf]
	tmp.error.report <- list()
	for (i in 1:length(tmp.errors)) tmp.error.report[[i]] <- list(Error=tmp.errors[[i]], Analysis=tmp.config[[i]])
	return(tmp.error.report)

} ### END getErrorReports
