par.uin <- function() {
	if (exists("is.R") && is.function(is.R) && is.R()) {
		par("pin") / {u <- par("usr"); c(diff(u[1:2]), diff(u[3:4]))}
	} else  par()$uin

}
