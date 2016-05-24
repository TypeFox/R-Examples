`plot.gp.list` <-
function(x, sds = 1, CI.at.point=FALSE, ...) {
	createWindow(x$numGPs)
	for (i in 1:x$numGPs) {
		main = x$name[i]
		plot(x[[i]], type = 1 , sds = sds, CI.at.point=CI.at.point, main = main)
	}
}

