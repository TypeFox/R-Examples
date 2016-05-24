`plotObservedEffectsDefault` <-
function(x, z, ...) {
	if (is.null(dim(x))) {
		plot(x,z, xlab = "p1")
	}
	else {
		createWindow(dim(x)[2])
		for (i in 1:dim(x)[2]) {
			plot(x[,i], z, xlab = paste("p", i))
		}
	}

}

