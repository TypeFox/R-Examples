`pmc` <-
function(x, ...) {
# proportion of trials misclassified by observer
#  w/ respect to estimated scale
	if (!(class(x) %in% c("mlds", "mlbs")))
		stop("x must be of class mlds or mlbs")
	pred <- predict(x)
	if (x$method == "glm") 
		resp <- x$obj$data$resp else
		resp <- x$data$resp
	pmc <- 1 - sum(((pred > 0) == resp))/length(resp)
	pmc		
	}

