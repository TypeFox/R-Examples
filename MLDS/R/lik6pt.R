`lik6pt` <-
function(x, Six.Pts, ...) {	
	p <- if (x$method == "glm")
		c(predict(x, make.ix.mat(Six.Pts$A, length(x$pscale)), 
			type = "response"),
		  predict(x, make.ix.mat(Six.Pts$B, length(x$pscale)), 
			type = "response"),
		  predict(x, make.ix.mat(Six.Pts$E, length(x$pscale)), 
			type = "response") ) else
		c(predict(x, Six.Pts$A, type = "response"),
		  predict(x, Six.Pts$B, type = "response"),
		  predict(x, Six.Pts$E, type = "response")) 
	grt <- with(Six.Pts, c(A$resp, B$resp, E$resp))
	grt %*% log(p) + (1 - grt) %*% log(1 - p)
}

