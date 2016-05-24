panel.psyfun <- function(x, y, n, lnk = "logit", ...) {
	xy.glm <- glm(cbind(n * y, n * (1 - y)) ~ x, 
			binomial(lnk))
	rr <- lattice::current.panel.limits()$xlim
	xx <- seq(rr[1], rr[2], len = 100)
	yy <- predict(xy.glm, data.frame(x = xx),  
			type = "response")
	lattice::panel.lines(xx, yy,  ...)
	}
