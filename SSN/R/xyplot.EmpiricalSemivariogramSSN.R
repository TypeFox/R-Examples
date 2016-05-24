xyplot.EmpiricalSemivariogramSSN <- function(object, ...)
{
		plot(c(1,1),type = "n")
		t.background <- trellis.par.get("background")
		t.background$col <- "white"
		trellis.par.set("background", t.background)
		xyplot(gamma ~ distance | azimuth, data=object, pch = 19, col = "black",
			as.table = TRUE, ...)
}

