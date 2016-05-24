"set.panel" <-

function(m = 1., n = 1.)



{



	par(mfrow = c(m, n))



	cat("plot window will lay out plots in a", m, "by", n, 



		"matrix ", fill = TRUE)



	invisible()



}

