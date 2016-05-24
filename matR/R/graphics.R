
#-----------------------------------------------------------------------------------------
#  functions in support of graphical rendering of analyses.
#-----------------------------------------------------------------------------------------

parAuto <- function (par.name, n) {
	switch (par.name,
		pch = {
			if (missing (n)) {
				pchByName ["circle.open"]
			} else if (n <= length (pchSets)) {
				unname (pchByName [pchSets [[n]]])
			} else if (n <= length (LETTERS)) {
				LETTERS [1:n]
			} else
				stop ("automatic \'pch\' values are insufficient; specify explicit values") },
		col = {
			if (missing (n)) {
				"grey50"
			} else if (n <= length (colSets)) {
				colSets [[n]]
			} else if (n <= 101) {
				paste0 ("grey", floor (seq (0, 85, len=n)))
			} else if (n <= length (colors ())) {
				colors () [floor (seq (1, length (colors ()), len=n))]
			} else
				stop ("automatic \'col\' values are insufficient; specify explicit values") },
		cex = {
			if (missing (n)) {
				0.8
			} else {
				min <- 1.0 / n
				max <- 2 - 1.0 / n
				seq (min, max, len=n)
			} },
		stop ("automatic \'", par.name, "\' values are unavailable; specify explicit values"))
	}

pchByName <- c(
	square.open = 0,
	circle.open = 1,
	triangle.open = 2,
	diamond.open = 5,
	invert.triangle.open = 6,
	square.solid = 15,
	circle.solid = 16,
	triangle.solid = 17,
	diamond.solid = 18,
	square.fill = 22,
	circle.fill = 21,
	triangle.fill = 24,
	diamond.fill = 23,
	invert.triangle.fill = 25,						
	plus = 3,
	times = 4,
	star = 8)

pchSets <- list(
	c("circle.solid"),
	c("square.solid", "triangle.solid"),
	c("square.solid", "triangle.solid", "circle.solid"),
	c("square.open", "triangle.open", "circle.open", "diamond.open"),
	c("square.open", "triangle.open", "circle.open", "diamond.open", "invert.triangle.open"),
	c("square.solid", "triangle.solid", "circle.solid", "square.open", "triangle.open", "circle.open"),
	c("square.solid", "triangle.solid", "circle.solid", "square.open", "triangle.open", "circle.open", "diamond.open"),
	c("square.solid", "triangle.solid", "circle.solid", "square.open", "triangle.open", "circle.open", "plus", "times"),
	c("square.solid", "triangle.solid", "circle.solid", "square.open", "triangle.open", "circle.open", "plus", "times", "star"))

colSets <- list(
	c ("slateblue"),
	c ("red", "slateblue"),
	c ("red", "slateblue", "chocolate4"),
	c ("red", "slateblue", "chocolate4", "darkorange"),
	c ("red", "slateblue", "chocolate4", "darkorange", "olivedrab4"))
