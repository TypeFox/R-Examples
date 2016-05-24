xyplot(Sepal.Length ~ Sepal.Width, groups = Species, data = iris, 
	main = "Some Iris Data",
	sub = "(R. A. Fisher analysized this data in 1936)",
	xlab = "sepal width (cm)",
	ylab = "sepal length (cm)",
	alpha = .5,        
	auto.key = list(columns = 3))   

