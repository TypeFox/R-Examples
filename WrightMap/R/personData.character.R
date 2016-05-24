personData.character <- function(thetas, p.type = NULL,...) {
	model <- CQmodel(p.est = thetas, p.type = p.type)
	return(personData(model))
}