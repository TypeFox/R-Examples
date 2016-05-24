load.model <-
function(file) {
	name <- load(file=file);
	model <- get(name);
	return(model);
}