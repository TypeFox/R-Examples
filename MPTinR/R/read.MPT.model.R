
.read.MPT.model <- function(model.filename) {
	whole <- readLines(model.filename)
	whole <- gsub("#.*", "", whole)
	model <- vector("list", length(whole))
	c2 <- 1
	c3 <- 1
	s.flag <- FALSE
	for (c1 in 1:length(whole)) {
		if (!(grepl("^[[:space:]]*$", whole[c1]))) {
			s.flag <- TRUE
			model[[c2]][c3] <- parse(text = whole[c1])[1]
			c3 <- c3 + 1
			fin <- c2
		}
		else {
			if (s.flag == TRUE) c2 <- c2 + 1
			c3 <- 1
			s.flag <- FALSE
		}
	}
	return (model[1:fin])
}
