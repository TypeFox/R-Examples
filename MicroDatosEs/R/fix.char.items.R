###################################################################
# cjgb
# 20140508
# Fixes problems with character items in data.sets from package memisc
###################################################################

fix.char.items <- function(x){
	where <- which(sapply(x, class) == "character.item")

	for(i in where)
		attributes(x[[i]])$value.labels@values <- gsub("'", "", attributes(x[[i]])$value.labels@values)

	x
}

