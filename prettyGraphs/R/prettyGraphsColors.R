prettyGraphsColors <-
function(){
	possibleColors = colors()
	possibleColors = possibleColors[-grep("white",possibleColors)]
	possibleColors = possibleColors[-grep("black",possibleColors)]
	possibleColors = possibleColors[-grep("snow",possibleColors)]
	possibleColors = possibleColors[-grep("ivory",possibleColors)]
	possibleColors = possibleColors[-grep("azure",possibleColors)]
	possibleColors = possibleColors[-grep("gray",possibleColors)]
	possibleColors = possibleColors[-grep("grey",possibleColors)]
	possibleColors = possibleColors[-grep("seashell",possibleColors)]
	possibleColors = possibleColors[-grep("steel",possibleColors)]
	possibleColors = possibleColors[-grep("light",possibleColors)]
	possibleColors = possibleColors[which(colSums(col2rgb(possibleColors)) < 650)]
	return(possibleColors)
}
