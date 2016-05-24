stats.BOLD <- 
function(taxon){
	space <- grep(" ", taxon)
	if(length(space) > 0) taxon[space] <- sapply(space, function(x) paste("\"", taxon[x], "\"", sep=""))
	taxon <- paste(taxon, collapse="%20")
	taxon <- gsub(" ", "%20", taxon)
	URL <- paste("http://v3.boldsystems.org/index.php/Public_SearchTerms?query=", taxon, sep="")
	res <- scan(file = URL, what = "", sep = "\n", quiet = TRUE)
	res <- res[grep("totalRecords", res)]
	as.numeric(gsub("[^[:digit:]]", "", res))
}