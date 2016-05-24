get.rhea.all <-
function() {
	
	url = "ftp://ftp.ebi.ac.uk/pub/databases/rhea/biopax//rhea-biopax_full.owl.gz"
	tmpdest = tempfile(pattern = "rhea")
	download.file(url, destfile = tmpdest)
	owl = readLines(tmpdest)
  
	result = .parse.rhea(owl)
  	return(result)
}
