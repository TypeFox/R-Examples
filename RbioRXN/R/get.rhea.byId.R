get.rhea.byId <-
function(rheaId) {
	if(length(rheaId) == 0) {
		cat('Please enter correct Rhea ID...\n')
	}
	
	result = data.frame()
	for(i in rheaId) {
		urlBase = 'http://www.ebi.ac.uk/rhea/rest/1.0/ws/reaction/biopax2/%s'
		url = sprintf(urlBase, i)

		# Gathering XML from Rhea REST

		h = basicTextGatherer()
		curlPerform(url = url, writefunction = h$update)

		biopax = h$value()
		biopax = unlist(strsplit(biopax, '\n'))

		parsedBiopax = .parse.rhea(biopax)	
		parsedBiopax = as.data.frame(parsedBiopax, stringsAsFactors = FALSE)
		result = rbind.fill(result, parsedBiopax)
	}
	result[is.na(result)] = ''
	return(result)
}
