get.rhea.byCompound <-
function(rheaCompound) {
  if(length(rheaCompound) == 0) {
    cat('Please enter correct query...\n')
  }

  urlBase = 'http://www.ebi.ac.uk/rhea/rest/1.0/ws/reaction/biopax2?q=%s'
  url = sprintf(urlBase, rheaCompound)

  # Gathering XML from Rhea REST

  h = basicTextGatherer()
  curlPerform(url = url, writefunction = h$update)

  biopax = h$value()
	biopax = unlist(strsplit(biopax, '\n'))

	rheaId = c()
	for(i in biopax) {
		if(grepl('<uri>', i)) {
			regexp = '(<uri>http://www.ebi.ac.uk/rhea/rest/1.0/ws/reaction/biopax2/)(.*)(</uri>)'
			id = sub(regexp, '\\2', i)
			rheaId = c(rheaId, id)
		}
	}
	rheaId = trim(rheaId)

	# Parsing

  parsedBiopax = get.rhea.byId(rheaId)
  return(parsedBiopax)
}
