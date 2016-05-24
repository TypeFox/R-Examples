get.kegg.all <-
function() {
	cmp <- keggList("compound")
	reactionEntry = keggList("reaction")

	cmpId = names(cmp)
	cmpId = sub('cpd:', '', cmpId)

	reactionEntry = names(reactionEntry)
	reactionEntry = sub('rn:', '', reactionEntry)
	
	keggReaction = get.kegg.byId(reactionEntry)
	keggReaction[is.na(keggReaction)] = ""

	keggCompound = get.kegg.byId(cmpId)
	keggCompound[is.na(keggCompound)] = ""
	
	# reference
	referIndex = grep('.+', keggReaction$REFERENCE)
	referId = keggReaction[grep('.+', keggReaction$REFERENCE), 'ENTRY']
	referIdUnique = unique(keggReaction[grep('.+', keggReaction$REFERENCE), 'ENTRY'])

	redundantIndex = c()
	for(i in referIdUnique) {
  	    index = grep(i, referId)
  	    index = referIndex[index[-1]]
  	    redundantIndex = c(redundantIndex, index)
	}

    if(length(redundantIndex) > 0) {
	    keggReaction_unique = keggReaction[-redundantIndex,]
    } else {
        keggReaction_unique = keggReaction
    }
	
	result = list()
	result[['reaction']] = keggReaction_unique
	result[['compound']] = keggCompound
    cat('# of reactions:', nrow(keggReaction_unique), '\n')
    cat('# of compounds:', nrow(keggCompound), '\n')
	return(result)
}
