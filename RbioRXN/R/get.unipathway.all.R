get.unipathway.all <-
function(oboPath='') {
	if(oboPath == '') {
		 url = 'http://www.grenoble.prabi.fr/obiwarehouse/download/unipathway/public/unipathway.obo'
 		tmpdest = tempfile(pattern = "unipathway")
 		download.file(url, destfile=tmpdest) # download from url in 'tmpdest'
		obo = readLines(tmpdest)
	} else {
		obo = readLines(oboPath)
	}

	# Parsing
	cat('Start parsing.', length(obo), 'lines being processed\n')

	unipathwayRow = list()
	unipathway = data.frame()

	for(i in obo[7:length(obo)]) {
  	if(grepl('^\\[', i)) {
    	type = gsub('\\[', '', i)
    	type = gsub('\\]', '', type)
    	unipathwayRow[['type']] = type
  	} else if(grepl('^$', i)) {
    	if(unipathwayRow[['type']] == 'Term') {
      	unipathwayRow = as.data.frame(unipathwayRow, stringsAsFactors=FALSE)
      	unipathway = rbind.fill(unipathway, unipathwayRow)
    	}
    	unipathwayRow = list()
  	} else {
    	entry = unlist(strsplit(i, ': '))
    	unipathwayRow[[entry[1]]] = paste(unipathwayRow[[entry[1]]], entry[[2]], sep='///')
  	}
	}

	unipathway[is.na(unipathway)] = ''
	unipathway$def = sub('^"', '', unipathway$def)
	unipathway$def = sub('" \\[.*\\]$', '', unipathway$def)
	unipathway$synonym = sub('^"', '', unipathway$synonym)
	unipathway$synonym = sub('" \\[', ' [', unipathway$synonym)
	for(i in names(unipathway)) {
  	unipathway[[i]] = sub('^///', '', unipathway[[i]])
	}

	# Post-processing
	cat('Removing redundancy.\n')
	
	# In the field 'xref', I parsed database (KEGG, Rhea, MetaCyc) and ID and create new column with each database.

	unipathway_reaction = unipathway[unipathway$namespace %in% c('reaction', 'enzymatic_reaction'),]

	regexp = '(.+:[[:alnum:][:punct:]]+)( .+)'

	unipathway2 = data.frame()
	unipathway2Row = list()
	for(i in 1:nrow(unipathway_reaction)) {
		unipathway2Row = unipathway_reaction[i,]
  	xref = unipathway_reaction[i, 'xref']
		if(length(xref) > 0) {
	  	xref2 = unlist(strsplit(xref, '///'))
  		xref3 = sub(regexp, '\\1', xref2)
  		xref5 = list()
  		for(j in xref3) {
    		xref4 = unlist(strsplit(j, ':'))
    		unipathway2Row[[xref4[1]]] = paste(unipathway2Row[[xref4[1]]], xref4[2], sep='///')
  		}
		}
  	unipathway2Row = as.data.frame(unipathway2Row, stringsAsFactors = FALSE)
  	unipathway2 = rbind.fill(unipathway2, unipathway2Row)
  	unipathway2Row = list()
	}

	unipathway2$xref = NULL

	unipathway3 = unique(unipathway2)

	unipathway3[is.na(unipathway3)] = ""

	for(i in names(unipathway3)) {
  	unipathway3[[i]] = sub('^///','',unipathway3[[i]])
  	unipathway3[[i]] = sub('///$','',unipathway3[[i]])
	}

	# Remove master
	unipathway3 = unipathway3[-c(1,2),]

	# GO uniprot
	unipathway3$GO = sub(' ".*', '', unipathway3$GO)

	# Create new fields
	unipathway3[['direction']] = ''
	unipathway3[['part_of']] = ''
	unipathway3[['compoundId']] = ''
	unipathway3[['compoundName']] = ''
	unipathway3[['equation']] = ''
	unipathway3[['enzName']] = ''

	# Parse part_of
	unipathway_relationship = build.subtable(unipathway3, 'id', 'relationship', '///')
	unipathway_part_of = unipathway_relationship[grep('part_of', unipathway_relationship$relationship),]

	regexp = '(^part_of )(.*)( \\{cardinality=)(.*)(direction=")(.+)("\\})'

	unipathway_part_of[['part_of']] = sub(regexp, '\\2', unipathway_part_of$relationship)
	unipathway_part_of[['direction']] = sub(regexp, '\\6', unipathway_part_of$relationship)

	unipathway_part_of$direction[grep('^part_of', unipathway_part_of$direction)] = ''
	unipathway_part_of$part_of[grep('^part_of', unipathway_part_of$part_of)] = ''

	for(i in 1:nrow(unipathway_part_of)) {
		unipathway3[unipathway3$id == unipathway_part_of$id[i], 'part_of'] = paste(unipathway3[unipathway3$id == unipathway_part_of$id[i], 'part_of'], unipathway_part_of[i,'part_of'], sep='///')
		unipathway3[unipathway3$id == unipathway_part_of$id[i], 'direction'] = paste(unipathway3[unipathway$id == unipathway_part_of$id[i], 'direction'], unique(unipathway_part_of[i,'direction']), sep='///')
	}

	## parse has_input, has_output
	has_compound = unipathway_relationship[grep('^has_.+put_compound', unipathway_relationship$relationship),]
	has_compound = data.frame(has_compound)

	regexp = '^(has_.+put_compound )(.+)( \\{cardinality="[0-9]"\\} ! .+ ! )(.+)'

	has_compound[['compoundId']] = sub(regexp, '\\2', has_compound$relationship)
	has_compound[['compoundName']] = sub(regexp, '\\4', has_compound$relationship)

	for(i in 1:nrow(has_compound)) {
  	unipathway3[unipathway3$id == has_compound$id[i], 'compoundId'] = paste(unipathway3[unipathway3$id == has_compound$id[i], 'compoundId'], has_compound[i,'compoundId'], sep='///')
 	 unipathway3[unipathway3$id == has_compound$id[i], 'compoundName'] = paste(unipathway3[unipathway3$id == has_compound$id[i], 'compoundName'], has_compound[i,'compoundName'], sep='///')
	}

	index_part_of = grep('UPa:UER', unipathway3$part_of)

	for(i in index_part_of) {
  	part_of = unipathway3[i,'part_of']
  	part_of = unlist(strsplit(part_of, '///'))
  	unipathway3[i,'EC'] = paste(unique(unipathway3[unipathway3$id %in% part_of, 'EC']), collapse='///')
  	unipathway3[i,'GO'] = paste(unique(unipathway3[unipathway3$id %in% part_of, 'GO']), collapse='///')
  	unipathway3[i,'UNIPROT'] = paste(unique(unipathway3[unipathway3$id %in% part_of, 'UNIPROT']), collapse='///')
  	unipathway3[i,'enzName'] = paste(unique(unipathway3[unipathway3$id %in% part_of, 'name']), collapse='///')
  	unipathway3[i,'equation'] = paste(unique(unipathway3[unipathway3$id %in% part_of, 'def']), collapse='///')  
	}

	for(i in names(unipathway3)) {
 		unipathway3[[i]] = sub('^///','',unipathway3[[i]])
  	unipathway3[[i]] = sub('///$','',unipathway3[[i]])
	}

	ind_remove = grep('UPa:UER', unipathway3$id)
	unipathway4 = unipathway3[-ind_remove,]

	unipathway4$equation = sub('&gt;', '>', unipathway4$equation)
	unipathway4$equation = sub('^\\"', '', unipathway4$equation)
	unipathway4$equation = sub('\\.$', '', unipathway4$equation)

	unipathway5 = unipathway4[,c('id', 'enzName', 'equation', 'KEGG', 'RHEA', 'METACYC', 'EC', 'GO', 'UNIPROT', 'compoundId', 'compoundName')]

	# UniPathway compound
	cat('processing compound\n')
	unipathway_compound = unipathway[unipathway$namespace == 'compound',]
	unipathway_compound = unipathway_compound[-1,]

	unipathway_compound2 = data.frame()
	unipathway_compound2Row = list()
	for(i in 1:nrow(unipathway_compound)) {
  	unipathway_compound2Row = unipathway_compound[i,]
  	xref = unipathway_compound[i, 'xref']
  	if(length(xref) > 0) {
    	xref2 = unlist(strsplit(xref, '///'))
    	xref3 = sub(regexp, '\\1', xref2)
    	xref5 = list()
    	for(j in xref3) {
      	xref4 = unlist(strsplit(j, ':'))
      	unipathway_compound2Row[[xref4[1]]] = xref4[2]
    	}
  	}
  	unipathway_compound2Row = as.data.frame(unipathway_compound2Row, stringsAsFactors = FALSE)
  	unipathway_compound2 = rbind.fill(unipathway_compound2, unipathway_compound2Row)
  	unipathway_compound2Row = list()
	}

	unipathway_compound_sub = build.subtable(unipathway_compound2, 'id', 'synonym')

	for(i in 1:nrow(unipathway_compound_sub)) {
        if(grepl('RELATED FORMULA', unipathway_compound_sub[i,'synonym'])) {
            unipathway_compound_sub[i, 'type'] = 'formula'
        } else if(grepl('RELATED InChI', unipathway_compound_sub[i,'synonym'])) {
            unipathway_compound_sub[i, 'type'] = 'inchi'
        } else if(grepl('\\[KEGG\\]$', unipathway_compound_sub[i,'synonym'])) {
            unipathway_compound_sub[i, 'type'] = 'synonym.kegg'
        }
    }

	unipathway_compound_sub$synonym = gsub(' RELATED.*', '', unipathway_compound_sub$synonym)
	unipathway_compound_sub$synonym = gsub('"', '', unipathway_compound_sub$synonym)
	unipathway_compound_sub$synonym = gsub(' \\[KEGG\\]$', '', unipathway_compound_sub$synonym)

	for(i in 1:nrow(unipathway_compound_sub)) {
  	    unipathway_compound2[unipathway_compound2$id == unipathway_compound_sub[i,'id'], unipathway_compound_sub[i,'type']] = unipathway_compound_sub[i,'synonym']
	}

	unipathway_compound3 = unipathway_compound2[,c('id', 'type', 'name', 'KEGG', 'CHEBI', 'formula', 'inchi', 'synonym.kegg')]

	unipathway_compound3[is.na(unipathway_compound3)] = ''

	result = list()
	result[['reaction']] = unipathway5
	result[['compound']] = unipathway_compound3
	return(result)
}
