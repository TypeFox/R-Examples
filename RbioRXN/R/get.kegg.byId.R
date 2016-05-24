get.kegg.byId <-
function(keggId) {
    kegg = data.frame()
    i = 1
    while(i <= length(keggId)) {

        cat('processing', keggId[i], '\n')
        query <- keggGet(keggId[i:(i+9)])

        for(l in 1:length(query)) {

            keggRow = query[[l]]
    
            for(j in names(keggRow)) {
                if(j == 'DBLINKS') {
                    for(k in 1:length(keggRow$DBLINKS)) {
                        db = unlist(strsplit(keggRow$DBLINKS[k], ': '))[1]
                        id = unlist(strsplit(keggRow$DBLINKS[k], ': '))[2]
                        keggRow[[db]] = id
                    }
                } else if (j == 'PATHWAY') {
                for(k in 1:length(keggRow$PATHWAY)) {
                    keggRow$PATHWAY[k] = paste(names(keggRow$PATHWAY[k]), keggRow$PATHWAY[k], sep=': ')
                }
                keggRow$PATHWAY = paste(keggRow$PATHWAY, collapse='///')
                } else if (j == 'REFERENCE') {
                    keggRow$REFERENCE = paste(keggRow$REFERENCE[[1]]$REFERENCE, collapse='///')
                } else {
                    if(length(keggRow[[j]]) > 1) {
                        keggRow[[j]] = paste(keggRow[[j]], collapse='///')
                    }
                }
            }
            keggRow[['DBLINKS']] = NULL
            keggRow = as.data.frame(keggRow, stringsAsFactors=FALSE)
            kegg = rbind.fill(kegg, keggRow)
            kegg[is.na(kegg)] = ''
        }
        i = i + 10 
    }
    return(kegg)
}
