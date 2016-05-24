exportSynonyms <- function(file){
	syn <- IgetSynonyms(getOption("odm.handler"))
	export.df <- data.frame()
	for(theTab in unique(syn$tab)){
		sel <- syn$tab == theTab
		sel.id <- syn[sel,2]
		db.entry <- getMetadata(table=theTab, EXACT=TRUE, ID=sel.id)
		if(theTab == "QualityControlLevel"){
			db.column <- which(tolower(names(db.entry)) %in% c("code"))

		} else {
			db.column <- which(tolower(names(db.entry)) %in% c("name", "term", "description", "title"))
		}
		stopifnot(length(db.column)==1)
		export.df <- rbind(export.df, data.frame(phrase=syn[sel,1], table=syn[sel,3], key=db.entry[,db.column]))
	}
	write.table(export.df, file=file, row.names=FALSE)
}
