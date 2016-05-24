addSynonym <- function(table, phrase, id){
	#Replace ' by ` in phrases
	phraseA <- gsub("'", "`", phrase)

	if(length(answer <- IgetSynonymID(getOption("odm.handler"), table=table, phrase=phraseA))>0) return(answer)
	id  <- getID(table, id, remove.special.character=FALSE)
	IaddSynonym(getOption("odm.handler"), table=table, phrase=phraseA, id=id)
	return(id)
}
