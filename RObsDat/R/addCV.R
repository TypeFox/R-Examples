addCV <- function(table, term, definition){
	stopifnot(length(term)==length(definition))
	for(i in seq(along=term)){
		definition[i] <- gsub("'", "`", definition[i])
		definition[i] <- gsub("\"", "`", definition[i])
		if(NROW(existing <- IgetCV(getOption("odm.handler"), table=table, term=term[i], definition=definition[i], exact=TRUE))>0){
			warning(paste("Skipping existing entry:", term[i]))
			return()
		}
		warning(paste("Extending CV", table, " which should not be necessary. Please propose new term to CUASHI at http://his.cuahsi.org/mastercvreg/", sep=""))

		IaddCV(getOption("odm.handler"), table=table, term=term[i], definition=definition[i])
	}
}
