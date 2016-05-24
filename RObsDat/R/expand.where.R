expand.where <- function(w.o, var, var.name, exact=TRUE, isnumeric=FALSE){
	

  if(exact){
		like = " = "
		perc = ""
	} else {
		like = " like "
		perc = "%"
	}
	if(isnumeric){
		quotation = ""
	} else {
		quotation = "'"
	}
  
	if(!is.null(var)){
		if(isnumeric & !is.numeric(var)){ 
			warning ("Numeric variables need to be numeric if isnumeric is set. Trying to convert")
			var.old <- var
			var <- as.numeric(var.old)
			print(paste("From", var.old, "to", var))
		} 

		if(!all(is.na(var))){
			if(length(var>1)){
				orterm <- ""
				theor <- ""
				for(i in seq(along=var)){
					orterm <- paste(orterm, theor, var.name, like, quotation, perc, var[i], perc, quotation, sep="")
					theor <- " OR "
				}
				w.o$where.clause <- paste(w.o$where.clause, w.o$the.and, " ( ", orterm, ' ) ', sep="")
			} else {
				w.o$where.clause <- paste(w.o$where.clause, w.o$the.and, var.name, like, quotation, perc, var, perc, quotation, sep="")
			}
			w.o$the.and <- " AND "
		}
	}
	return(w.o)
}
