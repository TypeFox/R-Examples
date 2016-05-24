parse.item <- function(obj){
	varname <- names(obj);
	datevariables <- c("urn:ohmage:context:utc_timestamp", "urn:ohmage:context:timestamp");
	textvariables <- c("urn:ohmage:survey:privacy_state","urn:ohmage:user:id","urn:ohmage:context:location:latitude","urn:ohmage:context:location:longitude");
	factorvariables <- c("urn:ohmage:survey:id")
	
	if(varname %in% datevariables){
		return(parse.date(obj));		
	} else if(varname %in% textvariables) {
		return(sapply(obj[[1]]$values, parsevector))
	} else if(varname %in% factorvariables) {
		return(as.factor(sapply(obj[[1]]$values, parsevector)))
	} else if(substring(varname, 1, 17) == "urn:ohmage:prompt") {
		return(parse.prompt(obj));
	} else {
		stop("Dont know how to parse column: ", varname)
	}
}
