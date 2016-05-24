addQualityControlLevel <- function(ID, Code, Definition, Explanation=""){
	if(NROW((existing <- IgetQualityControlLevel(getOption("odm.handler"), ID=ID, exact=TRUE)))>0){
		warning(paste("Existing entry with ID:", ID, "Skiping import!"))
		return(existing)
	}
	IaddQualityControlLevel(getOption("odm.handler"), ID=ID, Code=Code, Definition=Definition, Explanation=Explanation)

}
