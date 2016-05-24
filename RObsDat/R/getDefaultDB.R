getDefaultDB <- function(file="RODM.db"){
	odm.close()
	if(!is.na(file.info(file)$size)) {
		try(getID("Site", value="test"), silent=TRUE)
	} else {
		file=system.file(file, package="RObsDat")
		file.copy(file, ".")
		try(getID("Site", value="test"), silent=TRUE)
	}
}
