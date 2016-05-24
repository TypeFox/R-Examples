addSpatialReferences <- function(ID, SRSID, Name, IsGeographic, Notes){
	stopifnot(length(ID)==length(Name))
	stopifnot(length(ID)==length(SRSID))
	stopifnot(length(ID)==length(Notes))
	stopifnot(length(ID)==length(IsGeographic))
	stopifnot(is.logical(IsGeographic))
	if(grepl("'", Notes)){
		Notes <- gsub("'", "`", Notes)
		warning("Replacing ' by ` in Notes to avoid conflict with Queries")
	}
	for(i in seq(along=ID)){
		#check existing
		if(NROW(IgetSpatialReference(getOption("odm.handler"), ID=ID[i], SRSName=Name[i], SRSID=SRSID[i], IsGeographic=IsGeographic, Notes=Notes[i],exact=TRUE))>0){
			warning(paste("Skiping existing entry:", Name[i]))
			next
		}
		warning(paste("Extending SpatialReferences table which should not be necessary. Please propose new term to CUASHI at http://his.cuahsi.org/mastercvreg/", sep=""))
		IaddSpatialReference(getOption("odm.handler"), ID=ID[i], SRSName=Name[i], SRSID=SRSID[i], IsGeographic=IsGeographic, Notes=Notes[i])
	}
}
