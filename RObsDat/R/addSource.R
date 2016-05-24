addSource <- function(Organization, SourceDescription, SourceLink=NULL, ContactName=rep("Unknown", length(Organization)), Phone=rep("Unknown", length(Organization)), Email=rep("Unknown", length(Organization)), Address=rep("Unknown", length(Organization)), City=rep("Unknown", length(Organization)), State=rep("Unknown", length(Organization)), ZipCode=rep("Unknown", length(Organization)), Citation=rep("Unknown", length(Organization)), Metadata=rep("Unknown", length(Organization))){
	#optional: SourceLink
	stopifnot(length(Organization) == length(SourceDescription))
	stopifnot(length(Metadata) == length(SourceDescription))

	#checking for existing entries 
	for(i in seq(along=Organization)){
		if(NROW(existing <- getMetadata("Source",Organization=Organization[i], Description=SourceDescription[i], Citation=Citation[i], exact=TRUE))>0){
			warning(paste("Skipping existing ISOMetadata entry:", SourceDescription[i]))
			next
		}

		#check from reference table
		if(!is.null(Metadata[i])){
			MetadataID <- getID("ISOMetadata",Metadata[i])
		} else {
			MetadataID <- NA 
		}

		IaddSource(getOption("odm.handler"), Organization=Organization[i], SourceDescription= SourceDescription[i], SourceLink= SourceLink[i], ContactName= ContactName[i], Phone= Phone[i], Email= Email[i], Address= Address[i], City= City[i], State= State[i], ZipCode= ZipCode[i], Citation= Citation[i], Metadata=MetadataID)
	}

}
