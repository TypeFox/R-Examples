addISOMetadata <- function(TopicCategory="Unknown", Title="Unknown", Abstract="Unknown", ProfileVersion="Unknown", MetadataLink=NULL){

	theUnexisting <- rep(TRUE, length(Title))
	if(NROW(existing <- getMetadata("ISOMetadata",TopicCategory=TopicCategory, Title=Title, exact=TRUE))>0){
		warning(paste("Skiping existing ISOMetadata entry:", paste(existing[,3], collapse="; ")))
		theUnexisting <- !(Title %in% existing[,3])
	}
	if(!any(theUnexisting)){
		return()
	}
	#check from reference table
	if(!is.null(TopicCategory)){
		stopifnot(length(TopicCategory) == length(Title))
		TheTopicCategory <- getID("TopicCategory",TopicCategory[theUnexisting])
	} else {
		TheTopicCategory <- rep(NA, sum(theUnexisting))
	}
	if(Abstract=="Unknown") Abstract <- rep("Unknwon", length(Title))
	if(ProfileVersion=="Unknown") ProfileVersion <- rep("Unknwon", length(Title))
	if(Abstract=="Unknown") Abstract <- rep("Unknwon", length(Title))
	
	IaddISOMetadata(getOption("odm.handler"), TopicCategory=TheTopicCategory, Title= Title[theUnexisting], Abstract= Abstract[theUnexisting], ProfileVersion= ProfileVersion[theUnexisting], MetadataLink= MetadataLink[theUnexisting])

}
