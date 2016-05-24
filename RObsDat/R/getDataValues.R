
#S3-Objekte
getDataValues <- function(ID=NULL, from=NULL, to=NULL,  tz=c("global", "UTC", "GMT", "0", "local"), Site=NULL, Variable=NULL, Offset=NULL, OffsetType=NULL, CensorCode=NULL, Qualifier=NULL, Method=NULL, Source=NULL, Sample=NULL, DerivedFromID=NULL, QualityControlLevel=NULL, VersionID=NULL, VersionDate=NULL,show.deleted=FALSE, all.ID=FALSE ){

  # if all parameters are ID's they can use so. When not the ID's have to established
	if(all.ID){
    SiteID <-  Site
		VariableID <-  Variable
		OffsetTypeID <-  OffsetType
		QualifierID <-  Qualifier
		MethodID <-  Method
		SourceID <-  Source
		SampleID <-  Sample
		QualityControlLevelID <-  QualityControlLevel
	} 
  
  else {    
		SiteID <- getID("Site", Site)
		VariableID <- getID("Variable", Variable)
		OffsetTypeID <- getID("OffsetType", OffsetType)
		QualifierID <- getID("Qualifier", Qualifier)
		MethodID <- getID("Method", Method)
		SourceID <- getID("Source", Source)
		SampleID <- getID("Sample", Sample)
		QualityControlLevelID <- getID("QualityControlLevel", QualityControlLevel)
	}
	
	all.args <- list(SiteID=SiteID, VariableID=VariableID, Offset=Offset, OffsetTypeID=OffsetTypeID, CensorCode=CensorCode, QualifierID=QualifierID, MethodID=MethodID, SourceID=SourceID, SampleID=SampleID, DerivedFromID=DerivedFromID, QualityControlLevelID=QualityControlLevelID, VersionID=VersionID, VersionDate=VersionDate)
	#check if only single constraints are used -> necessary
	# because we need to expand it in addDataValues
	# in order to provide smart treatment of meta data
	count.unique <- sapply(all.args, function(x) length(unique(x)))
	if(all(count.unique <= 1)){
		
		the.unique <- sapply(all.args, function(x) unique(x))
		for(i in seq(along=the.unique)){
			if(!is.null(the.unique[[i]])){
				assign(names(the.unique)[i], value=the.unique[[i]])
			}
		}
	} else if(sum(count.unique > 1)==1) {
    
    the.long <- which(count.unique > 1)
		for(i in seq(along=all.args)[-the.long]){
			the.unique <- unique(all.args[[i]])
			if(!is.null(the.unique)){
				assign(names(all.args)[i], value=the.unique)
			}
		}
	} else {
        
		todo("Smart processing of multiple arguments is missing")
		browser()
	}



	# Datensaetze mit groesster VersionsID und mit ValidUntilID <= aktuell verlangte Version existiert
	old.entry <- NULL
	if(!is.null(VersionDate)){
		if(!is.null(VersionID)){
			stop("Not possible to use both VersionID and VersionDate")
		}
		todo("implement VersionDate")
		
	}

  # get a older Data-Version. (only old/replaced information)
	if(!is.null(VersionID)){	  
		old.entry <- IgetOldDataValues(options("odm.handler")[[1]], ID=as.numeric(ID), from=from, to=to, tz=tz, SiteID=SiteID, VariableID=VariableID, Offset=Offset, OffsetTypeID=OffsetTypeID, CensorCode=CensorCode, QualifierID=QualifierID, MethodID=MethodID, SourceID=SourceID, SampleID=SampleID, DerivedFromID=DerivedFromID, QualityControlLevelID=QualityControlLevelID, VersionID)
	}

	if(show.deleted){
		deleted.entry <- IgetDeletedDataValues(options("odm.handler")[[1]], ID=as.numeric(ID), from=from, to=to, tz=tz, SiteID=SiteID, VariableID=VariableID, Offset=Offset, OffsetTypeID=OffsetTypeID, CensorCode=CensorCode, QualifierID=QualifierID, MethodID=MethodID, SourceID=SourceID, SampleID=SampleID, DerivedFromID=DerivedFromID, QualityControlLevelID=QualityControlLevelID)
		#merge with old data

		if(!is.null(deleted.entry)){
			if(!is.null(old.entry)){
					replace.by.delete.data <- old.entry$ValueID %in% deleted.entry$ValueID
					old.entry$DataValue[replace.by.delete.data] = deleted.entry$DataValue	
				}
			else {
				old.entry <-  deleted.entry
			}
			}
		}		
	
	# latest dataset
	
	entry <- IgetDataValues(options("odm.handler")[[1]], ID=as.numeric(ID), from=from, to=to, tz=tz, SiteID=SiteID, VariableID=VariableID, Offset=Offset, OffsetTypeID=OffsetTypeID, CensorCode=CensorCode, QualifierID=QualifierID, MethodID=MethodID, SourceID=SourceID, SampleID=SampleID, DerivedFromID=DerivedFromID, QualityControlLevelID=QualityControlLevelID)


	if(!is.null(old.entry)){
		# all updated and deleted data have to replaced.
		new.data <- entry
		
		if(show.deleted == FALSE){
			deleted.entry <- IgetDeletedDataValues(options("odm.handler")[[1]], ID=as.numeric(ID), from=from, to=to, tz=tz, SiteID=SiteID, VariableID=VariableID, Offset=Offset, OffsetTypeID=OffsetTypeID, CensorCode=CensorCode, QualifierID=QualifierID, MethodID=MethodID, SourceID=SourceID, SampleID=SampleID, DerivedFromID=DerivedFromID, QualityControlLevelID=QualityControlLevelID)
			if(nrow(deleted.entry)>0){
			deleted.entry$DataValue <- 'd'
			replace.by.delete.data <- old.entry$ValueID %in% deleted.entry$ValueID
			old.entry$DataValue[replace.by.delete.data] = deleted.entry$DataValue	
			}		
		}	
		
		replace.by.old.data <- entry$ValueID %in% old.entry$ValueID
		# we need all columns except VersionID for merge
		columns <- colnames(old.entry) %in% (c("ValueID", "DataValue", "ValueAccuracy", "LocalDateTime", "UTCOffset", "DateTimeUTC", "SiteID", "VariableID", "OffsetValue", "OffsetTypeID", "CensorCode", "QualifierID", "MethodID", "SourceID", "SampleID", "DerivedFromID", "QualityControlLevelID"))
		old.entry = old.entry[,columns]
		new_data = rbind(new.data[!replace.by.old.data,], old.entry)

		#change <NA> to NA	#unfortunately 'd' from line 100 is also removed.
		new_data$DataValue = as.numeric(new_data$DataValue)
		
		# sort to ValueID
		entry <- new_data[order(new_data$ValueID),]
	}
		
	stfdf_entry = restructureDataResult(entry, tz=tz)
	return(stfdf_entry)

}

