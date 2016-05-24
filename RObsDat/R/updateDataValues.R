updateDataValues <- function(getDataResult, reason=NULL){
	stopifnot(class(getDataResult)=="inherited_stfdf")

	inDB <- getDataValues(ID=unique(as.vector(as.matrix(getDataResult@ValueIDs@data))))

	if(any(dim(inDB@data)!=dim(getDataResult@data))){
		stop(paste("error while comparing data:\n dimension of data in DB: ", paste(dim(inDB@data), collapse=", "), "dimension of data passed to function: ", paste(dim(getDataResult@data),collapse=", ")))
	}
	# get only updated data sets
	comp <- getDataResult == inDB			# == doesn't work on STFDF. Therefore class_stfdf contains a method named ==

	#update the single records
	if(all(comp)){
		warning("Nothing to update")
		return()
	}
  
	to.update <- getDataResult@ValueIDs@data[comp]
	IarchiveDataValues(getOption("odm.handler"),ValueID=to.update, reason)
		
		
	for(rec.id in to.update){
		the.row <- which(getDataResult@ValueIDs@data==rec.id)
		
		metaID <- getDataResult@MetadataRel@data[the.row,]
		
		#make sure that ID is processed correctly in Iupdate if it is NA
		todo("Correct handling of tz in updateDataValues")
		the.tz <- "GMT"
		obj <- getOption("odm.handler")
	    selectMeta = which(getDataResult@Metadata$MetaId == metaID)
	   
		IupdateDataValues(obj,ValueID=rec.id, 
				localDateTime=index(getDataResult@time[the.row,]),
				value=sv(coredata(getDataResult@data[the.row,])),				
				valueAccuracy=sv(getDataResult@Metadata$Valueaccuracy[selectMeta]),			
				TZ = the.tz,
				SiteID=getID("Site",getDataResult@Metadata$site[selectMeta]),				
				VariableID=getID("Variable",getDataResult@Metadata$variable[selectMeta]),				
				Offset=sv(getDataResult@Metadata$offsetvalue[selectMeta]),					
				OffsetTypeID=getID("OffsetType",getDataResult@Metadata$offsettype[selectMeta]),				
				CensorCode=sv(getDataResult@Metadata$censorcode[selectMeta]),				
				QualifierID=getID("Qualifier",getDataResult@Metadata$qualifier[selectMeta]),				
				MethodID=getID("Method",getDataResult@Metadata$method[selectMeta]),			
				SourceID=getID("Source",getDataResult@Metadata$source[selectMeta]),				
				SampleID=getID("Sample",getDataResult@Metadata$sample[selectMeta]),				
				DerivedFromID= sv(getDataResult@DerivedFromIDs@data[selectMeta]),		
				QualityControlLevelID=getID("QualityControlLevel",getDataResult@Metadata$qualitycontrollevel[selectMeta])
				
				)				
								
	}
}
