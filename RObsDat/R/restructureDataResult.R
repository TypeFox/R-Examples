restructureDataResult <- function(to.ret, value.numeric=TRUE, tz=c("global", "UTC", "GMT", "0", "local")){
  
  
  tz <- match.arg(tz)
  if(NROW(to.ret)>0) {
    #to be compatible with postgres
    colnames(to.ret) <- tolower(colnames(to.ret))
    #to be able to compare also NA
    to.ret[is.na(to.ret)] <- -9999999999
    #sort according to data types and create to xts objects, one numeric, one string
    data.types <- sapply(to.ret, is.numeric)
    the.numerics <- which(data.types)
    the.char <- which(!data.types)
    
    
    if(tz=="local"){
      index.col <- which(names(to.ret)=="localdatetime")
      tzname <- to.ret$utcoffset
      todo("Find a solution to timezone representation. 
           The mapping Offset (as suggested in ODM1.1) to
           a specific timezone is not unique, thus currently
           the information get's lost from the database")
    } else {
      index.col <- which(names(to.ret)=="datetimeutc")
      tzname <- "GMT"
    }
    the.char <- the.char[the.char!=index.col]

    time.order.by <- chr2date(to.ret[,index.col], tz=tzname)
    time.order.by.unique <- chr2date(sort(unique(to.ret[,index.col])), tz=tzname) # chr2date macht es zum POSIX
   
    metadata.id <- rep(NA, NROW(to.ret)) # create empty vector
	non.metadata.columns <- colnames(to.ret) %in% tolower(c("ValueID", "DataValue", "LocalDateTime", "DateTimeUTC", "UTCOffset", "DerivedFromID", "VersionID"))
    metadata <- unique(to.ret[,!non.metadata.columns])
    metadata.plain <- id2name(metadata)
	metadata.id <- rep(NROW(metadata.plain), NROW(to.ret))
    
    #convert NA back
    to.ret[to.ret == -9999999999] <- NA
    
    #create the spatialPoint
      siteData = getMetadata("Site", ID=unique(to.ret$siteid))
      sp = cbind(siteData$Latitude, siteData$Longitude)
      row.names(sp) = siteData$Name

    sp = SpatialPoints(sp)

    # create Spacetime-Columnnames (consists of variablenname and MetadataID)
    variablenname = unique(metadata.plain$variable)
    st_colnames <- NULL
    	
    for(variable in variablenname){
        st_colnames <- c(st_colnames,variable)
    }
	MetaId = unique(metadata.id)
	
	# create STFDF-Objects-Attachments
    stfdfMain <- createST(sp=sp, siteID=siteData$ID, location=to.ret$siteid, timeobject=time.order.by.unique, timelong=time.order.by, variables = st_colnames, thedata=to.ret$datavalue)
	stfdfValueIDs <- createST(sp=sp, siteID=siteData$ID, location=to.ret$siteid, timeobject=time.order.by.unique, timelong=time.order.by, variables = st_colnames, thedata=to.ret$valueid)
    stfdfDerivedFromIDs <- createST(sp=sp, siteID=siteData$ID, location=to.ret$siteid, timeobject=time.order.by.unique, timelong=time.order.by, variables = st_colnames, thedata=to.ret$derivedfromid)
    stfdfMetaRelation <- createST(sp=sp, siteID=siteData$ID, location=to.ret$siteid, timeobject=time.order.by.unique, timelong=time.order.by, variables = st_colnames, thedata=metadata.id)
    
	
	# Content Metadata
    metadataTable <- cbind(MetaId, metadata.plain)

	 # create the hole objects in an inherited_stfdf
	inherited_stfdf_object = inherited_stfdf(sp = stfdfMain@sp, time = stfdfMain@time, data = stfdfMain@data, endtime = stfdfMain@endTime,
									ValueIDs=stfdfValueIDs,
									DerivedFromIDs=	stfdfDerivedFromIDs,
									MetadataRel = stfdfMetaRelation,
									Metadata=metadataTable)    
	return(inherited_stfdf_object)
  }
}