setOldClass("request")
setClass(Class = "MovebankLogin",
	 contains="request",
	 validity = function(object){
		 if(nchar(object$headers['user'])==0 || nchar(object$headers['password']==0))
			 return(TRUE)
	 }
	 )
## Browsing Movebank data base
setGeneric("movebankLogin", function(username, password,...) standardGeneric("movebankLogin"))
setMethod(f="movebankLogin", 
	  signature=c(username="character", password="character"), 
	  definition = function(username, password){
		  return(new("MovebankLogin", add_headers(user=username, password=password)))
	  })

setMethod(f="movebankLogin", 
	  signature=c(username="character", password="missing"),
	  definition = function(username, password){
		  pwd<-readline("password:")
		  return(movebankLogin(username=username, password=pwd))
	  })

setMethod(f="movebankLogin", 
	  signature=c(username="missing", password="missing"),
	  definition = function(username, password){
		  user<-readline("username:")
		  return(movebankLogin(username=user))
	  })


##construct URLs and download from Movebank
setGeneric("getMovebank", function(entity_type, login,...) standardGeneric("getMovebank"))
setMethod(f="getMovebank", 
	  signature=c(entity_type="character", login="MovebankLogin"), 
	  definition = function(entity_type, login, ...){
		  tmp <- list(...)
		  url <- paste("https://www.movebank.org/movebank/service/direct-read?entity_type=",entity_type  ,sep="")
		  if(length(tmp)!=0){
			  tmp <- lapply(tmp, paste, collapse='%2C')
			  url <- paste(url, sep="&",paste(names(tmp),tmp, collapse="&", sep="="))
		  }
		  f<-GET(url, config = login)
		  if(entity_type=='event'){cols<-c(location_long='numeric', location_lat='numeric')}else{cols<-NA}
		  
		  data <- read.csv(textConnection(content(f, as='text', encoding = "UTF-8")), colClasses=cols)
		  
		  if(grepl(pattern="You.may.only.download.it.if.you.agree.with.the.terms", x=names(data)[1])) stop("You need a permission to access this data set. Go to www.movebank.org and accept the license terms when downloading the data set (you only have to do this once per data set).")
		  if (grepl(pattern="X.html..head..title.Apache.Tomcat", capture.output(data)[1])) stop("It looks like you are not allowed to download this data set. Or there is a sensor for which no attributes are available.")
		  if (grepl(pattern="are.not.available.for.download", capture.output(data)[1])) stop("You have no permission to download this data set.")
		  if (any(grepl(pattern="503 Service Temporarily Unavailable", unlist(head(data))))) stop("Movebank is (temporarily) unavailable")
		  return(data)
	  })

setMethod(f="getMovebank", 
	  signature=c(entity_type="character", login="missing"), 
	  definition = function(entity_type, login, ...){
		  d<-movebankLogin()
		  getMovebank(entity_type=entity_type, login=d,...)
	  })


#names of the studies
setGeneric("searchMovebankStudies", function(x,login) standardGeneric("searchMovebankStudies"))
setMethod(f="searchMovebankStudies", 
	  signature=c(x="character",login="MovebankLogin"), 
	  definition = function(x,login){  
		  data <- getMovebank("study", login, sort="name", attributes="id%2Cname%2Ci_am_owner%2Ci_can_see_data%2Cthere_are_data_which_i_cannot_see")
		  res <- as.character(data$name)[grepl(x,data$name,useBytes=TRUE)]
		  #    names(res) <- paste("##### Results for ",x,"#####")
		  if(length(res)>0){return(res)}else{"No study matches your search criteria"}
	  })

setMethod(f="searchMovebankStudies", 
	  signature=c(x="character",login="missing"), 
	  definition = function(x,login){  
		  login=movebankLogin()
		  searchMovebankStudies(x=x,login=login)
	  })



#get all study names
setGeneric("getMovebankStudies", function(login) standardGeneric("getMovebankStudies"))
setMethod(f="getMovebankStudies", 
	  signature=c(login="missing"), 
	  definition = function(login){
		  login <- movebankLogin()
		  getMovebankStudies(login=login)
	  })

setMethod(f="getMovebankStudies", 
	  signature=c(login="MovebankLogin"), 
	  definition = function(login){
		  data <- getMovebank("study", login, sort="name", attributes="id%2Cname%2Ci_am_owner%2Ci_can_see_data%2Cthere_are_data_which_i_cannot_see")
		  return(data$name)
	  })


#names of the sensors
setGeneric("getMovebankSensors", function(study, login) standardGeneric("getMovebankSensors"))
setMethod(f="getMovebankSensors", 
	  signature=c(study="ANY",login="missing"), 
	  definition = function(study,login){
		  login <- movebankLogin()
		  getMovebankSensors(study=study, login=login)
	  })
setMethod(f="getMovebankSensors", 
	  signature=c(study="missing",login="missing"), 
	  definition = function(study,login){
		  login <- movebankLogin()
		  getMovebankSensors(login=login)        
	  })

setMethod(f="getMovebankSensors", 
	  signature=c(study="missing",login="MovebankLogin"), 
	  definition = function(study,login){
		  data <- getMovebank("tag_type", login)
		  return(data)
	  })

setMethod(f="getMovebankSensors", 
	  signature=c(study="character",login="MovebankLogin"), 
	  definition = function(study,login){      
		  study <- getMovebankID(study, login)
		  callGeneric()
	  })
setMethod(f="getMovebankSensors", 
	  signature=c(study="numeric",login="MovebankLogin"), 
	  definition = function(study,login){      
		  data <- getMovebank("sensor", login, tag_study_id=study)
		  return(data)
	  })



setGeneric("getMovebankSensorsAttributes", function(study, login) standardGeneric("getMovebankSensorsAttributes"))
setMethod(f="getMovebankSensorsAttributes", 
	  signature=c(study="character",login="MovebankLogin"), 
	  definition = function(study,login){
		  study<-getMovebankID(study, login )
		  callGeneric()
	  })
setMethod(f="getMovebankSensorsAttributes", 
	  signature=c(study="numeric",login="MovebankLogin"), 
	  definition = function(study,login){
		  data <- getMovebank("sensor", login, tag_study_id=study)
		  studySensors <- unique(data$sensor_type_id)
		  data2 <- lapply(studySensors, function(y, login, study) {try(getMovebank("study_attribute", login, study_id=study, sensor_type_id=y), silent=T)} ,login=login, study=study)
		  data2<-data2[(lapply(data2, class))!='try-error']
		  return(as.data.frame(do.call(rbind, data2)))
	  })

###all or a certain ID
setGeneric("getMovebankID", function(study, login) standardGeneric("getMovebankID"))
setMethod(f="getMovebankID", 
	  signature=c(study="character", login="missing"), 
	  definition = function(study, login){
		  login <- movebankLogin()
		  getMovebankID(study=study, login=login)
	  })

setMethod(f="getMovebankID", 
	  signature=c(study="character", login="MovebankLogin"), 
	  definition = function(study=NA, login){
		  data <- getMovebank("study", login, sort="name", attributes="id%2Cname%2Ci_am_owner%2Ci_can_see_data%2Cthere_are_data_which_i_cannot_see")
		  if (is.na(study)) {
			  return(data[ ,c("id","name")])
		  } else {
			  studyNUM <- data[gsub(" ","", data$name)==gsub(" ","", study),c("id")] #get rid of all spaces to avoid miss matching between different spaced words
			  if (length(studyNUM)>1) stop(paste("There was more than one study with the name:",study))
			  return(studyNUM)
		  }
	  })



###retrieving information of a certain study
setGeneric("getMovebankStudy", function(study, login) standardGeneric("getMovebankStudy"))
setMethod(f="getMovebankStudy", 
	  signature=c(study="numeric", login="MovebankLogin"),
	  definition = function(study, login){
		  data <- getMovebank("study", login, id=study)
		  return(data)
	  })
setMethod(f="getMovebankStudy", 
	  signature=c(study="character", login="MovebankLogin"),
	  definition = function(study, login){
		  study<- getMovebankID(study, login)
		  callGeneric()
	  })

setMethod(f="getMovebankStudy", 
	  signature=c(study="ANY", login="missing"),
	  definition = function(study, login){
		  login <- movebankLogin()
		  getMovebankStudy(study=study,login=login)
	  })


##get all animals with their IDs
setGeneric("getMovebankAnimals", function(study, login) standardGeneric("getMovebankAnimals"))
setMethod(f="getMovebankAnimals",
	  c(study="character", login="MovebankLogin"),
	  definition = function(study, login){  
		  study  <- getMovebankID(study,login)
		  callGeneric()
	  })
setMethod(f="getMovebankAnimals",
	  c(study="numeric", login="MovebankLogin"),
	  definition = function(study, login){  
		  tags <- getMovebank(entity_type="sensor", login, tag_study_id=study)
		  animalID <- getMovebank("individual", login, study_id=study, attributes="id%2Clocal_identifier")
		  deploymentID <- getMovebank("deployment", login=login, study_id=study, attributes="individual_id%2Ctag_id%2Cid")
		  #  if (nrow(deploymentID)==0) warning("There are no deployment IDs available!")
		  names(deploymentID)  <- c("individual_id", "tag_id", "deployment_id")   
		  if (nrow(tags)!=0){ 
			  tagdep <- merge.data.frame(x=tags, y=deploymentID, by.x="tag_id", by.y="tag_id", all=TRUE) #skipping tags that have no deployment
			  tagdepid <- merge.data.frame(x=tagdep, y=animalID, by.x="individual_id", by.y="id", all.y=TRUE)[,-3]#skipping the column of the movebank internal tag id
			  colnames(tagdepid) <- c("individual_id", "tag_id", "sensor_type_id", "deployment_id", "animalName")
			  #if (any(apply(deploymentID[,1:2], 2, FUN=duplicated)))##if there are multiple deployments: idData$local_identifier  <-  paste(localID_deploymentID)
			  if (any(duplicated(tagdepid$individual_id)|duplicated(tagdepid$tag_id))){
				  tagdepid$animalName <- paste(tagdepid$animalName, tagdepid$deployment_id, sep="_")
				  names(tagdepid)  <- c("individual_id", "tag_id", "sensor_type_id", "deployment_id", "animalName_deployment")}
			  return(tagdepid) 
		  } else {
			  return(merge.data.frame(x=deploymentID, y=animalID, by.x="individual_id", by.y="id", all.y=TRUE))
		  }
	  })

setMethod(f="getMovebankAnimals",
	  c(study="ANY", login="missing"),
	  definition = function(study, login){
		  login <- movebankLogin()
		  getMovebankAnimals(study=study,login=login)
	  })

setGeneric("getMovebankData", function(study,animalName,login, ...) standardGeneric("getMovebankData"))

setMethod(f="getMovebankData", 
	  signature=c(study="ANY",animalName="ANY", login="missing"),
	  definition = function(study, animalName, login=login, ...){ 
		  login <- movebankLogin()
		  callGeneric()
	  })

setMethod(f="getMovebankData", 
	  signature=c(study="character",animalName="ANY", login="MovebankLogin"),
	  definition = function(study, animalName, login, ...){ 
		  study <- getMovebankID(study, login) 
		  callGeneric()
	  })

setMethod(f="getMovebankData", 
	  signature=c(study="numeric",animalName="missing", login="MovebankLogin"),
	  definition = function(study, animalName, login, ...){ 
		  d<- getMovebank("individual", login=login, study_id=study, attributes=c('id'))$id
		  getMovebankData(study=study, login=login, ..., animalName=d)
	  })
setMethod(f="getMovebankData", 
	  signature=c(study="numeric",animalName="character", login="MovebankLogin"),
	  definition = function(study, animalName, login, ...){ 
		  d<- getMovebank("individual", login=login, study_id=study, attributes=c('id','local_identifier'))
		  animalName<-d[as.character(d$local_identifier)%in%animalName,'id']
		  callGeneric()
	  })
setMethod(f="getMovebankData", 
	  signature=c(study="numeric",animalName="numeric", login="MovebankLogin"),
	  definition = function(study, animalName, login, removeDuplicatedTimestamps=F,...){ 
		  idData <- getMovebank("individual", login=login, study_id=study, id=animalName)         
		  ##which deployments are imporant
		  deploymentID <- getMovebank("deployment", login=login, study_id=study, attributes=c("individual_id","tag_id","id"), individual_id=animalName)
		  ##which track Data are important
		  new <- merge.data.frame(deploymentID, idData, by.x="individual_id", by.y="id")

		  sensorTypes <- getMovebank("tag_type", login=login)
		  rownames(sensorTypes)<-sensorTypes$id
		  locSen <- sensorTypes[as.logical(sensorTypes$is_location_sensor),"id"] #reduce track to location only sensors & only the correct animals
		  
		  attribs <- unique(c(as.character(getMovebankSensorsAttributes(study, login)$short_name),"sensor_type_id","deployment_id",'event_id'))

		  trackDF <- getMovebank("event", login=login, study_id=study, attributes = attribs , deployment_id=new$id, sensor_type_id=locSen)
		  if(nrow(trackDF)==0){
			  stop('No records found for this individual/study combination')
		  }
		  trackDF$timestamp<-as.POSIXct(strptime(as.character(trackDF$timestamp), format = "%Y-%m-%d %H:%M:%OS",tz="UTC"), tz="UTC")
		  if(any(tapply(trackDF$sensor_type_id, trackDF$deployment_id, length)!=1)){# data comes in ordered by sensor but needs to be ordered by timestamp
			  trackDF <- trackDF[ with(trackDF, order(trackDF$deployment_id, timestamp)) , ]  
		  }
		  outliers<- is.na(trackDF$location_long)
		  if('algorithm_marked_outlier'%in%names(trackDF))
			  outliers[trackDF$algorithm_marked_outlier=="true"]<-T
		  if('manually_marked_outlier'%in%names(trackDF)){
			  outliers[trackDF$manually_marked_outlier=="true"]<-T
			  outliers[trackDF$manually_marked_outlier=="false"]<-F
		  }
		  if(all(outliers))
			  stop("There not observed records for this study/individual")

		  spdf<-SpatialPointsDataFrame(trackDF[!outliers,c('location_long','location_lat')], data=trackDF[!outliers,], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"), match.ID=T)
		  if(any(duplicated(new$local_identifier)))
			  new$local_identifier<-paste0(new$local_identifier,'_', new$id)
		  rownames(new)<-validNames(new$local_identifier)
		  id<-paste(format(trackDF$timestamp,"%Y %m %d %H %M %OS4"),trackDF$deployment_id, trackDF$sensor_type_id)
		  if(any(s<-id[outliers]%in%id[!outliers]))
		  {
			  warning("There are timestamps ",sum(s)," in the unused data that are also in the real data, those records are omitted")
			  outliers[outliers][s]<-F
		  }
		  if(any(s<-duplicated(id[outliers])))
		  {
			  warning("There are ",sum(s), " duplicated timestamps in the unused that those will be removed")
			  outliers[outliers][s]<-F
		  }

		  new<-new[new$id %in% unique(spdf$deployment_id),]
		  new<-new[order(new$id),]
		  trackId<-droplevels(factor(spdf$deployment_id, labels=rownames(new), levels=new$id))
		  unUsed<-new('.unUsedRecordsStack', dataUnUsedRecords=trackDF[outliers,],timestampsUnUsedRecords=trackDF$timestamp[outliers], 
			      sensorUnUsedRecords= sensorTypes[as.character(trackDF[outliers,'sensor_type_id']),'name'],
			      trackIdUnUsedRecords=factor(trackDF[outliers,'deployment_id'], labels=rownames(new), levels=new$id))

		  if(any(!(s<-(as.character(unUsed@trackIdUnUsedRecords) %in% levels(trackId)))))
		  {
			  warning('Omiting individual(s) (n=',length(unique(unUsed@trackIdUnUsedRecords[!s])), ') that have only unUsedRecords')
			  unUsed<-unUsed[s,]
		  }
		  unUsed@trackIdUnUsedRecords<-factor(unUsed@trackIdUnUsedRecords, levels=levels(trackId))

		  if(removeDuplicatedTimestamps){
			  message("removeDupilcatedTimestamps was set to true we strongly suggest against it and that the problem is solved before because there is no systematic to which locations are removed. This can for example be done by marking them outlier in movebank.")
			  dupsDf<-(data.frame(format(spdf$timestamp,"%Y %m %d %H %M %OS4"), spdf$sensor_type_id, trackId))
			  dups<-duplicated(dupsDf)
			  spdf<-spdf[!dups,]
			  trackId<-trackId[!dups]
			  warning(sum(dups)," location(s) is/are removed by removeDuplicatedTimestamps")
		  }

		  res<-new("MoveStack", spdf, timestamps=spdf$timestamp, 
			   sensor=sensorTypes[as.character(spdf$sensor_type_id),'name'],
			   unUsed,trackId=trackId,
			   idData=new[as.character(unique(trackId)),])
		  if(length(n.locs(res))==1)
			  res<-as(res, 'Move')
		  return(res)
	  })
