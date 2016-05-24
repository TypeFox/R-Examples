setGeneric("move", function(x, y, time, data, proj=as.character(NA), ...) standardGeneric("move"))
setMethod(f = "move", 
	  signature = c(x="character",y='missing',time='missing', data='missing', proj='missing'), 
	  definition = function(x, ...){
		  if(!file.exists(x))
			  stop("x should be a file on disk but it cant be found")
		  if(grepl(".zip", x))
		  {
			files<-as.character(unzip(x,list=T)$Name)
		  	if(1!=sum(rd<-(files=='readme.txt'))| length(files)!=2)
				stop('zip file not as expected')
			m<-move(unz(x, files[!rd]),...)
			m@license<-paste(m@license,readLines(con<-unz(x , files[rd])), collapse='\n')
			close(con)
			return(m)
		  }else{
			  return(move(file(x),...))
		  }

	  })
setMethod(f = "move", 
	  signature = c(x="connection",y='missing',time='missing', data='missing', proj='missing'), 
	  definition = function(x, removeDuplicatedTimestamps=F,...){
		  if(version$major=='3' & version$minor=='1.0'){# exception to type convert back then that got reverted, setting colclasses causes problems with quoted downloads from envData
		  df <- read.csv(x, header=TRUE, sep=",", dec=".", stringsAsFactors=T, colClasses=c(location.long='numeric', location.lat='numeric'))
		  }else{
		  df <- read.csv(x, header=TRUE, sep=",", dec=".", stringsAsFactors=T)
		  }
		  if (!all(c("timestamp", 
			     "location.long",  
			     "location.lat", 
#			     "study.timezone", 
#			     "study.local.timestamp", 
			     "sensor.type", 
			     "individual.local.identifier", 
			     "individual.taxon.canonical.name")%in%colnames(df))){
			  stop("The entered file does not seem to be from Movebank. Please use the alternative import function.")
		  }

		  if(any(dups<-duplicated( do.call('paste',c(df[duplicated(df$timestamp)|duplicated(df$timestamp, fromLast=T),names(df)!="event.id"], list(sep="__")))))){
			  #first find atleast the ones where the timestamp (factor) is duplicated
			  warning("Exact duplicate records removed (n=",sum(dups),") (movebank allows them but the move package can't deal with them)")
			  df <- df[!duplicated( do.call('paste',c(df[,names(df)!="event.id"], list(sep="__")))),]
			  # cant use dups here since it that uses the optimization of only looking at timestamps first
		  }

		  df$timestamp <- as.POSIXct(strptime(as.character(df$timestamp), format = "%Y-%m-%d %H:%M:%OS",tz="UTC"), tz="UTC") # need to make character out of it to ensure milli seconds are considerd
		  if("study.local.timestamp"%in% names(df)){
		   df$study.local.timestamp <- as.POSIXct(strptime(df$study.local.timestamp, format="%Y-%m-%d %H:%M:%OS"))
		  }

		  if(any(tapply(df$sensor.type, df$individual.local.identifier, length)!=1)){
			  df <- df[with(df, order(df$individual.local.identifier, timestamp)), ]  
		  }
		  df$individual.local.identifier<-as.factor( df$individual.local.identifier)
		  levels(df$individual.local.identifier) <- validNames(levels(factor(df$individual.local.identifier))) #changing names to 'goodNames' skipping spaces

		  if("visible" %in% colnames(df))
		  {
			  v<-df$visible=='false'
		  }else{
			  v<-F
		  }
		  unUsed<-is.na(df$location.long)|is.na(df$location.lat)|v| is.na(df$individual.local.identifier)| df$individual.local.identifier==''
		  sensor<-df$sensor.type
		  timestamps<-df$timestamp
		  individual.local.identifier<-df$individual.local.identifier
		  uniquePerID <- unlist(lapply(df,  function(x,y){all(tapply(x,y,function(x){length(unique(x))})==1)}, y=factor(df$individual.local.identifier)))
		  idData <- subset(df, select=names(uniquePerID[uniquePerID]), !duplicated(df$individual.local.identifier))
		  rownames(idData)<-idData$individual.local.identifier
		  df<-df[,!(names(df)%in%unique(c('sensor.type','timestamps', colnames(idData))))]
		  unUsedDf<-df[unUsed,]
		  unUsedRecords<-new('.unUsedRecords', dataUnUsedRecords=unUsedDf, timestampsUnUsedRecords=timestamps[unUsed], sensorUnUsedRecords=sensor[unUsed])
		  if(stk<-nrow(idData)!=1){
			  unUsedRecords<-new('.unUsedRecordsStack', unUsedRecords, trackIdUnUsedRecords=individual.local.identifier[unUsed])
		  }
		  if(any(s<-!(rownames(idData) %in% as.character(unique(individual.local.identifier[!unUsed])))))
		  {
			  warning('omiting ',sum(s),' individual(s) because they do not have observation data')
			  unUsedRecords<-unUsedRecords[(as.character(unUsedRecords@trackIdUnUsedRecords) %in% rownames(idData)[!s]),]
			  unUsedRecords@trackIdUnUsedRecords<-factor(unUsedRecords@trackIdUnUsedRecords)
			  idData<-idData[!s,]
		  }

		  df<-df[!unUsed,]
		  coordinates(df)<- ~location.long+location.lat
		  proj4string(df)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

		  track<-new('.MoveTrack', df, timestamps=timestamps[!unUsed], sensor=sensor[!unUsed], idData=idData)
	          individual.local.identifier<-factor(individual.local.identifier[!unUsed])
		  if(removeDuplicatedTimestamps){
			message("removeDupilcatedTimestamps was set to true we strongly suggest against it and that the problem is solved before because there is no systematic to which locations are removed. This can for example be done by marking them outlier in movebank.")
			dupsDf<-(data.frame(format(timestamps(track),"%Y %m %d %H %M %OS4"), track@sensor))
			dupsDfUn<-(data.frame(format(timestamps(unUsedRecords),"%Y %m %d %H %M %OS4"), unUsedRecords@sensorUnUsedRecords))
			
			if(stk){
				dupsDf<-data.frame(id=individual.local.identifier, dupsDf)
				dupsDfUn<-data.frame(id=unUsedRecords@trackIdUnUsedRecords, dupsDfUn)
			}
			dups<-duplicated(dupsDf)
			track<-track[!dups,]
			individual.local.identifier<-individual.local.identifier[!dups]
			warning(sum(dups)," location(s) is/are removed by removeDuplicatedTimestamps")
			unUsedRecords<-unUsedRecords[inNormRe<-!(apply(dupsDfUn, 1, paste, collapse='_')%in%apply(dupsDf, 1, paste, collapse='_')),T]
      warning(sum(!inNormRe), " location(s) is/are removed by removeDuplicatedTimestamps from the un used records")
		  }

		  if(stk){
			  return(new('MoveStack', track, unUsedRecords, trackId=individual.local.identifier))
		  }else{
			  return(new('Move',track, unUsedRecords))
		  }
	  }
	  )

#if non-Movebank data are used, table is new defined 
setMethod(f="move",
	  signature=c(x="numeric", y="numeric", time="POSIXct", data="missing", proj="ANY"),
	  definition = function(x,y,time,data,proj, ...){
		  data<-data.frame(x,y,time)
		  move(x=x,y=y,time=time,proj=proj,data=data,...)
	  }
)

setMethod(f="move",
          signature=c(x="numeric", y="numeric", time="POSIXct", data="data.frame", proj="character"),
          definition = function(x,y,time,data,proj, ...){
            move(x=x,y=y,time=time,proj=CRS(proj),data=data,...)
          }
)

setMethod(f="move",
          signature=c(x="numeric", y="numeric", time="POSIXct", data="data.frame", proj="missing"),
          definition = function(x,y,time,data,proj, ...){
            move(x=x,y=y,time=time,proj=CRS(),data=data,...)
          }
)

setMethod(f="move",
	  signature=c(x="numeric", y="numeric", time="POSIXct", data="data.frame", proj="CRS"),
	  definition = function(x,y,time,data,proj,sensor='unknown',animal='unnamed', ...){
		  data$location.long <- x
		  data$location.lat <- y
		  data$timestamp <- time
		  data$individual.local.identifier <- animal
		  data$sensor <- factor(sensor)
      
		  if(any(is.na(data$location.long)|is.na(data$location.lat))){ 
		    warning("There were NA locations detected and omitted. Currently they are not stored in unusedrecords")
		    data <- data[!(is.na(data$location.long)|is.na(data$location.lat)), ]
		  }
		  data$individual.local.identifier<-factor(data$individual.local.identifier, levels = unique(data$individual.local.identifier))
		  levels(data$individual.local.identifier) <- validNames(levels(factor(data$individual.local.identifier))) #changing names to 'goodNames' skipping spaces
		  
		  if(length(unique(data$individual.local.identifier))>1 & any(unique(as.character(data$individual.local.identifier))==""))
		  {
		    warning("Omitting locations that have an empty local identifier (n=",sum(tmp<-as.character(data$individual.local.identifier)==""),"). Most likely the tag was not deployed") 
		    data <- data[!tmp,]
		    data$individual.local.identifier <- factor(data$individual.local.identifier)
		  }
		  ids <- as.list(as.character(unique(data$individual.local.identifier)))
		  uniquePerID<-unlist(lapply(colnames(data),  function(x,y,data){nrow(unique(data[,c(x,'individual.local.identifier')]))},data=data))==sum(!duplicated(data$individual.local.identifier))
		  names(uniquePerID)<-colnames(data)
		  uniquePerID["sensor"] <- FALSE
		  idData <- subset(data, select=names(uniquePerID[uniquePerID]), !duplicated(data$individual.local.identifier))
		  
		  if(length(names(idData))!=1)# dont shorten it because we need something
		    idData<-subset(idData, select=names(idData)!="individual.local.identifier")
		  
		  if(length(unique(idData$citation))>0) 
		  {
		    if(length(unique(idData$citation))>1) 
		      warning("There were more than one citation for this study found! Only using the first.")
        citations <- as.character(unique(idData$citation))[1]
		  } else {
        citations <- character()
		  }
		  
		  rownames(idData) <- unique(data$individual.local.identifier)
		  auxData <- data.frame(data[names(data)[!names(data)%in%c("location.lat", "location.long","timestamp", colnames(idData))]])
		  
		  if (ncol(auxData)==0) auxData <- data.frame(auxData, empty=NA)
	  
		  spdf <- SpatialPointsDataFrame(
		    coords = cbind(data$location.long,data$location.lat),
		    data = auxData, 
		    proj4string = proj,
		    match.ID = TRUE)
			  
		  if (length(ids)==1){
		    return(new("Move", 
		               timestamps = data$timestamp, 
		               sensor = data$sensor,
		               sensorUnUsedRecords=factor(levels=levels(data$sensor)),
		               spdf, 
		               citation = citations,
		               idData = idData
		    ))
		  } else {
		    trackId<-factor(data$individual.local.identifier)
		    return(new("MoveStack", 
		               spdf, 
		               idData = idData,
		               sensor = data$sensor,
		               sensorUnUsedRecords=factor(levels=levels(data$sensor)),
		               timestamps = data$timestamp, 
		               citation = citations,
		               trackId = trackId,
		               trackIdUnUsedRecords=factor(levels=levels(trackId))))
		  }
	  }
	  )
setMethod(
  f = "move",
  signature = c(
    x = "ltraj",y = 'missing',time = 'missing', data = 'missing', proj = 'missing'
  ),
  definition = function(x, ...) {
    if (length(x) == 1) {
      return(as(x,"Move"))
    }else{
      return(as(x,"MoveStack"))
    }
  }
)