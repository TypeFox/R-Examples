context('move')
test_that('move',
	  {
		  expect_error(move())
		  #expect_error(move(x=1:10, y=1:10, time=as.POSIXct(1:10, origin="1970-1-1"))) #we now allow projection to be NA
		  data <- move(fileUri<-system.file("extdata","leroy.csv.gz",package="move"))
		  expect_equal(n.indiv(data),1)
		  cvs<-na.omit(read.csv(fileUri)[,c('location.long','location.lat')])
		  expect_equivalent(coordinates(data),as.matrix(cvs))
		  data2<-move(file(fileUri))
		  data@dateCreation<- data2@dateCreation
		  expect_equal(data, data2)
		  expect_message(data3<-move(file(fileUri), removeDuplicatedTimestamps=T),'removeDupilcatedTimestamps was set to true we strongly suggest against it and that the problem is solved before because there is no systematic to which locations are removed. This can for example be done by marking them outlier in movebank'
				 )
		  data3@dateCreation<- data2@dateCreation
		  expect_equal(data, data3)
		  expect_warning( suppressMessages(data3<-move(file(fileUri), removeDuplicatedTimestamps=T)))
		  fileR<-system.file("extdata","leroy.csv.gz",package="move")
		  if((class(p<-pipe(''))!='pipeWin32')[1])
		  {  if(class(try(move(pipe(paste0('zcat ',fileR, " | sed ''")))))=='Move'){
			  dataR<-move(pipe(paste0('zcat ',fileR)))
			  f<-c( "\"\",\"event.id\",\"timestamp\",\"location.long\",\"location.lat\",\"behavioural.classification\",\"eobs.battery.voltage\",\"eobs.fix.battery.voltage\",\"eobs.horizontal.accuracy.estimate\",\"eobs.key.bin.checksum\",\"eobs.speed.accuracy.estimate\",\"eobs.start.timestamp\",\"eobs.status\",\"eobs.temperature\",\"eobs.type.of.fix\",\"eobs.used.time.to.get.fix\",\"ground.speed\",\"heading\",\"height.above.ellipsoid\",\"manually.marked.outlier\",\"visible\",\"sensor.type\",\"individual.taxon.canonical.name\",\"tag.local.identifier\",\"individual.local.identifier\",\"study.name\",\"utm.easting\",\"utm.northing\",\"utm.zone\",\"study.timezone\",\"study.local.timestamp\"",
			       "\"\",\"event.id\",\"timestamp\",\"location.long\",\"location.lat\",\"behavioural.classification\",\"eobs.battery.voltage\",\"eobs.fix.battery.voltage\",\"eobs.horizontal.accuracy.estimate\",\"eobs.key.bin.checksum\",\"eobs.speed.accuracy.estimate\",\"eobs.start.timestamp\",\"eobs.status\",\"eobs.temperature\",\"eobs.type.of.fix\",\"eobs.used.time.to.get.fix\",\"ground.speed\",\"heading\",\"height.above.ellipsoid\",\"manually.marked.outlier\",\"visible\",\"sensor.type\",\"individual.taxon.canonical.name\",\"tag.local.identifier\",\"individual.local.identifier\",\"study.name\",\"utm.easting\",\"utm.northing\",\"utm.zone\",\"study.timezone\",\"study.local.timestamp\"",
			       "\"1\",44528397,\"2010-02-09 16:07:21.000\",NA,NA,NA,3686,3468,NA,2210166656,NA,\"2010-02-09 16:05:21.000\",\"D\",27,0,120,NA,NA,NA,\"\",\"true\",\"gps\",\"Martes pennanti\",1016,\"Ricky T\",\"Urban fisher GPS tracking\",NA,NA,\"\",\"Eastern Standard Time\",\"2010-02-09 11:07:21.000\"",                                                                                                                                                                                                                                                                                                                                                                                                    
			       "\"1\",44528397,\"2010-02-09 16:07:21.000\",NA,NA,NA,3686,3468,NA,2210166656,NA,\"2010-02-09 16:05:21.000\",\"D\",27,0,120,NA,NA,NA,\"\",\"true\",\"gps\",\"Martes pennanti\",1016,\"Ricky T\",\"Urban fisher GPS tracking\",NA,NA,\"\",\"Eastern Standard Time\",\"2010-02-09 11:07:21.000\"",                                                                                                                                                                                                                                                                                                                                                                                                    
			       "\"2\",44528398,\"2010-02-09 17:01:23.000\",-73.9042594,42.8418856,NA,3662,3457,13.57,2236237607,1.64,\"2010-02-09 17:00:00.000\",\"A\",3,3,83,0.12,292.95,80.3,\"\",\"true\",\"gps\",\"Martes pennanti\",1016,\"Ricky T\",\"Urban fisher GPS tracking\",589541.097246732,4743838.94122749,\"18N\",\"Eastern Standard Time\",\"2010-02-09 12:01:23.000\"",                                                                                                                                                                                                                                                                                                                                         
			       "\"2\",44528398,\"2010-02-09 17:01:23.000\",-73.9042594,42.8418856,NA,3662,3457,13.57,2236237607,1.64,\"2010-02-09 17:00:00.000\",\"A\",3,3,83,0.12,292.95,80.3,\"\",\"true\",\"gps\",\"Martes pennanti\",1016,\"Ricky T\",\"Urban fisher GPS tracking\",589541.097246732,4743838.94122749,\"18N\",\"Eastern Standard Time\",\"2010-02-09 12:01:23.000\"")
			  ff<-textConnection(f)
			  class(ff)<-'connection'
			  expect_error(move(ff))# location is character due to duplicated header
			  expect_warning(move(pipe(paste0('zcat ',fileR,"|  sed '1!p;'"))),'Exact duplicate records removed')
			  dataR2<-suppressWarnings(move(pipe(paste0('zcat ',fileR,"|  sed '1!p;'"))))
			  dataR@dateCreation<- dataR2@dateCreation
			  rownames(dataR@data)<- as.character(rownames(dataR2@data))
			  rownames(dataR2@data)<- as.character(rownames(dataR2@data))
			  rownames(dataR@dataUnUsedRecords)<- as.character(rownames(dataR2@dataUnUsedRecords))
			  rownames(dataR2@dataUnUsedRecords)<- as.character(rownames(dataR2@dataUnUsedRecords))
			  rownames(dataR2@coords)<- rownames(dataR@coords)
			  expect_equal(dataR, dataR2)
			  expect_error(dataR3<-(move(pipe(paste0('zcat ',fileR,"| sed '48p; 48s/Martes/Maggg/'")))))
			  expect_error(move(pipe(paste0('zcat ',fileR,"| sed '48p; 48s/\"A\"/\"B\"/'")), removeDuplicatedTimestamps=F))
			  expect_error(move(pipe(paste0('zcat ',fileR,"| sed '56p; 56s/8986618/8886618/'")), removeDuplicatedTimestamps=F))
			  expect_equal(n.locs(move(pipe(paste0('zcat ',fileR,"| sed '55p; 55s/49.001/49.002/'")), removeDuplicatedTimestamps=F)), n.locs(dataR2)+1)
			  expect_warning(suppressMessages(dataR3<-move(pipe(paste0('zcat ',fileR,"| sed '53p; 53s/10.49/10.48/'")), removeDuplicatedTimestamps=T)),
					 ". location.s. is/are removed by removeDuplicatedTimestamps")
			  dataR@dateCreation<- dataR3@dateCreation
			  rownames(dataR@data)<- as.character(rownames(dataR3@data))
			  rownames(dataR3@data)<- as.character(rownames(dataR3@data))
			  rownames(dataR@dataUnUsedRecords)<- as.character(rownames(dataR3@dataUnUsedRecords))
			  rownames(dataR3@dataUnUsedRecords)<- as.character(rownames(dataR3@dataUnUsedRecords))
			  rownames(dataR3@coords)<- rownames(dataR@coords)
			  expect_equal(dataR, dataR3)
}
		  }
		  close(p)
		  expect_is(move(1:2,1:2, Sys.time()+1:2,animal=c('b','a')),'MoveStack')
		  expect_is(move(1:2,1:2, Sys.time()+1:2,animal=c(2:1)),'MoveStack')
		  tt<-Sys.time()+1:2
		  expect_equal(  move(1:2,1:2,tt ,proj="+proj=longlat +ellps=WGS84") , move(1:2,1:2, tt,proj=CRS("+proj=longlat +ellps=WGS84"))  )
		  expect_equal(  move(1:2,1:2,tt ,data=data.frame(3:4),proj="+proj=longlat +ellps=WGS84") , move(1:2,1:2, tt,data=data.frame(3:4),proj=CRS("+proj=longlat +ellps=WGS84"))  )

	  }
)
