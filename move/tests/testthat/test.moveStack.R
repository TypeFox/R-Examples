context('moveStack')
test_that('moveStack',
{
	at<-as.POSIXct(1:10, origin='1970-1-1',tz="GMT")
	bt<-as.POSIXct(5+1:10, origin='1970-1-1',tz="GMT")       
	attr(bt,'tzone')<-NULL
	attr(at,'tzone')<-NULL
	a<-move(x=1:10,y=1:10,time=at,proj=CRS('+proj=longlat +ellps=WGS84'))
	b<-move(x=1:10,y=1:10,time=bt,proj=CRS('+proj=longlat +ellps=WGS84'), 
		animal="a")
	expect_identical(coordinates(a), coordinates(b))
	#	DEACTIVATED("Need to look what we want here")
	bb<-split(d<-moveStack(list(a,b)))
  expect_true(validObject(d))
	aa<-list(unnamed=a,a=b)
	row.names(aa[[2]])<-1:10
	row.names(bb[[2]])<-1:10
	bb<-lapply(bb, function(x){x@idData$individual.local.identifier<-factor(x@idData$individual.local.identifier); return(x)})
	expect_equal(bb,aa)# one problem seems to be moveStack does not deal with missed fixes, the other the rownames of the data frame
	row.names(d@idData)<-sub('a','A A', row.names(d@idData))
	expect_error(new('MoveStack', d , trackId=factor(sub('a','A A', as.character(d@trackId))), idData=d@idData))# track ids are no good names
	expect_error(validObject(d))#validity check needs to fail because of changed rownames

	a<-move(x=1:10,y=1:10,time=at,proj=CRS('+proj=longlat +ellps=WGS84'),animal="AAA")
	a2<-move(x=1:10,y=1:10,time=at,proj=CRS('+proj=longlat +ellps=WGS84'),animal="AAA")  
	expect_warning(uuu<-moveStack(list(a,a2)),'Detected duplicated names. Renamed the duplicated individuals accordingly.')# warn about duplicate ids
	expect_equal(n.indiv(uuu), 2)

	projection(a2) <- CRS(as.character(NA))
	expect_error(moveStack(list(a,a2,a2)))
	expect_error(moveStack(list(a,a,a2)))
	expect_error(moveStack(list(a2,a2,a)))
	expect_error(moveStack(list(a,a2,a)))
	expect_error(moveStack(list(a2,a,a2)))


	m<-lapply(1:5, function(x){
	  m<-move(rnorm(5), rnorm(5), Sys.time()+1:5)
	  idData(m)<-data.frame(groupID=sample(letters,1), groupDay=round(rnorm(1)))
	  m
	  })
	expect_warning(mm<-moveStack(m))
	expect_true(all(grepl('^............: ',capture.output(print(mm)))))
  expect_equal(kk<-capture.output(print(mm)), kl<-capture.output(show(mm)))
	expect_equal(unique(timeLag(mm,units='mins')), list(c(1,1,1,1)/60))
  expect_identical(names(idData(mm)),c('groupID','groupDay'))
}
)