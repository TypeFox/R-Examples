context('subset')
test_that('subSet',
{
	a<-move(x=1:10,y=1:10,time=as.POSIXct(1:10, origin='1970-1-1'),proj=CRS('+proj=longlat +ellps=WGS84'))
	b<-move(x=1:10,y=1:10,time=as.POSIXct(1:10, origin='1970-1-1'),proj=CRS('+proj=longlat +ellps=WGS84'),animal='b') 
	s<-moveStack(list(a,j=b))
	sl<-moveStack(list(l=a,j=b, k=a,p=b ))
	bb<-s[['j']]
	names(bb@sensor)<-NULL
	names(bb@timestamps)<-NULL
	rownames(bb@data)<-(1:10)
	rownames(bb@coords)<-NULL
	bb@idData<-b@idData
	#need to make sure the results are more equal
	attr(b@timestamps,"tzone")<-NULL
	expect_equal(b, bb)
	expect_equal(s[['j']], split(s)[['j']]) 
	expect_equal(s[[1]], s[[c(T,F)]])
	expect_equal(s[[2]], s[['j']])
	ss<-s[[T]]
#	ss@bbox<-bbox(s)
#	attributes(ss@coords)$dimnames[[2]]<- attributes(s@coords)$dimnames[[2]]
	expect_equal(s,ss)
	expect_equal(sl[[idData(sl,T,'individual.local.identifier')=='b']], sl[[c(2,4)]])
	expect_equal(sl[[c('j','p')]], sl[[c(2,4)]])
#	for(i in slotNames(b))
#	{
#		message(i)
#		try(expect_equal(slot(b,i), slot(bb,i)))
#	}
	

})
