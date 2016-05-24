context('citation')
test_that('citations',{
	a<-move(x=1:10,y=1:10,time=as.POSIXct(1:10, origin='1970-1-1'))
	expect_equal(citations(a), character())
	citations(a)<-"bla"
	expect_equal(citations(a),"bla")
	expect_error(citations(a)<-factor("a"))
	expect_warning(citations(a)<-c("a","a"))
}
)
