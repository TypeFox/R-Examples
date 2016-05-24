context('move ade')
test_that('moveAde',
{
data(ricky)
data<-ricky
	dataSp<-spTransform(data, center=T)
	dataLtraj<-as(dataSp, 'ltraj')
	dataBack<-as(dataLtraj, 'Move')
	expect_equivalent(coordinates(dataSp), coordinates(dataBack))
	expect_equal(timestamps(dataSp), timestamps(dataBack))
	spLtraj<-move2ade(dataSp)
	#suppressMessages(require(adehabitatLT))
	a<-as(adehabitatLT::simm.crw(1:100),'Move')
	expect_true(validObject(a))
	aa<-as(a, 'ltraj')
	ma<-move(aa)
	ma$sensor<-NULL
	idData(ma)<-idData(a)
	expect_equal(a,ma)
	expect_true(validObject(a))
	a<-as(adehabitatLT::simm.crw(1:100, id=gl(25,4)),'MoveStack')
	expect_true(validObject(a))
	aa<-as(a, 'ltraj')
	expect_true(validObject(aa))
	ma<-move(aa)
	ma$sensor<-NULL
	idData(ma)<-idData(a)
	ma@trackId<-trackId(a)
	levels(ma@trackIdUnUsedRecords)<-levels(trackId(unUsedRecords(a)))
	expect_equal(a,ma)
}
)
