context('idData')
test_that('idData',{
data(leroy)
  expect_equal(class(idData(leroy)),'data.frame')
	expect_equal(idData(leroy), idData(leroy,1))
	expect_equal(idData(leroy), idData(leroy,T))
	expect_equal(class(idData(leroy,T,'sensor.type')),'factor')
	expect_equal(class(idData(leroy,T,'sensor.type', drop=F)),'data.frame')
	dataOld<-leroy
	idData(leroy,1)<-idData(leroy)
	expect_equal(leroy, dataOld)
})
