context('spTransform')
test_that('spTransfrom',{
data(leroy)
  data2<-spTransform(x=d<-spTransform(leroy,center=T), CRSobj=(proj4string(leroy)))
	leroy@coords.nrs=numeric(0)# sptransfrom does this for some reason in line 123 of project.R
	expect_equal(leroy,data2)
	expect_error(spTransform(d, center=T))# somethign not long lat cant go to aeqd
})
