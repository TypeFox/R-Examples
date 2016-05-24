context('extent calculations')
test_that('extent calculations',{
data(leroy)
  expect_equal(
		extent(leroy),
		extent(move:::.extcalc(leroy,0))
	)
	expect_equal(
		extent(leroy)*2,
		extent(move:::.extcalc(leroy,0.5))
	)
	expect_equal(
		extent(leroy)*1:2,
		extent(move:::.extcalc(leroy,c(0,0.5)))
	)
})
test_that('extent calculations between bgb and dbbmm',{
	data<-move(x=rep(2:3,9)/5, y=rep(2:3,9)/5, Sys.time()+1:18)# needs to be divisable by 50, and specifc rounding because otherwise one cell gets added in the ud calc of brownian
  expect_message(x<-brownian.bridge.dyn(data, location.error=.1, ext=7, dimSize=50, window.size=11, margin=5),'Computa')
	expect_warning(xd<-dynBGB(data, locErr=.1, ext=7, dimSize=50, windowSize=11, margin=5),'Brownian motion assumed, because no direction could be calculated')
  expect_equal(
		raster(xd),
		raster(x)
	)
	expect_message(x<-brownian.bridge.dyn(data, location.error=.1, ext=7, raster=.020, window.size=11, margin=5),'Computa')
	expect_warning(xd<-dynBGB(data, locErr=.1, ext=7, raster=.020, windowSize=11, margin=5),'Brownian motion assumed, because no direction could be calculated')
	expect_equal(
		raster(xd)
		,
		raster(x)
	)
})
