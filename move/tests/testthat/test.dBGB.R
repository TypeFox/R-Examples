# echo '`require(testthat); require(raster);require(geosphere); suppressMessages(require(rgdal)); t<-sapply(list.files("R", full.names=TRUE), source); test_dir("tests/testthat")' | R -q --vanilla
context("dynBGB")
test_that("dbgb vs dbbmm variance",{
	 vBBMM<-move:::brownian.motion.variance(x<-c(1,1.5,3,3.5,5,5.5,7,7.5,9), y<-c(0,.5,0,.5,0,.5,0,.5,0), time.lag=rep(1,9), location.error=rep(.1,9))
	 yy<-move(x,y, Sys.time()+1:9*60)
	 vBGB<-move:::BGBvar(yy, locErr=rep(.1,9))
	 expect_equal(unlist(vBBMM$cll), unlist(vBGB[3]), check.names=F, tolerance=1e-7)
	 expect_equal(vBBMM$BMvar, prod(vBGB[1:2]), tolerance=2e-5)
})
test_that('dbgb error handling in relation to projections',{
data(leroy)
  # 	expect_warning(
    expect_error(
		     suppressWarnings(dynBGB(leroy[1:40,], locErr=.0005, raster=.00150, ext=17.3, margin=15, windowSize=31))
		     ,'You can not use longitude latitude projection for this function. To transform your coordinates use the spTransform function')
# 		     ,'Optimized to zero')
	r <- raster(nrows=100, ncols=100, xmn=0, xmx=10)
	expect_error(
		     dynBGB(spTransform(leroy[1:40,], center=T), locErr=.0005, raster=r, ext=17.3, margin=15, windowSize=31),'The projection of the raster and the Move object are not equal'
		     )# equal projection
})
test_that("dbgb vs dbbmm",{
	 x<-12:42
	 xx<-floor(x/2)+(y<-((x/2)%%1)/6)
	 m<-move(x+y,y, Sys.time()+x)
#	 l<-(x+max(x))/(max(x)*5)
	 l<-rep(.1, length(x))
	 r<-raster(extent(m)+c(-22,1.5))
	 res(r)<-.025
	 tss<-.00110523/60
	 bgb<-dynBGB(m,r ,locErr=l, windowSize=23, margin=9,  timeStep=(tss<-.00110523/60) ) 
	 suppressMessages(bbmm<-brownian.bridge.dyn(m,raster=bgb ,location.error=l, window.size=23, margin=9 , time.step=tss))
#	 par(mfrow=c(2:1));plot(bgb); lines(m); points(m);plot(bbmm); lines(m); points(m)
#	 plot(bgb-bbmm); lines(m)
#	 plot((bgb-bbmm)/(bgb+bbmm), zlim=c(-.1,.1)); lines(m)
#	 plot( bgb@var@orthSd*bgb@var@paraSd, bbmm@DBMvar@means); abline(0,1)
	 bbmm@method<-'dynBGB'
	 expect_equal(as(bgb, '.UD'), as(bbmm, '.UD'), tolerance=tss*60)
})
test_that("deltaParaOrth",{
	  expect_equivalent(move:::deltaParaOrth(cbind(0,1), cbind(1,2), cbind(0,1)), cbind(0,0))
	  expect_equivalent(move:::deltaParaOrth(cbind(0,1), cbind(1,2), cbind(0,0)), cbind(-sqrt(.5),sqrt(.5)))
	  expect_equivalent(move:::deltaParaOrth(cbind(0,0), cbind(2,2), cbind(1,1)), cbind(sqrt(2),0))
	  suppressWarnings(expect_equivalent(move:::deltaParaOrth(cbind(0,0), cbind(0,0), cbind(1,1)), cbind(1,1)))
	  expect_warning(move:::deltaParaOrth(cbind(0,0), cbind(0,0), cbind(1,1)), "Brownian motion assumed, because no direction could be calculated")
	  
})
test_that("dyn bgb basics",{
data(leroy)
  dataC<-spTransform(leroy, center=T)
	  resUd<-5.3
	  ud<-dynBGB(dataC[1:45,], windowSize=31, margin=15, locErr=4, raster=resUd, ext=9)
	  ud2<-dynBGB(dynBGBvariance(dataC[1:45,], windowSize=31, margin=15, locErr=l<-rep(4, n.locs(dataC))), raster=resUd, ext=9, locErr=l)
	  expect_is(ud, 'dynBGB')
	  expect_is(ud, '.UD')
	  expect_equal(res(ud), resUd[c(1,1)])
	  expect_equal(proj4string(ud), proj4string(dataC))
	  expect_equal(ud2, ud)
})
test_that('work with time step with varying loc err',{
	  l<-dynBGB(m<-new('dBGBvariance',move(0:1,0:1, Sys.time()+0:1), margin=1, windowSize=1, segInterest=c(T,F), paraSd=1:0/2, orthSd=1:0/2, nEstim=2:3), raster=.004, ext=2, locErr=(le<-c(.01,.03)), timeStep=.4/60)
	  expect_equal(sum(values(p<-crop(l,e<-extent(.5,.5,.5,.5)+.45))) ,1/3, tolerance=1.4e-7)
	  a<-.5
	  s<-sqrt((1/60)*a*(1-a)*m@paraSd[1]* m@orthSd[1]+a^2*le[2]^2+(1-a)^2*le[1]^2)
	  pp<-calc(rasterFromXYZ(rasterToPoints(p)[,c(1,2,1,2)]), function(x){prod(dnorm(x, .5, sd=s))* prod(res(p))/3})
	  expect_equal(values(p), values(pp), tolerance=1e-5)
	  expect_less_than(sum(values(abs(p-pp))),2e-6)
	  expect_equal(sum(values(p<-crop(l,e<-extent(.1,.1,.1,.1)+.4))) ,1/3, tolerance=1e-7)
	  a<-.1
	  s<-sqrt((1/60)*a*(1-a)*m@paraSd[1]* m@orthSd[1]+a^2*le[2]^2+(1-a)^2*le[1]^2)
	  pp<-calc(rasterFromXYZ(rasterToPoints(p)[,c(1,2,1,2)]), function(x){prod(dnorm(x, .1, sd=s))* prod(res(p))/3})
	  expect_equal(values(p), values(pp), tolerance=1e-5)
	  expect_less_than(sum(values(abs(p-pp))),2.2e-6)
})
test_that('work with time step var orth para',{
	  r<-raster(extent(c(-.5,1.54,-1.123,1)))
	  res(r)<-.0053
	  l<-dynBGB(m<-new('dBGBvariance',move(0:2/2,c(0,0,0), Sys.time()+0:2), margin=1, windowSize=1, segInterest=c(T,T,F), paraSd=2:0/2, orthSd=1:3/3, nEstim=rep(2,3)), raster=r, ext=2, locErr=le<-.01, timeStep=1.6/60)
	  expect_equal(sum(values(p<-crop(l,e<-extent(.1,.1,0,0)+.6))) ,1/2, tolerance=2e-8)
	  a<-.2
	  s<-sqrt((1/60)*a*(1-a)*c(m@paraSd[1], m@orthSd[1])^2+a^2*le^2+(1-a)^2*le^2)
	  pp<-calc(rasterFromXYZ(rasterToPoints(p)[,c(1,2,1,2)]), function(x){prod(dnorm(x, c(.1,0), sd=s))})
	  expect_equal(values(p), values(pp)* prod(res(pp))*1/2, tolerance=7e-7)
	  expect_less_than(sum(values(abs(p-pp*prod(res(pp))/2))),7e-7)

	  expect_equal(sum(values(p<-crop(l,e<-extent(.9,.9,0,0)+.6))) ,1/2, tolerance=2e-8)
	  a<-.8
	  s<-sqrt((1/60)*a*(1-a)*c(m@paraSd[2], m@orthSd[2])^2+a^2*le^2+(1-a)^2*le^2)
	  pp<-calc(rasterFromXYZ(rasterToPoints(p)[,c(1,2,1,2)]), function(x){prod(dnorm(x, c(.9,0), sd=s))})
	  expect_equal(values(p), values(pp)* prod(res(pp))*1/2, tolerance=3e-6)
	  expect_less_than(sum(values(abs(p-pp*prod(res(pp))/2))),15e-7)
})
