context('smallFunctions')
test_that('extcalc',{
  expect_error(
    brownian.bridge.dyn(SpatialPointsDataFrame(cbind(1:3,1:3), data.frame(1:3)), dimSize=40,ext=1:3, location.error=3)
    ,'The ext argument must be a vector of 1, 2 or 4 numbers')

  #expect_identical(brownian.bridge.dyn(SpatialPointsDataFrame(cbind(1:3,1:3), data.frame(1:3)),
  #                             c(0,1)),c(1,3,-1,5))
}
)

test_that('trackId',
{
  m<-moveStack(list(
    move(1:4,1:4, Sys.time()+1:4, animal='a'),
    move(1:6,1:6, Sys.time()+1:6, animal='b')
    ))
  expect_identical(trackId(m) ,structure(c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L,2L,2L), .Label = c("a","b"), class = "factor"))
expect_equivalent(n.locs(m), as.array(c(a=4,b=6)))
}
  )
test_that('unUsed',{
data(leroy)
  expect_is(unUsedRecords(leroy),".unUsedRecords")
  expect_true(validObject(unUsedRecords(leroy)))
  expect_error(unUsedRecords(leroy)<-T,'Selection length does not match with number of locations')
})

test_that('linemidpoint',{
  a<-destPoint(4:5,5,1000)
  b<-destPoint(unlist(a), 123,2000)
  d<-destPoint(unlist(b), 13,1000)
  spdf<-SpatialPointsDataFrame(rbind(4:5,a,b,d),data = data.frame(a=4:7),proj4string = CRS('+proj=longlat +ellps=WGS84'),match.ID = F)
  
  u<-move:::lineMidpoint(spdf)
  uu<-move:::lineMidpoint(spdf[2:3,])
  expect_equal(u,uu)
  })

test_that('spatialLines',{
  data(leroy)
  crds<-coordinates(leroy)
 # rownames(crds)<-NULL # coordinates methods for lines removes rownames
    expect_equal(crds, coordinates(as(leroy,'SpatialLines'))[[1]][[1]])
    expect_equal(crds, coordinates(as(leroy,'SpatialLinesDataFrame'))[[1]][[1]])
  
data(fishers)


expect_equal(lapply(lapply(split(fishers), coordinates), function(x) {x}),lapply(coordinates(spldf<-as(fishers,'SpatialLinesDataFrame')),'[[',1))
  expect_equal(idData(fishers), data.frame(spldf))
  expect_equal(as(fishers,'SpatialLines'), as(spldf,'SpatialLines'))
  })

test_that('data.frame',{
data(leroy)  
    d<-slot(leroy,name = 'data')
    expect_equivalent(d, as(leroy,'data.frame')[,names(d)])
    data(fishers)
    dd<-slot(fishers,name = 'data')
    expect_equivalent(dd, as(fishers,'data.frame')[,names(dd)])
    
  })
