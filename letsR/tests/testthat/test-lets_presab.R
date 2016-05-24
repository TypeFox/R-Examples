context("Test for lets.presab")

data(Phyllomedusa)

test_that("lets.presab return a correct PresenceAbsence object", {
  skip_on_cran()
  PAM <- lets.presab(Phyllomedusa, xmn=-93, xmx=-29, ymn= -57, ymx=15,
                     resol=1, remove.cells=TRUE, remove.sp=TRUE, show.matrix=FALSE,
                     crs=CRS("+proj=longlat +datum=WGS84"), cover=0, presence=NULL,
                     origin=NULL, seasonal=NULL, count=FALSE)
  
  expect_equal(class(PAM), "PresenceAbsence")
  expect_equal(class(PAM[[1]]), "matrix")
  expect_true(inherits(PAM[[2]], "RasterLayer"))
  expect_equal(class(PAM[[3]]), "character")
  
  
})



test_that("lets.presab return a correct PresenceAbsence object for the world", {
  skip_on_cran()
  
  PAM <- lets.presab(Phyllomedusa, resol=5, remove.cells=TRUE,
                     remove.sp=TRUE, show.matrix=FALSE,
                     crs=CRS("+proj=longlat +datum=WGS84"),
                     cover=0, presence=NULL, origin=NULL,
                     seasonal=NULL, count=FALSE)
  
  expect_equal(class(PAM), "PresenceAbsence")
  expect_equal(class(PAM[[1]]), "matrix")
  expect_true(inherits(PAM[[2]], "RasterLayer"))
  expect_equal(class(PAM[[3]]), "character")
  
  
})



test_that("lets.presab return a correct PresenceAbsence object (count=TRUE)", {
  skip_on_cran()
  
  PAM <- lets.presab(Phyllomedusa, xmn=-93, xmx=-29, ymn= -57, ymx=15,
                     resol=1, remove.cells=TRUE, remove.sp=TRUE, show.matrix=FALSE,
                     crs=CRS("+proj=longlat +datum=WGS84"), cover=0, presence=NULL,
                     origin=NULL, seasonal=NULL, count=TRUE)
  
  expect_equal(class(PAM), "PresenceAbsence")
  expect_equal(class(PAM[[1]]), "matrix")
  expect_true(inherits(PAM[[2]], "RasterLayer"))
  expect_equal(class(PAM[[3]]), "character")
  
  
})


test_that("lets.presab return a correct PresenceAbsence object, cover=0.2", {
  skip_on_cran()
  
  PAM <- lets.presab(Phyllomedusa, xmn=-93, xmx=-29, ymn= -57, ymx=15,
                     resol=1, remove.cells=TRUE, remove.sp=TRUE, show.matrix=FALSE,
                     crs=CRS("+proj=longlat +datum=WGS84"), cover=0.2, presence=NULL,
                     origin=NULL, seasonal=NULL, count=FALSE)
  
  expect_equal(class(PAM), "PresenceAbsence")
  expect_equal(class(PAM[[1]]), "matrix")
  expect_true(inherits(PAM[[2]], "RasterLayer"))
  expect_equal(class(PAM[[3]]), "character")
  
  
})


test_that("lets.presab return a correct PresenceAbsence object, remove.sp=FALSE", {
  skip_on_cran()
  
  PAM <- lets.presab(Phyllomedusa, xmn=-93, xmx=-29, ymn= -57, ymx=15,
                     resol=1, remove.cells=TRUE, remove.sp=FALSE, show.matrix=FALSE,
                     crs=CRS("+proj=longlat +datum=WGS84"), cover=0.2, presence=NULL,
                     origin=NULL, seasonal=NULL, count=FALSE)
  
  expect_equal(class(PAM), "PresenceAbsence")
  expect_equal(class(PAM[[1]]), "matrix")
  expect_true(inherits(PAM[[2]], "RasterLayer"))
  expect_equal(class(PAM[[3]]), "character")
  
  
  response <- summary(PAM)
  expect_true(response$Specieswithoutanypresence > 0)
  
})


test_that("lets.presab return a correct PresenceAbsence object, remove.cells=FALSE", {
  skip_on_cran()
  
  PAM <- lets.presab(Phyllomedusa, xmn=-93, xmx=-29, ymn= -57, ymx=15,
                     resol=1, remove.cells=FALSE, remove.sp=TRUE, show.matrix=FALSE,
                     crs=CRS("+proj=longlat +datum=WGS84"), cover=0, presence=NULL,
                     origin=NULL, seasonal=NULL, count=FALSE)
  
  expect_equal(class(PAM), "PresenceAbsence")
  expect_equal(class(PAM[[1]]), "matrix")
  expect_true(inherits(PAM[[2]], "RasterLayer"))
  expect_equal(class(PAM[[3]]), "character")
  
    
  response <- summary(PAM)
  expect_true(response$Cellswithoutanypresence > 0)
  
})



test_that("lets.presab new projection grid", {
  skip_on_cran()
  PAM2 <- lets.presab(Phyllomedusa, xmn = -3500000, xmx = 20000000,
                      ymn = -9000000, ymx = 0, 
                      crs.grid =  CRS("+init=epsg:26978"),
                      res = 100000, count = TRUE)
  
  expect_equal(class(PAM), "PresenceAbsence")
  expect_equal(class(PAM[[1]]), "matrix")
  expect_true(inherits(PAM[[2]], "RasterLayer"))
  expect_equal(class(PAM[[3]]), "character")  
  
})