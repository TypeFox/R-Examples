context("Test for lets.presab.points")

# Species "occurrence" data
species <- c(rep("sp1", 100), rep("sp2", 100),
             rep("sp3", 100), rep("sp4", 100),
             rep("sp5", 1))
x <- runif(400, min = -69, max = -51)
y <- runif(400, min = -23, max = -4)
x <- c(x, 1000)
y <- c(y, 1000)
xy <- cbind(x, y)



test_that("lets.presab.points return a correct PresenceAbsence object", {
  

  PAM <- lets.presab.points(xy, species, xmn = -93,
                            xmx = -29, ymn = -57, 
                            ymx = 15)
  
  expect_equal(class(PAM), "PresenceAbsence")
  expect_equal(class(PAM[[1]]), "matrix")
  expect_true(inherits(PAM[[2]], "RasterLayer"))
  expect_equal(class(PAM[[3]]), "character")
  
  
})


test_that("lets.presab.points return a correct PresenceAbsence object for the world", {
  
  
  PAM <- lets.presab.points(xy, species)
  
  expect_equal(class(PAM), "PresenceAbsence")
  expect_equal(class(PAM[[1]]), "matrix")
  expect_true(inherits(PAM[[2]], "RasterLayer"))
  expect_equal(class(PAM[[3]]), "character")
  
  
})



test_that("lets.presab.points return a correct PresenceAbsence object (count=TRUE)", {
  
  
  PAM <- lets.presab.points(xy, species, xmn = -93,
                            xmx = -29, ymn = -57, 
                            ymx = 15, count = TRUE)
  
  expect_equal(class(PAM), "PresenceAbsence")
  expect_equal(class(PAM[[1]]), "matrix")
  expect_true(inherits(PAM[[2]], "RasterLayer"))
  expect_equal(class(PAM[[3]]), "character")
  
  
})




test_that("lets.presab.points return a correct PresenceAbsence object, remove.sp=FALSE", {
  
  
  PAM <- lets.presab.points(xy, species, xmn = -93,
                            xmx = -29, ymn = -57, 
                            ymx = 15, remove.sp = FALSE)
  
  expect_equal(class(PAM), "PresenceAbsence")
  expect_equal(class(PAM[[1]]), "matrix")
  expect_true(inherits(PAM[[2]], "RasterLayer"))
  expect_equal(class(PAM[[3]]), "character")
  
  response <- summary(PAM)
  expect_true(response$Specieswithoutanypresence > 0)
  
})


test_that("lets.presab.points return a correct PresenceAbsence object, remove.cells=FALSE", {
  
  
  PAM <- lets.presab.points(xy, species, xmn = -93,
                            xmx = -29, ymn = -57, 
                            ymx = 15, remove.cells = FALSE)
  
  expect_equal(class(PAM), "PresenceAbsence")
  expect_equal(class(PAM[[1]]), "matrix")
  expect_true(inherits(PAM[[2]], "RasterLayer"))
  expect_equal(class(PAM[[3]]), "character")
  
  
  response <- summary(PAM)
  expect_true(response$Cellswithoutanypresence > 0)
  
})
