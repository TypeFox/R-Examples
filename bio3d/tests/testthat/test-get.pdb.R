context("Testing get.pdb()")

test_that("get.pdb() works properly", {
  skip_on_cran()
  
   ids <- c("1tag", "1tnd")  # Gt
   tmp <- tempdir()
   files <- get.pdb(ids, tmp, verbose=FALSE)
   expect_identical(files, paste(tmp, "/", ids, ".pdb", sep=""))
   expect_warning(get.pdb("3c7kxxx", tmp, verbose=FALSE))
   expect_warning(get.pdb("1tag", tmp, verbose=FALSE))
   
   files <- get.pdb("1as0", tmp, verbose=FALSE, gzip=TRUE)
   expect_identical(files, paste(tmp, "/1as0.pdb", sep=""))
#   expect_error(get.pdb("aaaa", tmp, verbose=FALSE))
})

test_that("get.pdb() with ncore>1 works properly", {
  skip_on_cran()
     
   ids <- c("1tag", "1tnd", "3v00", "1got")
   tmp1 <- paste(tempdir(), "1", sep="")
   tmp2 <- paste(tempdir(), "2", sep="")
   time1 <- system.time(r1 <- get.pdb(ids, tmp1, ncore=1, verbose=FALSE))["elapsed"]
   time2 <- system.time(r2 <- get.pdb(ids, tmp2, ncore=NULL, verbose=FALSE))["elapsed"]
   expect_identical(r2, paste(tmp2, "/", ids, ".pdb", sep=""))
   expect_identical(list.files(tmp1), list.files(tmp2))
#   cat("Speed up by ", round((time1-time2)/time2*100, 1), "%", sep="")
#   if(getOption("cores") > 1)
#      expect_true(time1 > time2)
   unlink(tmp1, recursive=TRUE)
   unlink(tmp2, recursive=TRUE)
})
