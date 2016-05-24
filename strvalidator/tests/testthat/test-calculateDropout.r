context("calculateDropout")

################################################################################
# TODO LIST
# TODO: Update when dropout method is fixed.
# TODO: ...

################################################################################
# CHANGE LOG
# 07.12.2015: Fixed reference sample name subsetting bug.
# 
# require(testthat)
# test_dir("inst/tests/")
# test_file("tests/testthat/test-calculateDropout.r")
# test_dir("tests/testthat")

test_that("calculateDropout", {

  # Get test data.
  data(set4)
  data(ref4)
  
  # TEST 01 -------------------------------------------------------------------
  # Test that analysis of one sample works with double ref allele entries.

  # Get sample.
  testSample <- set4[set4$Sample.Name == "18-A2.14",]

  res <- calculateDropout(data=testSample, ref=ref4, ignore.case=TRUE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))

  # Check that expected columns and rows exist.  
  #expect_that(ncol(res), equals(7)) # TODO: Number of columns not fixed yet.
  expect_that(nrow(res), equals(24))
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  expect_true(any(grepl("Dropout", names(res))))
  expect_true(any(grepl("Rfu", names(res))))
  expect_true(any(grepl("Heterozygous", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Dropout)))
  expect_false(any(is.na(res$Heterozygous)))
  expect_true(any(is.na(res$Allele)))
  expect_true(any(is.na(res$Height)))
  expect_true(any(is.na(res$Rfu)))
  
  # Check result.
  expect_that(sum(res$Dropout==0), equals(16))
  expect_that(sum(res$Dropout==1), equals(7))
  expect_that(sum(res$Dropout==2), equals(1))

  expect_that(sum(res$Heterozygous==0), equals(2))
  expect_that(sum(res$Heterozygous==1), equals(22))
  
  expect_that(sum(res$Height, na.rm=TRUE), equals(5954))
  expect_that(sum(is.na(res$Height)), equals(1))
  
  expect_that(sum(res$Rfu, na.rm=TRUE), equals(1647))
  expect_that(sum(is.na(res$Rfu)), equals(17))
  
  # TEST 02 -------------------------------------------------------------------
  # Test that analysis of one sample works with single ref allele entries.
  
  # Get sample.
  testSample <- set4[set4$Sample.Name == "37-F2.8",]
  
  res <- calculateDropout(data=testSample, ref=ref4, ignore.case=TRUE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns and rows exist.  
  #expect_that(ncol(res), equals(7))  # TODO: Number of columns not fixed yet.
  expect_that(nrow(res), equals(18))
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  expect_true(any(grepl("Dropout", names(res))))
  expect_true(any(grepl("Rfu", names(res))))
  expect_true(any(grepl("Heterozygous", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Dropout)))
  expect_false(any(is.na(res$Heterozygous)))
  expect_true(any(is.na(res$Allele)))
  expect_true(any(is.na(res$Height)))
  expect_true(any(is.na(res$Rfu)))
  
  # Check result.
  expect_that(sum(res$Dropout==0), equals(5))
  expect_that(sum(res$Dropout==1), equals(4))
  expect_that(sum(res$Dropout==2), equals(9))
  
  expect_that(sum(res$Heterozygous==0), equals(4))
  expect_that(sum(res$Heterozygous==1), equals(14))
  
  expect_that(sum(res$Height, na.rm=TRUE), equals(2936))
  expect_that(sum(is.na(res$Height)), equals(9))
  
  expect_that(sum(res$Rfu, na.rm=TRUE), equals(802))
  expect_that(sum(is.na(res$Rfu)), equals(14))
  
  # TEST 03 -------------------------------------------------------------------
  # Test that analysis of one sample works with double ref allele entries
  # and ignore case = TRUE.
  
  # Get sample.
  testSample <- set4[set4$Sample.Name == "09-BC8",]
  
  res <- calculateDropout(data=testSample, ref=ref4, ignore.case=TRUE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns and rows exist.  
  #expect_that(ncol(res), equals(7))  # TODO: Number of columns not fixed yet.
  expect_that(nrow(res), equals(20))
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  expect_true(any(grepl("Dropout", names(res))))
  expect_true(any(grepl("Rfu", names(res))))
  expect_true(any(grepl("Heterozygous", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Dropout)))
  expect_false(any(is.na(res$Heterozygous)))
  expect_true(any(is.na(res$Allele)))
  expect_true(any(is.na(res$Height)))
  expect_true(any(is.na(res$Rfu)))
  
  # Check result.
  expect_that(sum(res$Dropout==0), equals(7))
  expect_that(sum(res$Dropout==1), equals(5))
  expect_that(sum(res$Dropout==2), equals(8))
  
  expect_that(sum(res$Heterozygous==0), equals(1))
  expect_that(sum(res$Heterozygous==1), equals(19))
  
  expect_that(sum(res$Height, na.rm=TRUE), equals(1532))
  expect_that(sum(is.na(res$Height)), equals(8))
  
  expect_that(sum(res$Rfu, na.rm=TRUE), equals(684))
  expect_that(sum(is.na(res$Rfu)), equals(15))
  
  
  # TEST 04 -------------------------------------------------------------------
  # Test that analysis of one sample works with double ref allele entries
  # and ignore case = FALSE.

  # Get sample.
  testSample <- set4[set4$Sample.Name == "09-BC8",]
  
  res <- calculateDropout(data=testSample, ref=ref4, ignore.case=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns and rows exist.  
  #expect_that(ncol(res), equals(7))  # TODO: Number of columns not fixed yet.
  expect_that(nrow(res), equals(0))
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  expect_true(any(grepl("Dropout", names(res))))
  expect_true(any(grepl("Rfu", names(res))))
  expect_true(any(grepl("Heterozygous", names(res))))
  
  # TEST 05 -------------------------------------------------------------------
  # Test that analysis of multiple samples works with double ref allele entries,
  # and ignore case = TRUE.

  # Get samples.
  testSample <- set4[set4$Sample.Name == "09-BC8" | set4$Sample.Name == "10-bc9",]
  
  res <- calculateDropout(data=testSample, ref=ref4, ignore.case=TRUE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns and rows exist.  
  #expect_that(ncol(res), equals(7))  # TODO: Number of columns not fixed yet.
  expect_that(nrow(res), equals(37))
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  expect_true(any(grepl("Dropout", names(res))))
  expect_true(any(grepl("Rfu", names(res))))
  expect_true(any(grepl("Heterozygous", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Dropout)))
  expect_false(any(is.na(res$Heterozygous)))
  expect_true(any(is.na(res$Allele)))
  expect_true(any(is.na(res$Height)))
  expect_true(any(is.na(res$Rfu)))
  
  # Check result.
  expect_that(sum(res$Dropout==0), equals(7))
  expect_that(sum(res$Dropout==1), equals(7))
  expect_that(sum(res$Dropout==2), equals(23))
  
  expect_that(sum(res$Heterozygous==0), equals(2))
  expect_that(sum(res$Heterozygous==1), equals(35))
  
  expect_that(sum(res$Height, na.rm=TRUE), equals(1746))
  expect_that(sum(is.na(res$Height)), equals(23))
  
  expect_that(sum(res$Rfu, na.rm=TRUE), equals(898))
  expect_that(sum(is.na(res$Rfu)), equals(30))

  
  # TEST 06 -------------------------------------------------------------------
  # Test that analysis of multiple samples work with double ref allele entries,
  # and ignore case = FALSE.

  # Get samples.
  testSample <- set4[set4$Sample.Name == "09-BC8" | set4$Sample.Name == "10-bc9",]
  
  res <- calculateDropout(data=testSample, ref=ref4, ignore.case=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns and rows exist.  
  #expect_that(ncol(res), equals(7))  # TODO: Number of columns not fixed yet.
  expect_that(nrow(res), equals(17))
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  expect_true(any(grepl("Dropout", names(res))))
  expect_true(any(grepl("Rfu", names(res))))
  expect_true(any(grepl("Heterozygous", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Dropout)))
  expect_false(any(is.na(res$Heterozygous)))
  expect_true(any(is.na(res$Allele)))
  expect_true(any(is.na(res$Height)))
  expect_true(any(is.na(res$Rfu)))
  
  # Check result.
  expect_that(sum(res$Dropout==0), equals(0))
  expect_that(sum(res$Dropout==1), equals(2))
  expect_that(sum(res$Dropout==2), equals(15))
  
  expect_that(sum(res$Heterozygous==0), equals(1))
  expect_that(sum(res$Heterozygous==1), equals(16))
  
  expect_that(sum(res$Height, na.rm=TRUE), equals(214))
  expect_that(sum(is.na(res$Height)), equals(15))
  
  expect_that(sum(res$Rfu, na.rm=TRUE), equals(214))
  expect_that(sum(is.na(res$Rfu)), equals(15))

  
  # TEST 07 -------------------------------------------------------------------
  # Test with missing markers.
  
  # Get samples.
  testSample <- set4[set4$Sample.Name == "09-BC8" | set4$Sample.Name == "10-bc9",]
  # Remove all rows without an allele.
  testSample <- testSample[!is.na(testSample$Allele),]
  
  res <- calculateDropout(data=testSample, ref=ref4, ignore.case=TRUE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns and rows exist.  
  #expect_that(ncol(res), equals(7))  # TODO: Number of columns not fixed yet.
  expect_that(nrow(res), equals(37))
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  expect_true(any(grepl("Dropout", names(res))))
  expect_true(any(grepl("Rfu", names(res))))
  expect_true(any(grepl("Heterozygous", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Dropout)))
  expect_false(any(is.na(res$Heterozygous)))
  expect_true(any(is.na(res$Allele)))
  expect_true(any(is.na(res$Height)))
  expect_true(any(is.na(res$Rfu)))
  
  # Check result.
  expect_that(sum(res$Dropout==0), equals(7))
  expect_that(sum(res$Dropout==1), equals(7))
  expect_that(sum(res$Dropout==2), equals(23))
  
  expect_that(sum(res$Heterozygous==0), equals(2))
  expect_that(sum(res$Heterozygous==1), equals(35))
  
  expect_that(sum(res$Height, na.rm=TRUE), equals(1746))
  expect_that(sum(is.na(res$Height)), equals(23))
  
  expect_that(sum(res$Rfu, na.rm=TRUE), equals(898))
  expect_that(sum(is.na(res$Rfu)), equals(30))

  # TEST 08 -------------------------------------------------------------------
  # Test with "OL" alleles and no dropout.
  
  # Get samples.
  testSample <- set4[set4$Sample.Name == "03-A2.1",]
  
  res <- calculateDropout(data=testSample, ref=ref4, ignore.case=TRUE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns and rows exist.  
  #expect_that(ncol(res), equals(7))  # TODO: Number of columns not fixed yet.
  expect_that(nrow(res), equals(32))
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  expect_true(any(grepl("Dropout", names(res))))
  expect_true(any(grepl("Rfu", names(res))))
  expect_true(any(grepl("Heterozygous", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Dropout)))
  expect_false(any(is.na(res$Heterozygous)))
  expect_false(any(is.na(res$Allele)))
  expect_false(any(is.na(res$Height)))
  expect_true(any(is.na(res$Rfu)))
  
  # Check result.
  expect_that(sum(res$Dropout==0), equals(32))
  expect_that(sum(res$Dropout==1), equals(0))
  expect_that(sum(res$Dropout==2), equals(0))
  
  expect_that(sum(res$Heterozygous==0), equals(2))
  expect_that(sum(res$Heterozygous==1), equals(30))
  
  expect_that(sum(res$Height, na.rm=TRUE), equals(583706))
  expect_that(sum(is.na(res$Height)), equals(0))
  
  expect_that(sum(res$Rfu, na.rm=TRUE), equals(0))
  expect_that(sum(is.na(res$Rfu)), equals(32))
  
  # TEST 09 -------------------------------------------------------------------
  # Test reference subsetting bug.
  # This test gives error in version 1.5.2 and earlier because
  # reference alleles from 'A2' and 'A2B' get mixed.

  # Get sample.
  testSample <- set4[set4$Sample.Name == "18-A2.14",]
  
  # Rename F2 to A2B wich match with A2
  testRef <- ref4
  testRef[testRef$Sample.Name == "F2", ]$Sample.Name <- "A2B"
  
  res <- calculateDropout(data=testSample, ref=testRef, ignore.case=TRUE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns and rows exist.  
  #expect_that(ncol(res), equals(7)) # TODO: Number of columns not fixed yet.
  expect_that(nrow(res), equals(24))
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  expect_true(any(grepl("Dropout", names(res))))
  expect_true(any(grepl("Rfu", names(res))))
  expect_true(any(grepl("Heterozygous", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Dropout)))
  expect_false(any(is.na(res$Heterozygous)))
  expect_true(any(is.na(res$Allele)))
  expect_true(any(is.na(res$Height)))
  expect_true(any(is.na(res$Rfu)))
  
  # Check result.
  expect_that(sum(res$Dropout==0), equals(16))
  expect_that(sum(res$Dropout==1), equals(7))
  expect_that(sum(res$Dropout==2), equals(1))
  
  expect_that(sum(res$Heterozygous==0), equals(2))
  expect_that(sum(res$Heterozygous==1), equals(22))
  
  expect_that(sum(res$Height, na.rm=TRUE), equals(5954))
  expect_that(sum(is.na(res$Height)), equals(1))
  
  expect_that(sum(res$Rfu, na.rm=TRUE), equals(1647))
  expect_that(sum(is.na(res$Rfu)), equals(17))
  
  
})