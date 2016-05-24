context("addSize")

################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG
# 26.08.2014: Added test for scrampled markers (test 05 and test 06) [Issue#5].
# 07.05.2014: Added test for 'ESX17' (test 03 and test 04).
# 02.03.2014: First tests for 'addSize'.
# 
# test_dir("inst/tests/")
# test_file("tests/testthat/test-addSize.r")
# test_dir("tests/testthat")

test_that("addSize", {

  # Get test data.
  data(set2)

  # Correct marker order.
  scrampled <- rbind(set2[7:8,],set2[1:6,], set2[9:16,])
  
  # Get kit information for 'bins=TRUE'.
  kitBins <- getKit("SGMPlus", what="Size")
  
  # Get kit information for 'bins=FALSE'.
  kitCalc <- getKit("SGMPlus", what="Offset")
  
  # Get test data.
  data(ref4)
  # Extract one sample.
  ref4 <- ref4[ref4$Sample.Name=="A2",]
  
  # Get kit information for 'bins=TRUE'.
  kitBins2 <- getKit("ESX17", what="Size")
  
  # Get kit information for 'bins=FALSE'.
  kitCalc2 <- getKit("ESX17", what="Offset")
  
  
  # TEST 01 -------------------------------------------------------------------
  # Test adding size using bins=TRUE.
  
  # Analyse dataframe.
  res <- addSize(data=set2, kit=kitBins, bins=TRUE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))

  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res))
  expect_true("Marker" %in% names(res))
  expect_true("Allele" %in% names(res))
  expect_true("Height" %in% names(res))
  expect_true("Dye" %in% names(res))
  expect_true("Size" %in% names(res))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_true(any(is.na(res$Allele)))
  expect_true(any(is.na(res$Height)))
  expect_false(any(is.na(res$Dye)))
  expect_true(any(is.na(res$Size)))
  
  # Check result.
  expect_that(res$Size[1], equals(126))
  expect_that(res$Size[2], equals(138))
  expect_that(res$Size[3], equals(169))
  expect_that(res$Size[4], equals(258))
  expect_that(res$Size[5], equals(266))
  expect_that(res$Size[6], equals(305))
  expect_that(res$Size[7], equals(107))
  expect_that(res$Size[8], equals(113))
  expect_that(res$Size[9], equals(144))
  expect_that(res$Size[10], equals(211))
  expect_that(res$Size[11], equals(as.numeric(NA)))
  expect_that(res$Size[12], equals(301))
  expect_that(res$Size[13], equals(130))
  expect_that(res$Size[14], equals(173))
  expect_that(res$Size[15], equals(189))
  expect_that(res$Size[16], equals(247))
  expect_that(res$Size[17], equals(126))
  expect_that(res$Size[18], equals(138))
  expect_that(res$Size[19], equals(169))
  expect_that(res$Size[20], equals(258))
  expect_that(res$Size[21], equals(266))
  expect_that(res$Size[22], equals(305))
  expect_that(res$Size[23], equals(107))
  expect_that(res$Size[24], equals(113))
  expect_that(res$Size[25], equals(144))
  expect_that(res$Size[26], equals(211))
  expect_that(res$Size[27], equals(as.numeric(NA)))
  expect_that(res$Size[28], equals(301))
  expect_that(res$Size[29], equals(130))
  expect_that(res$Size[30], equals(173))
  expect_that(res$Size[31], equals(189))
  expect_that(res$Size[32], equals(as.numeric(NA)))
  
  # TEST 02 -------------------------------------------------------------------
  # Test adding size using bins=FALSE.
  
  # Analyse dataframe.
  res <- addSize(data=set2, kit=kitCalc, bins=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res))
  expect_true("Marker" %in% names(res))
  expect_true("Allele" %in% names(res))
  expect_true("Height" %in% names(res))
  expect_true("Dye" %in% names(res))
  expect_true("Size" %in% names(res))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_true(any(is.na(res$Allele)))
  expect_true(any(is.na(res$Height)))
  expect_false(any(is.na(res$Dye)))
  expect_true(any(is.na(res$Size)))
  
  # Check result.
  expect_that(res$Size[1], equals(126))
  expect_that(res$Size[2], equals(138))
  expect_that(res$Size[3], equals(169))
  expect_that(res$Size[4], equals(258))
  expect_that(res$Size[5], equals(266))
  expect_that(res$Size[6], equals(305))
  expect_that(res$Size[7], equals(107))
  expect_that(res$Size[8], equals(116))
  expect_that(res$Size[9], equals(144))
  expect_that(res$Size[10], equals(211))
  expect_that(res$Size[11], equals(363))
  expect_that(res$Size[12], equals(301))
  expect_that(res$Size[13], equals(130))
  expect_that(res$Size[14], equals(173))
  expect_that(res$Size[15], equals(189))
  expect_that(res$Size[16], equals(247))
  expect_that(res$Size[17], equals(126))
  expect_that(res$Size[18], equals(138))
  expect_that(res$Size[19], equals(169))
  expect_that(res$Size[20], equals(258))
  expect_that(res$Size[21], equals(266))
  expect_that(res$Size[22], equals(305))
  expect_that(res$Size[23], equals(107))
  expect_that(res$Size[24], equals(116))
  expect_that(res$Size[25], equals(144))
  expect_that(res$Size[26], equals(211))
  expect_that(res$Size[27], equals(363))
  expect_that(res$Size[28], equals(301))
  expect_that(res$Size[29], equals(130))
  expect_that(res$Size[30], equals(173))
  expect_that(res$Size[31], equals(189))
  expect_that(res$Size[32], equals(as.numeric(NA)))
  
  # TEST 03 -------------------------------------------------------------------
  # Test adding size using bins=TRUE.
  
  # Analyse dataframe.
  res <- addSize(data=ref4, kit=kitBins2, bins=TRUE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res))
  expect_true("Marker" %in% names(res))
  expect_true("Allele" %in% names(res))
  expect_true("Size" %in% names(res))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Allele)))
  expect_false(any(is.na(res$Size)))
  
  # Check result.
  expect_that(res$Size[1], equals(81.79))
  expect_that(res$Size[2], equals(81.79))
  expect_that(res$Size[3], equals(122.88))
  expect_that(res$Size[4], equals(126.97))
  expect_that(res$Size[5], equals(168.62))
  expect_that(res$Size[6], equals(168.62))
  expect_that(res$Size[7], equals(215.03))
  expect_that(res$Size[8], equals(219.03))
  expect_that(res$Size[9], equals(311.67))
  expect_that(res$Size[10], equals(315.56))
  expect_that(res$Size[11], equals(101.28))
  expect_that(res$Size[12], equals(105.42))
  expect_that(res$Size[13], equals(146.92))
  expect_that(res$Size[14], equals(159.18))
  expect_that(res$Size[15], equals(227.63))
  expect_that(res$Size[16], equals(247.60))
  expect_that(res$Size[17], equals(289.13))
  expect_that(res$Size[18], equals(301.08))
  expect_that(res$Size[19], equals(90.14))
  expect_that(res$Size[20], equals(105.15))
  expect_that(res$Size[21], equals(153.06))
  expect_that(res$Size[22], equals(161.12))
  expect_that(res$Size[23], equals(215.57))
  expect_that(res$Size[24], equals(227.59))
  expect_that(res$Size[25], equals(291.31))
  expect_that(res$Size[26], equals(306.69))
  expect_that(res$Size[27], equals(103.29))
  expect_that(res$Size[28], equals(112.67))
  expect_that(res$Size[29], equals(152.82))
  expect_that(res$Size[30], equals(169.06))
  expect_that(res$Size[31], equals(218.83))
  expect_that(res$Size[32], equals(226.77))
  expect_that(res$Size[33], equals(321.53))
  expect_that(res$Size[34], equals(368.23))
  
  # TEST 04 -------------------------------------------------------------------
  # Test adding size using bins=FALSE.
  
  # Analyse dataframe.
  res <- addSize(data=ref4, kit=kitCalc2, bins=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res))
  expect_true("Marker" %in% names(res))
  expect_true("Allele" %in% names(res))
  expect_true("Size" %in% names(res))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Allele)))
  expect_false(any(is.na(res$Size)))
  
  # Check result.
  expect_that(res$Size[1], equals(82))
  expect_that(res$Size[2], equals(82))
  expect_that(res$Size[3], equals(122))
  expect_that(res$Size[4], equals(126))
  expect_that(res$Size[5], equals(169))
  expect_that(res$Size[6], equals(169))
  expect_that(res$Size[7], equals(215))
  expect_that(res$Size[8], equals(219))
  expect_that(res$Size[9], equals(313))
  expect_that(res$Size[10], equals(317))
  expect_that(res$Size[11], equals(101))
  expect_that(res$Size[12], equals(105))
  expect_that(res$Size[13], equals(147))
  expect_that(res$Size[14], equals(159))
  expect_that(res$Size[15], equals(228))
  expect_that(res$Size[16], equals(248))
  expect_that(res$Size[17], equals(289))
  expect_that(res$Size[18], equals(301))
  expect_that(res$Size[19], equals(90))
  expect_that(res$Size[20], equals(105))
  expect_that(res$Size[21], equals(153))
  expect_that(res$Size[22], equals(161))
  expect_that(res$Size[23], equals(216))
  expect_that(res$Size[24], equals(228))
  expect_that(res$Size[25], equals(292))
  expect_that(res$Size[26], equals(308))
  expect_that(res$Size[27], equals(103))
  expect_that(res$Size[28], equals(112))
  expect_that(res$Size[29], equals(153))
  expect_that(res$Size[30], equals(169))
  expect_that(res$Size[31], equals(219))
  expect_that(res$Size[32], equals(227))
  expect_that(res$Size[33], equals(321))
  expect_that(res$Size[34], equals(367))
  
  # TEST 05 -------------------------------------------------------------------
  # Test adding size when marker order is not correct and bins=TRUE.

  # Analyse dataframe.
  res <- addSize(data=scrampled, kit=kitBins, bins=TRUE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res))
  expect_true("Marker" %in% names(res))
  expect_true("Allele" %in% names(res))
  expect_true("Height" %in% names(res))
  expect_true("Dye" %in% names(res))
  expect_true("Size" %in% names(res))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Allele)))
  expect_false(any(is.na(res$Height)))
  expect_false(any(is.na(res$Dye)))
  expect_true(any(is.na(res$Size)))
  
  # Check result.
  expect_that(res$Size[1], equals(107))
  expect_that(res$Size[2], equals(113))
  expect_that(res$Size[3], equals(126))
  expect_that(res$Size[4], equals(138))
  expect_that(res$Size[5], equals(169))
  expect_that(res$Size[6], equals(258))
  expect_that(res$Size[7], equals(266))
  expect_that(res$Size[8], equals(305))
  expect_that(res$Size[9], equals(144))
  expect_that(res$Size[10], equals(211))
  expect_that(res$Size[11], equals(as.numeric(NA)))
  expect_that(res$Size[12], equals(301))
  expect_that(res$Size[13], equals(130))
  expect_that(res$Size[14], equals(173))
  expect_that(res$Size[15], equals(189))
  expect_that(res$Size[16], equals(247))
  
  # TEST 06 -------------------------------------------------------------------
  # Test adding size when marker order is not correct and bins=FALSE.
  
  # Analyse dataframe.
  res <- addSize(data=scrampled, kit=kitCalc, bins=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res))
  expect_true("Marker" %in% names(res))
  expect_true("Allele" %in% names(res))
  expect_true("Height" %in% names(res))
  expect_true("Dye" %in% names(res))
  expect_true("Size" %in% names(res))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Allele)))
  expect_false(any(is.na(res$Height)))
  expect_false(any(is.na(res$Dye)))
  expect_false(any(is.na(res$Size)))
  
  # Check result.
  expect_that(res$Size[1], equals(107))
  expect_that(res$Size[2], equals(116))
  expect_that(res$Size[3], equals(126))
  expect_that(res$Size[4], equals(138))
  expect_that(res$Size[5], equals(169))
  expect_that(res$Size[6], equals(258))
  expect_that(res$Size[7], equals(266))
  expect_that(res$Size[8], equals(305))
  expect_that(res$Size[9], equals(144))
  expect_that(res$Size[10], equals(211))
  expect_that(res$Size[11], equals(363))
  expect_that(res$Size[12], equals(301))
  expect_that(res$Size[13], equals(130))
  expect_that(res$Size[14], equals(173))
  expect_that(res$Size[15], equals(189))
  expect_that(res$Size[16], equals(247))
  
  
  
})