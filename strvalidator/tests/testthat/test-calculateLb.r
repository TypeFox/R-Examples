context("calculateLb")

################################################################################
# TODO LIST
# TODO: Test ignore.case.
# TODO: ...

################################################################################
# CHANGE LOG
# 26.12.2015: First version.
# 
# require(strvalidator)
# require(testthat)
# test_dir("inst/tests/")
# test_file("tests/testthat/test-calculateLb.r")
# test_dir("tests/testthat")

test_that("calculateLb", {

  # Get test data.
  data(set1)
  data(ref1)
  data(set2)
  data(ref2)
  
  # Extract sample PC2 and fix datasets.
  set1Fat <- set1[set1$Sample.Name == "PC2", ]
  set1 <- slim(data = set1Fat, fix = c("Sample.Name", "Marker", "Dye"), stack = c("Allele", "Height"))
  ref1 <- slim(data = ref1, fix = c("Sample.Name", "Marker"), stack = c("Allele"))
  
  # TEST 01 -------------------------------------------------------------------
  # Test that proportional Lb by dye can be calculated by dye for data including
  # one locus dropout.
  
  # Analyse dataframe.
  res <- calculateLb(data = set2, ref = NULL, option = "prop", by.dye = TRUE,
                     ol.rm = FALSE, sex.rm = FALSE, na = NULL,
                     ignore.case = TRUE, word = FALSE, exact = FALSE)

  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_false(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$TPH))
  expect_false(is.null(res$Peaks))
  expect_false(is.null(res$Dye))
  expect_false(is.null(res$TPPH))
  expect_false(is.null(res$Lb))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_true(any(is.na(res$TPH)))
  expect_false(any(is.na(res$Peaks)))
  expect_false(any(is.na(res$Dye)))
  expect_true(any(is.na(res$TPPH)))
  expect_true(any(is.na(res$Lb)))
  
  # Check result: Locus balance.
  expect_that(res$Lb[1], equals(862/3005))
  expect_that(res$Lb[2], equals(506/3005))
  expect_that(res$Lb[3], equals(914/3005))
  expect_that(res$Lb[4], equals(723/3005))
  expect_that(res$Lb[5], equals(1219/3363))
  expect_that(res$Lb[6], equals(619/3363))
  expect_that(res$Lb[7], equals(766/3363))
  expect_that(res$Lb[8], equals(759/3363))
  expect_that(res$Lb[9], equals(592/1760))
  expect_that(res$Lb[10], equals(743/1760))
  expect_that(res$Lb[11], equals(425/1760))
  
  expect_that(res$Lb[12], equals(440/1530))
  expect_that(res$Lb[13], equals(304/1530))
  expect_that(res$Lb[14], equals(439/1530))
  expect_that(res$Lb[15], equals(347/1530))
  expect_that(res$Lb[16], equals(638/1750))
  expect_that(res$Lb[17], equals(309/1750))
  expect_that(res$Lb[18], equals(402/1750))
  expect_that(res$Lb[19], equals(401/1750))
  expect_that(res$Lb[20], equals(as.numeric(NA)))
  expect_that(res$Lb[21], equals(as.numeric(NA)))
  expect_that(res$Lb[22], equals(as.numeric(NA)))
  
  # Check result: Total peak height.
  expect_that(res$TPH[1], equals(862))
  expect_that(res$TPH[2], equals(506))
  expect_that(res$TPH[3], equals(914))
  expect_that(res$TPH[4], equals(723))
  expect_that(res$TPH[5], equals(1219))
  expect_that(res$TPH[6], equals(619))
  expect_that(res$TPH[7], equals(766))
  expect_that(res$TPH[8], equals(759))
  expect_that(res$TPH[9], equals(592))
  expect_that(res$TPH[10], equals(743))
  expect_that(res$TPH[11], equals(425))
  
  expect_that(res$TPH[12], equals(440))
  expect_that(res$TPH[13], equals(304))
  expect_that(res$TPH[14], equals(439))
  expect_that(res$TPH[15], equals(347))
  expect_that(res$TPH[16], equals(638))
  expect_that(res$TPH[17], equals(309))
  expect_that(res$TPH[18], equals(402))
  expect_that(res$TPH[19], equals(401))
  expect_that(res$TPH[20], equals(284))
  expect_that(res$TPH[21], equals(362))
  expect_that(res$TPH[22], equals(as.numeric(NA)))

  # Check result: Total profile peak height.
  expect_that(unique(res$TPPH), equals(c(3005,3363,1760,1530,1750,NA)))
  

  # TEST 02 -------------------------------------------------------------------
  # Test that global proportional Lb can be calculated for data including
  # one locus dropout.
  
  res <- calculateLb(data = set2, ref = NULL, option = "prop", by.dye = FALSE,
                     ol.rm = FALSE, sex.rm = FALSE, na = NULL,
                     ignore.case = TRUE, word = FALSE, exact = FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_false(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$TPH))
  expect_false(is.null(res$Peaks))
  expect_false(is.null(res$Dye))
  expect_false(is.null(res$TPPH))
  expect_false(is.null(res$Lb))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_true(any(is.na(res$TPH)))
  expect_false(any(is.na(res$Peaks)))
  expect_false(any(is.na(res$Dye)))
  expect_true(any(is.na(res$TPPH)))
  expect_true(any(is.na(res$Lb)))
  
  # Check result: Locus balance.
  expect_that(res$Lb[1], equals(862/8128))
  expect_that(res$Lb[2], equals(506/8128))
  expect_that(res$Lb[3], equals(914/8128))
  expect_that(res$Lb[4], equals(723/8128))
  expect_that(res$Lb[5], equals(1219/8128))
  expect_that(res$Lb[6], equals(619/8128))
  expect_that(res$Lb[7], equals(766/8128))
  expect_that(res$Lb[8], equals(759/8128))
  expect_that(res$Lb[9], equals(592/8128))
  expect_that(res$Lb[10], equals(743/8128))
  expect_that(res$Lb[11], equals(425/8128))
  
  expect_that(res$Lb[12], equals(as.numeric(NA)))
  expect_that(res$Lb[13], equals(as.numeric(NA)))
  expect_that(res$Lb[14], equals(as.numeric(NA)))
  expect_that(res$Lb[15], equals(as.numeric(NA)))
  expect_that(res$Lb[16], equals(as.numeric(NA)))
  expect_that(res$Lb[17], equals(as.numeric(NA)))
  expect_that(res$Lb[18], equals(as.numeric(NA)))
  expect_that(res$Lb[19], equals(as.numeric(NA)))
  expect_that(res$Lb[20], equals(as.numeric(NA)))
  expect_that(res$Lb[21], equals(as.numeric(NA)))
  expect_that(res$Lb[22], equals(as.numeric(NA)))

  # Check result: Total peak height.
  expect_that(res$TPH[1], equals(862))
  expect_that(res$TPH[2], equals(506))
  expect_that(res$TPH[3], equals(914))
  expect_that(res$TPH[4], equals(723))
  expect_that(res$TPH[5], equals(1219))
  expect_that(res$TPH[6], equals(619))
  expect_that(res$TPH[7], equals(766))
  expect_that(res$TPH[8], equals(759))
  expect_that(res$TPH[9], equals(592))
  expect_that(res$TPH[10], equals(743))
  expect_that(res$TPH[11], equals(425))
  
  expect_that(res$TPH[12], equals(440))
  expect_that(res$TPH[13], equals(304))
  expect_that(res$TPH[14], equals(439))
  expect_that(res$TPH[15], equals(347))
  expect_that(res$TPH[16], equals(638))
  expect_that(res$TPH[17], equals(309))
  expect_that(res$TPH[18], equals(402))
  expect_that(res$TPH[19], equals(401))
  expect_that(res$TPH[20], equals(284))
  expect_that(res$TPH[21], equals(362))
  expect_that(res$TPH[22], equals(as.numeric(NA)))

  # Check result: Total profile peak height.
  expect_that(unique(res$TPPH), equals(c(8128,NA)))

    
  # TEST 03 -------------------------------------------------------------------
  # Test that normalized Lb by dye can be calculated for data including
  # one locus dropout.
  
  res <- calculateLb(data = set2, ref = NULL, option = "norm", by.dye = TRUE,
                     ol.rm = FALSE, sex.rm = FALSE, na=NULL,
                     ignore.case = TRUE, word = FALSE, exact = FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_false(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$TPH))
  expect_false(is.null(res$Peaks))
  expect_false(is.null(res$Dye))
  expect_false(is.null(res$MTPH))
  expect_false(is.null(res$Lb))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_true(any(is.na(res$TPH)))
  expect_false(any(is.na(res$Peaks)))
  expect_false(any(is.na(res$Dye)))
  expect_true(any(is.na(res$MTPH)))
  expect_true(any(is.na(res$Lb)))

  # Check result: Locus balance.
  expect_that(res$Lb[1], equals(862/914))
  expect_that(res$Lb[2], equals(506/914))
  expect_that(res$Lb[3], equals(914/914))
  expect_that(res$Lb[4], equals(723/914))
  expect_that(res$Lb[5], equals(1219/1219))
  expect_that(res$Lb[6], equals(619/1219))
  expect_that(res$Lb[7], equals(766/1219))
  expect_that(res$Lb[8], equals(759/1219))
  expect_that(res$Lb[9], equals(592/743))
  expect_that(res$Lb[10], equals(743/743))
  expect_that(res$Lb[11], equals(425/743))
  
  expect_that(res$Lb[12], equals(440/440))
  expect_that(res$Lb[13], equals(304/440))
  expect_that(res$Lb[14], equals(439/440))
  expect_that(res$Lb[15], equals(347/440))
  expect_that(res$Lb[16], equals(638/638))
  expect_that(res$Lb[17], equals(309/638))
  expect_that(res$Lb[18], equals(402/638))
  expect_that(res$Lb[19], equals(401/638))
  expect_that(res$Lb[20], equals(as.numeric(NA)))
  expect_that(res$Lb[21], equals(as.numeric(NA)))
  expect_that(res$Lb[22], equals(as.numeric(NA)))

  # Check result: Total peak height.
  expect_that(res$TPH[1], equals(862))
  expect_that(res$TPH[2], equals(506))
  expect_that(res$TPH[3], equals(914))
  expect_that(res$TPH[4], equals(723))
  expect_that(res$TPH[5], equals(1219))
  expect_that(res$TPH[6], equals(619))
  expect_that(res$TPH[7], equals(766))
  expect_that(res$TPH[8], equals(759))
  expect_that(res$TPH[9], equals(592))
  expect_that(res$TPH[10], equals(743))
  expect_that(res$TPH[11], equals(425))
  
  expect_that(res$TPH[12], equals(440))
  expect_that(res$TPH[13], equals(304))
  expect_that(res$TPH[14], equals(439))
  expect_that(res$TPH[15], equals(347))
  expect_that(res$TPH[16], equals(638))
  expect_that(res$TPH[17], equals(309))
  expect_that(res$TPH[18], equals(402))
  expect_that(res$TPH[19], equals(401))
  expect_that(res$TPH[20], equals(284))
  expect_that(res$TPH[21], equals(362))
  expect_that(res$TPH[22], equals(as.numeric(NA)))
  
  # Check result: Maximum total peak height.
  expect_that(unique(res$MTPH), equals(c(914,1219,743,440,638,NA)))

    
  # TEST 04 -------------------------------------------------------------------
  # Test that global normalized Lb can be calculated for data including
  # one locus dropout.
  
  res <- calculateLb(data = set2, ref = NULL, option = "norm", by.dye = FALSE,
                     ol.rm = FALSE, sex.rm = FALSE, na=NULL,
                     ignore.case = TRUE, word = FALSE, exact = FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_false(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$TPH))
  expect_false(is.null(res$Peaks))
  expect_false(is.null(res$Dye))
  expect_false(is.null(res$MTPH))
  expect_false(is.null(res$Lb))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_true(any(is.na(res$TPH)))
  expect_false(any(is.na(res$Peaks)))
  expect_false(any(is.na(res$Dye)))
  expect_true(any(is.na(res$MTPH)))
  expect_true(any(is.na(res$Lb)))
  
  # Check result: Locus balance.
  expect_that(res$Lb[1], equals(862/1219))
  expect_that(res$Lb[2], equals(506/1219))
  expect_that(res$Lb[3], equals(914/1219))
  expect_that(res$Lb[4], equals(723/1219))
  expect_that(res$Lb[5], equals(1219/1219))
  expect_that(res$Lb[6], equals(619/1219))
  expect_that(res$Lb[7], equals(766/1219))
  expect_that(res$Lb[8], equals(759/1219))
  expect_that(res$Lb[9], equals(592/1219))
  expect_that(res$Lb[10], equals(743/1219))
  expect_that(res$Lb[11], equals(425/1219))
  
  expect_that(res$Lb[12], equals(as.numeric(NA)))
  expect_that(res$Lb[13], equals(as.numeric(NA)))
  expect_that(res$Lb[14], equals(as.numeric(NA)))
  expect_that(res$Lb[15], equals(as.numeric(NA)))
  expect_that(res$Lb[16], equals(as.numeric(NA)))
  expect_that(res$Lb[17], equals(as.numeric(NA)))
  expect_that(res$Lb[18], equals(as.numeric(NA)))
  expect_that(res$Lb[19], equals(as.numeric(NA)))
  expect_that(res$Lb[20], equals(as.numeric(NA)))
  expect_that(res$Lb[21], equals(as.numeric(NA)))
  expect_that(res$Lb[22], equals(as.numeric(NA)))
  
  # Check result: Total peak height.
  expect_that(res$TPH[1], equals(862))
  expect_that(res$TPH[2], equals(506))
  expect_that(res$TPH[3], equals(914))
  expect_that(res$TPH[4], equals(723))
  expect_that(res$TPH[5], equals(1219))
  expect_that(res$TPH[6], equals(619))
  expect_that(res$TPH[7], equals(766))
  expect_that(res$TPH[8], equals(759))
  expect_that(res$TPH[9], equals(592))
  expect_that(res$TPH[10], equals(743))
  expect_that(res$TPH[11], equals(425))
  
  expect_that(res$TPH[12], equals(440))
  expect_that(res$TPH[13], equals(304))
  expect_that(res$TPH[14], equals(439))
  expect_that(res$TPH[15], equals(347))
  expect_that(res$TPH[16], equals(638))
  expect_that(res$TPH[17], equals(309))
  expect_that(res$TPH[18], equals(402))
  expect_that(res$TPH[19], equals(401))
  expect_that(res$TPH[20], equals(284))
  expect_that(res$TPH[21], equals(362))
  expect_that(res$TPH[22], equals(as.numeric(NA)))

  # Check result: Maximum total peak height.
  expect_that(unique(res$MTPH), equals(c(1219,NA)))
  
  
  # TEST 05 -------------------------------------------------------------------
  # Test that centred quantities Lb by dye can be calculated for data including
  # one locus dropout.
  
  res <- calculateLb(data = set2, ref = NULL, option = "cent", by.dye = TRUE,
                     ol.rm = FALSE, sex.rm = FALSE, na=NULL,
                     ignore.case = TRUE, word = FALSE, exact = FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_false(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$TPH))
  expect_false(is.null(res$Peaks))
  expect_false(is.null(res$Dye))
  expect_false(is.null(res$MPH))
  expect_false(is.null(res$Lb))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_true(any(is.na(res$TPH)))
  expect_false(any(is.na(res$Peaks)))
  expect_false(any(is.na(res$Dye)))
  expect_true(any(is.na(res$MPH)))
  expect_true(any(is.na(res$Lb)))
  
  # Check result: Locus balance.
  expect_that(res$Lb[1], equals((862-751.25)/sqrt(751.25)))
  expect_that(res$Lb[2], equals((506-751.25)/sqrt(751.25)))
  expect_that(res$Lb[3], equals((914-751.25)/sqrt(751.25)))
  expect_that(res$Lb[4], equals((723-751.25)/sqrt(751.25)))
  expect_that(res$Lb[5], equals((1219-840.75)/sqrt(840.75)))
  expect_that(res$Lb[6], equals((619-840.75)/sqrt(840.75)))
  expect_that(res$Lb[7], equals((766-840.75)/sqrt(840.75)))
  expect_that(res$Lb[8], equals((759-840.75)/sqrt(840.75)))
  expect_that(res$Lb[9], equals((592-(586+2/3))/sqrt(586+2/3)))
  expect_that(res$Lb[10], equals((743-(586+2/3))/sqrt(586+2/3)))
  expect_that(res$Lb[11], equals((425-(586+2/3))/sqrt(586+2/3)))
  
  expect_that(res$Lb[12], equals((440-382.5)/sqrt(382.5)))
  expect_that(res$Lb[13], equals((304-382.5)/sqrt(382.5)))
  expect_that(res$Lb[14], equals((439-382.5)/sqrt(382.5)))
  expect_that(res$Lb[15], equals((347-382.5)/sqrt(382.5)))
  expect_that(res$Lb[16], equals((638-437.5)/sqrt(437.5)))
  expect_that(res$Lb[17], equals((309-437.5)/sqrt(437.5)))
  expect_that(res$Lb[18], equals((402-437.5)/sqrt(437.5)))
  expect_that(res$Lb[19], equals((401-437.5)/sqrt(437.5)))
  expect_that(res$Lb[20], equals(as.numeric(NA)))
  expect_that(res$Lb[21], equals(as.numeric(NA)))
  expect_that(res$Lb[22], equals(as.numeric(NA)))
  
  # Check result: Total peak height.
  expect_that(res$TPH[1], equals(862))
  expect_that(res$TPH[2], equals(506))
  expect_that(res$TPH[3], equals(914))
  expect_that(res$TPH[4], equals(723))
  expect_that(res$TPH[5], equals(1219))
  expect_that(res$TPH[6], equals(619))
  expect_that(res$TPH[7], equals(766))
  expect_that(res$TPH[8], equals(759))
  expect_that(res$TPH[9], equals(592))
  expect_that(res$TPH[10], equals(743))
  expect_that(res$TPH[11], equals(425))
  
  expect_that(res$TPH[12], equals(440))
  expect_that(res$TPH[13], equals(304))
  expect_that(res$TPH[14], equals(439))
  expect_that(res$TPH[15], equals(347))
  expect_that(res$TPH[16], equals(638))
  expect_that(res$TPH[17], equals(309))
  expect_that(res$TPH[18], equals(402))
  expect_that(res$TPH[19], equals(401))
  expect_that(res$TPH[20], equals(284))
  expect_that(res$TPH[21], equals(362))
  expect_that(res$TPH[22], equals(as.numeric(NA)))

  # Check result: Mean peak height.
  expect_that(unique(res$MPH), equals(c(751.25,840.75,586+2/3,382.5,437.5,NA)))

  
  # TEST 06 -------------------------------------------------------------------
  # Test that global centred quantity Lb can be calculated for data including
  # one locus dropout.
  
  res <- calculateLb(data = set2, ref = NULL, option = "cent", by.dye = FALSE,
                     ol.rm = FALSE, sex.rm = FALSE, na=NULL,
                     ignore.case = TRUE, word = FALSE, exact = FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_false(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$TPH))
  expect_false(is.null(res$Peaks))
  expect_false(is.null(res$Dye))
  expect_false(is.null(res$MPH))
  expect_false(is.null(res$Lb))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_true(any(is.na(res$TPH)))
  expect_false(any(is.na(res$Peaks)))
  expect_false(any(is.na(res$Dye)))
  expect_true(any(is.na(res$MPH)))
  expect_true(any(is.na(res$Lb)))
  
  # Check result: Locus balance.
  expect_that(round(res$Lb[1],4), equals(round((862-738.9091)/sqrt(738.9091),4)))
  expect_that(round(res$Lb[2],4), equals(round((506-738.9091)/sqrt(738.9091),4)))
  expect_that(round(res$Lb[3],4), equals(round((914-738.9091)/sqrt(738.9091),4)))
  expect_that(round(res$Lb[4],4), equals(round((723-738.9091)/sqrt(738.9091),4)))
  expect_that(round(res$Lb[5],4), equals(round((1219-738.9091)/sqrt(738.9091),4)))
  expect_that(round(res$Lb[6],4), equals(round((619-738.9091)/sqrt(738.9091),4)))
  expect_that(round(res$Lb[7],4), equals(round((766-738.9091)/sqrt(738.9091),4)))
  expect_that(round(res$Lb[8],4), equals(round((759-738.9091)/sqrt(738.9091),4)))
  expect_that(round(res$Lb[9],4), equals(round((592-738.9091)/sqrt(738.9091),4)))
  expect_that(round(res$Lb[10],4), equals(round((743-738.9091)/sqrt(738.9091),4)))
  expect_that(round(res$Lb[11],4), equals(round((425-738.9091)/sqrt(738.9091),4)))

  expect_that(res$Lb[12], equals(as.numeric(NA)))
  expect_that(res$Lb[13], equals(as.numeric(NA)))
  expect_that(res$Lb[14], equals(as.numeric(NA)))
  expect_that(res$Lb[15], equals(as.numeric(NA)))
  expect_that(res$Lb[16], equals(as.numeric(NA)))
  expect_that(res$Lb[17], equals(as.numeric(NA)))
  expect_that(res$Lb[18], equals(as.numeric(NA)))
  expect_that(res$Lb[19], equals(as.numeric(NA)))
  expect_that(res$Lb[20], equals(as.numeric(NA)))
  expect_that(res$Lb[21], equals(as.numeric(NA)))
  expect_that(res$Lb[22], equals(as.numeric(NA)))
  
  # Check result: Total peak height.
  expect_that(res$TPH[1], equals(862))
  expect_that(res$TPH[2], equals(506))
  expect_that(res$TPH[3], equals(914))
  expect_that(res$TPH[4], equals(723))
  expect_that(res$TPH[5], equals(1219))
  expect_that(res$TPH[6], equals(619))
  expect_that(res$TPH[7], equals(766))
  expect_that(res$TPH[8], equals(759))
  expect_that(res$TPH[9], equals(592))
  expect_that(res$TPH[10], equals(743))
  expect_that(res$TPH[11], equals(425))
  
  expect_that(res$TPH[12], equals(440))
  expect_that(res$TPH[13], equals(304))
  expect_that(res$TPH[14], equals(439))
  expect_that(res$TPH[15], equals(347))
  expect_that(res$TPH[16], equals(638))
  expect_that(res$TPH[17], equals(309))
  expect_that(res$TPH[18], equals(402))
  expect_that(res$TPH[19], equals(401))
  expect_that(res$TPH[20], equals(284))
  expect_that(res$TPH[21], equals(362))
  expect_that(res$TPH[22], equals(as.numeric(NA)))
  
  # Check result: Mean peak height.
  expect_that(unique(res$MPH), equals(c(738.9091,NA)))
  
  
  # TEST 07 -------------------------------------------------------------------
  # Test that proportional Lb by dye can be calculated for unfiltered data.

  # Analyse dataframe.
  res <- calculateLb(data = set1, ref = NULL, option = "prop", by.dye = TRUE,
                     ol.rm = FALSE, sex.rm = FALSE, na = NULL,
                     ignore.case = TRUE, word = FALSE, exact = FALSE)

  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_false(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$TPH))
  expect_false(is.null(res$Peaks))
  expect_false(is.null(res$Dye))
  expect_false(is.null(res$TPPH))
  expect_false(is.null(res$Lb))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$TPH)))
  expect_false(any(is.na(res$Peaks)))
  expect_false(any(is.na(res$Dye)))
  expect_false(any(is.na(res$TPPH)))
  expect_false(any(is.na(res$Lb)))
  
  # Check result: Locus balance.
  expect_that(res$Lb[1], equals(7324/34486))
  expect_that(res$Lb[2], equals(7648/34486))
  expect_that(res$Lb[3], equals(7911/34486))
  expect_that(res$Lb[4], equals(6774/34486))
  expect_that(res$Lb[5], equals(4829/34486))
  expect_that(res$Lb[6], equals(5482/26214))
  expect_that(res$Lb[7], equals(8098/26214))
  expect_that(res$Lb[8], equals(6629/26214))
  expect_that(res$Lb[9], equals(6005/26214))
  expect_that(res$Lb[10], equals(5548/24822))
  expect_that(res$Lb[11], equals(4975/24822))
  expect_that(res$Lb[12], equals(7320/24822))
  expect_that(res$Lb[13], equals(6979/24822))
  expect_that(res$Lb[14], equals(16865/50136))
  expect_that(res$Lb[15], equals(10906/50136))
  expect_that(res$Lb[16], equals(10649/50136))
  expect_that(res$Lb[17], equals(11716/50136))
  
  # Check result: Total peak height.
  expect_that(res$TPH[1], equals(7324))
  expect_that(res$TPH[2], equals(7648))
  expect_that(res$TPH[3], equals(7911))
  expect_that(res$TPH[4], equals(6774))
  expect_that(res$TPH[5], equals(4829))
  expect_that(res$TPH[6], equals(5482))
  expect_that(res$TPH[7], equals(8098))
  expect_that(res$TPH[8], equals(6629))
  expect_that(res$TPH[9], equals(6005))
  expect_that(res$TPH[10], equals(5548))
  expect_that(res$TPH[11], equals(4975))
  expect_that(res$TPH[12], equals(7320))
  expect_that(res$TPH[13], equals(6979))
  expect_that(res$TPH[14], equals(16865))
  expect_that(res$TPH[15], equals(10906))
  expect_that(res$TPH[16], equals(10649))
  expect_that(res$TPH[17], equals(11716))
  
  # Check result: Total profile peak height.
  expect_that(unique(res$TPPH), equals(c(34486,26214,24822,50136)))
  
  
  # TEST 08 -------------------------------------------------------------------
  # Test that proportional Lb by dye can be calculated for unfiltered data
  # and off-ladder peaks removed.
  
  # Analyse dataframe.
  res <- calculateLb(data = set1, ref = NULL, option = "prop", by.dye = TRUE,
                     ol.rm = TRUE, sex.rm = FALSE, na = NULL,
                     ignore.case = TRUE, word = FALSE, exact = FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_false(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$TPH))
  expect_false(is.null(res$Peaks))
  expect_false(is.null(res$Dye))
  expect_false(is.null(res$TPPH))
  expect_false(is.null(res$Lb))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$TPH)))
  expect_false(any(is.na(res$Peaks)))
  expect_false(any(is.na(res$Dye)))
  expect_false(any(is.na(res$TPPH)))
  expect_false(any(is.na(res$Lb)))
  
  # Check result: Locus balance.
  expect_that(res$Lb[1], equals(7131/34293))
  expect_that(res$Lb[2], equals(7648/34293))
  expect_that(res$Lb[3], equals(7911/34293))
  expect_that(res$Lb[4], equals(6774/34293))
  expect_that(res$Lb[5], equals(4829/34293))
  expect_that(res$Lb[6], equals(5482/26214))
  expect_that(res$Lb[7], equals(8098/26214))
  expect_that(res$Lb[8], equals(6629/26214))
  expect_that(res$Lb[9], equals(6005/26214))
  expect_that(res$Lb[10], equals(5548/24822))
  expect_that(res$Lb[11], equals(4975/24822))
  expect_that(res$Lb[12], equals(7320/24822))
  expect_that(res$Lb[13], equals(6979/24822))
  expect_that(res$Lb[14], equals(16865/50136))
  expect_that(res$Lb[15], equals(10906/50136))
  expect_that(res$Lb[16], equals(10649/50136))
  expect_that(res$Lb[17], equals(11716/50136))
  
  # Check result: Total peak height.
  expect_that(res$TPH[1], equals(7131))
  expect_that(res$TPH[2], equals(7648))
  expect_that(res$TPH[3], equals(7911))
  expect_that(res$TPH[4], equals(6774))
  expect_that(res$TPH[5], equals(4829))
  expect_that(res$TPH[6], equals(5482))
  expect_that(res$TPH[7], equals(8098))
  expect_that(res$TPH[8], equals(6629))
  expect_that(res$TPH[9], equals(6005))
  expect_that(res$TPH[10], equals(5548))
  expect_that(res$TPH[11], equals(4975))
  expect_that(res$TPH[12], equals(7320))
  expect_that(res$TPH[13], equals(6979))
  expect_that(res$TPH[14], equals(16865))
  expect_that(res$TPH[15], equals(10906))
  expect_that(res$TPH[16], equals(10649))
  expect_that(res$TPH[17], equals(11716))
  
  # Check result: Total profile peak height.
  expect_that(unique(res$TPPH), equals(c(34293,26214,24822,50136)))
  
  
  # TEST 09 -------------------------------------------------------------------
  # Test that proportional Lb by dye can be calculated for unfiltered data,
  # with off-ladder peaks, and sex markers removed.
  
  # Analyse dataframe.
  res <- calculateLb(data = set1, ref = NULL, option = "prop", by.dye = TRUE,
                     ol.rm = TRUE, sex.rm = TRUE, na = NULL,
                     ignore.case = TRUE, word = FALSE, exact = FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_false(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$TPH))
  expect_false(is.null(res$Peaks))
  expect_false(is.null(res$Dye))
  expect_false(is.null(res$TPPH))
  expect_false(is.null(res$Lb))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$TPH)))
  expect_false(any(is.na(res$Peaks)))
  expect_false(any(is.na(res$Dye)))
  expect_false(any(is.na(res$TPPH)))
  expect_false(any(is.na(res$Lb)))
  
  # Check result: Locus balance.
  expect_that(res$Lb[1], equals(7648/27162))
  expect_that(res$Lb[2], equals(7911/27162))
  expect_that(res$Lb[3], equals(6774/27162))
  expect_that(res$Lb[4], equals(4829/27162))
  expect_that(res$Lb[5], equals(5482/26214))
  expect_that(res$Lb[6], equals(8098/26214))
  expect_that(res$Lb[7], equals(6629/26214))
  expect_that(res$Lb[8], equals(6005/26214))
  expect_that(res$Lb[9], equals(5548/24822))
  expect_that(res$Lb[10], equals(4975/24822))
  expect_that(res$Lb[11], equals(7320/24822))
  expect_that(res$Lb[12], equals(6979/24822))
  expect_that(res$Lb[13], equals(16865/50136))
  expect_that(res$Lb[14], equals(10906/50136))
  expect_that(res$Lb[15], equals(10649/50136))
  expect_that(res$Lb[16], equals(11716/50136))
  
  # Check result: Total peak height.
  expect_that(res$TPH[1], equals(7648))
  expect_that(res$TPH[2], equals(7911))
  expect_that(res$TPH[3], equals(6774))
  expect_that(res$TPH[4], equals(4829))
  expect_that(res$TPH[5], equals(5482))
  expect_that(res$TPH[6], equals(8098))
  expect_that(res$TPH[7], equals(6629))
  expect_that(res$TPH[8], equals(6005))
  expect_that(res$TPH[9], equals(5548))
  expect_that(res$TPH[10], equals(4975))
  expect_that(res$TPH[11], equals(7320))
  expect_that(res$TPH[12], equals(6979))
  expect_that(res$TPH[13], equals(16865))
  expect_that(res$TPH[14], equals(10906))
  expect_that(res$TPH[15], equals(10649))
  expect_that(res$TPH[16], equals(11716))
  
  # Check result: Total profile peak height.
  expect_that(unique(res$TPPH), equals(c(27162,26214,24822,50136)))
  
  
  # TEST 10 -------------------------------------------------------------------
  # Test that proportional Lb by dye can be calculated for filtered data.
  
  # Analyse dataframe.
  res <- calculateLb(data = set1, ref = ref1, option = "prop", by.dye = TRUE,
                     ol.rm = FALSE, sex.rm = FALSE, na = NULL,
                     ignore.case = TRUE, word = FALSE, exact = FALSE)

  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_false(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$TPH))
  expect_false(is.null(res$Peaks))
  expect_false(is.null(res$Dye))
  expect_false(is.null(res$TPPH))
  expect_false(is.null(res$Lb))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$TPH)))
  expect_false(any(is.na(res$Peaks)))
  expect_false(any(is.na(res$Dye)))
  expect_false(any(is.na(res$TPPH)))
  expect_false(any(is.na(res$Lb)))
  
  # Check result: Locus balance.
  expect_that(res$Lb[1], equals(7131/33079))
  expect_that(res$Lb[2], equals(7308/33079))
  expect_that(res$Lb[3], equals(7911/33079))
  expect_that(res$Lb[4], equals(6282/33079))
  expect_that(res$Lb[5], equals(4447/33079))
  expect_that(res$Lb[6], equals(5033/24651))
  expect_that(res$Lb[7], equals(7847/24651))
  expect_that(res$Lb[8], equals(6088/24651))
  expect_that(res$Lb[9], equals(5683/24651))
  expect_that(res$Lb[10], equals(4757/22933))
  expect_that(res$Lb[11], equals(4573/22933))
  expect_that(res$Lb[12], equals(7005/22933))
  expect_that(res$Lb[13], equals(6598/22933))
  expect_that(res$Lb[14], equals(16099/47097))
  expect_that(res$Lb[15], equals(9805/47097))
  expect_that(res$Lb[16], equals(10287/47097))
  expect_that(res$Lb[17], equals(10906/47097))
  
  # Check result: Total peak height.
  expect_that(res$TPH[1], equals(7131))
  expect_that(res$TPH[2], equals(7308))
  expect_that(res$TPH[3], equals(7911))
  expect_that(res$TPH[4], equals(6282))
  expect_that(res$TPH[5], equals(4447))
  expect_that(res$TPH[6], equals(5033))
  expect_that(res$TPH[7], equals(7847))
  expect_that(res$TPH[8], equals(6088))
  expect_that(res$TPH[9], equals(5683))
  expect_that(res$TPH[10], equals(4757))
  expect_that(res$TPH[11], equals(4573))
  expect_that(res$TPH[12], equals(7005))
  expect_that(res$TPH[13], equals(6598))
  expect_that(res$TPH[14], equals(16099))
  expect_that(res$TPH[15], equals(9805))
  expect_that(res$TPH[16], equals(10287))
  expect_that(res$TPH[17], equals(10906))
  
  # Check result: Total profile peak height.
  expect_that(unique(res$TPPH), equals(c(33079,24651,22933,47097)))
  

  # TEST 11 -------------------------------------------------------------------
  # Test that proportional Lb by dye can be calculated for filtered data
  # and off-ladder peaks removed.
  
  # Analyse dataframe.
  res <- calculateLb(data = set1, ref = ref1, option = "prop", by.dye = TRUE,
                     ol.rm = TRUE, sex.rm = FALSE, na = NULL,
                     ignore.case = TRUE, word = FALSE, exact = FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_false(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$TPH))
  expect_false(is.null(res$Peaks))
  expect_false(is.null(res$Dye))
  expect_false(is.null(res$TPPH))
  expect_false(is.null(res$Lb))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$TPH)))
  expect_false(any(is.na(res$Peaks)))
  expect_false(any(is.na(res$Dye)))
  expect_false(any(is.na(res$TPPH)))
  expect_false(any(is.na(res$Lb)))
  
  # Check result: Locus balance.
  expect_that(res$Lb[1], equals(7131/33079))
  expect_that(res$Lb[2], equals(7308/33079))
  expect_that(res$Lb[3], equals(7911/33079))
  expect_that(res$Lb[4], equals(6282/33079))
  expect_that(res$Lb[5], equals(4447/33079))
  expect_that(res$Lb[6], equals(5033/24651))
  expect_that(res$Lb[7], equals(7847/24651))
  expect_that(res$Lb[8], equals(6088/24651))
  expect_that(res$Lb[9], equals(5683/24651))
  expect_that(res$Lb[10], equals(4757/22933))
  expect_that(res$Lb[11], equals(4573/22933))
  expect_that(res$Lb[12], equals(7005/22933))
  expect_that(res$Lb[13], equals(6598/22933))
  expect_that(res$Lb[14], equals(16099/47097))
  expect_that(res$Lb[15], equals(9805/47097))
  expect_that(res$Lb[16], equals(10287/47097))
  expect_that(res$Lb[17], equals(10906/47097))
  
  # Check result: Total peak height.
  expect_that(res$TPH[1], equals(7131))
  expect_that(res$TPH[2], equals(7308))
  expect_that(res$TPH[3], equals(7911))
  expect_that(res$TPH[4], equals(6282))
  expect_that(res$TPH[5], equals(4447))
  expect_that(res$TPH[6], equals(5033))
  expect_that(res$TPH[7], equals(7847))
  expect_that(res$TPH[8], equals(6088))
  expect_that(res$TPH[9], equals(5683))
  expect_that(res$TPH[10], equals(4757))
  expect_that(res$TPH[11], equals(4573))
  expect_that(res$TPH[12], equals(7005))
  expect_that(res$TPH[13], equals(6598))
  expect_that(res$TPH[14], equals(16099))
  expect_that(res$TPH[15], equals(9805))
  expect_that(res$TPH[16], equals(10287))
  expect_that(res$TPH[17], equals(10906))
  
  # Check result: Total profile peak height.
  expect_that(unique(res$TPPH), equals(c(33079,24651,22933,47097)))
  
  
  # TEST 12 -------------------------------------------------------------------
  # Test that proportional Lb by dye can be calculated for filtered data,
  # with off-ladder peaks, and sex markers removed.
  
  # Analyse dataframe.
  res <- calculateLb(data = set1, ref = ref1, option = "prop", by.dye = TRUE,
                     ol.rm = TRUE, sex.rm = TRUE, na = NULL,
                     ignore.case = TRUE, word = FALSE, exact = FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_false(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$TPH))
  expect_false(is.null(res$Peaks))
  expect_false(is.null(res$Dye))
  expect_false(is.null(res$TPPH))
  expect_false(is.null(res$Lb))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$TPH)))
  expect_false(any(is.na(res$Peaks)))
  expect_false(any(is.na(res$Dye)))
  expect_false(any(is.na(res$TPPH)))
  expect_false(any(is.na(res$Lb)))
  
  # Check result: Locus balance.
  expect_that(res$Lb[1], equals(7308/25948))
  expect_that(res$Lb[2], equals(7911/25948))
  expect_that(res$Lb[3], equals(6282/25948))
  expect_that(res$Lb[4], equals(4447/25948))
  expect_that(res$Lb[5], equals(5033/24651))
  expect_that(res$Lb[6], equals(7847/24651))
  expect_that(res$Lb[7], equals(6088/24651))
  expect_that(res$Lb[8], equals(5683/24651))
  expect_that(res$Lb[9], equals(4757/22933))
  expect_that(res$Lb[10], equals(4573/22933))
  expect_that(res$Lb[11], equals(7005/22933))
  expect_that(res$Lb[12], equals(6598/22933))
  expect_that(res$Lb[13], equals(16099/47097))
  expect_that(res$Lb[14], equals(9805/47097))
  expect_that(res$Lb[15], equals(10287/47097))
  expect_that(res$Lb[16], equals(10906/47097))
  
  # Check result: Total peak height.
  expect_that(res$TPH[1], equals(7308))
  expect_that(res$TPH[2], equals(7911))
  expect_that(res$TPH[3], equals(6282))
  expect_that(res$TPH[4], equals(4447))
  expect_that(res$TPH[5], equals(5033))
  expect_that(res$TPH[6], equals(7847))
  expect_that(res$TPH[7], equals(6088))
  expect_that(res$TPH[8], equals(5683))
  expect_that(res$TPH[9], equals(4757))
  expect_that(res$TPH[10], equals(4573))
  expect_that(res$TPH[11], equals(7005))
  expect_that(res$TPH[12], equals(6598))
  expect_that(res$TPH[13], equals(16099))
  expect_that(res$TPH[14], equals(9805))
  expect_that(res$TPH[15], equals(10287))
  expect_that(res$TPH[16], equals(10906))
  
  # Check result: Total profile peak height.
  expect_that(unique(res$TPPH), equals(c(25948,24651,22933,47097)))
  
  
  # TEST 13 -------------------------------------------------------------------
  # Test when a marker is missing from a sample.

  # Remove TH01 from the first sample.
  setMissing <- set2[!(set2$Marker == "TH01" & set2$Sample.Name == "SampleA01"),]

  # Test that proportional Lb by dye can be calculated for filtered data,
  # with off-ladder peaks removed and markers are missing.
  
  # This work since ref is provided (data is filtered and missing markers added).
  res <- calculateLb(data = setMissing, ref = ref2, option = "prop", by.dye = TRUE,
                     ol.rm = TRUE, sex.rm = FALSE, na = NULL,
                     ignore.case = TRUE, word = FALSE, exact = FALSE)

  # Same as above but global Lb.
  res <- calculateLb(data = setMissing, ref = ref2, option = "prop", by.dye = FALSE,
                     ol.rm = TRUE, sex.rm = FALSE, na = NULL,
                     ignore.case = TRUE, word = FALSE, exact = FALSE)
  
  # This will not work since one marker is missing in the first sample.  
  expect_that(calculateLb(data = setMissing, ref = NULL, option = "prop", by.dye = TRUE,
                     ol.rm = TRUE, sex.rm = FALSE, na = NULL,
                     ignore.case = TRUE, word = FALSE, exact = FALSE),
              throws_error())
  
  # Same as above but global Lb.
  expect_that(calculateLb(data = setMissing, ref = NULL, option = "prop", by.dye = FALSE,
                          ol.rm = TRUE, sex.rm = FALSE, na = NULL,
                          ignore.case = TRUE, word = FALSE, exact = FALSE),
              throws_error())
  
  
  # TEST 14 -------------------------------------------------------------------
  # Test when all markers in one dye channel is missing from a sample.
  
  # Remove yellow dye channel from the first sample.
  setMissing <- set2[!(set2$Dye == "Y" & set2$Sample.Name == "SampleA01"),]
  
  # Test that proportional Lb by dye can be calculated for filtered data,
  # with off-ladder peaks removed and markers are missing.
  
  # This work since ref is provided (data is filtered and missing markers added).
  res <- calculateLb(data = setMissing, ref = ref2, option = "prop", by.dye = TRUE,
                     ol.rm = TRUE, sex.rm = FALSE, na = NULL,
                     ignore.case = TRUE, word = FALSE, exact = FALSE)
  
  # Same as above but global Lb.
  res <- calculateLb(data = setMissing, ref = ref2, option = "prop", by.dye = FALSE,
                     ol.rm = TRUE, sex.rm = FALSE, na = NULL,
                     ignore.case = TRUE, word = FALSE, exact = FALSE)
  
  # This will not work since one marker is missing in the first sample.  
  expect_that(calculateLb(data = setMissing, ref = NULL, option = "prop", by.dye = TRUE,
                          ol.rm = TRUE, sex.rm = FALSE, na = NULL,
                          ignore.case = TRUE, word = FALSE, exact = FALSE),
              throws_error())
  
  # Same as above but global Lb.
  expect_that(calculateLb(data = setMissing, ref = NULL, option = "prop", by.dye = FALSE,
                          ol.rm = TRUE, sex.rm = FALSE, na = NULL,
                          ignore.case = TRUE, word = FALSE, exact = FALSE),
              throws_error())
  
  
  # TEST 15 -------------------------------------------------------------------
  # Test that 'fat' data throws error.
  
  # Expect error due to 'fat' data.
  expect_error(calculateLb(data = set1Fat))
  
})