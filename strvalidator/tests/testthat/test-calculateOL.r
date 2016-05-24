context("calculateOL")

################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG
# 20.01.2014: All tests working.
# 
# test_dir("inst/tests/")
# test_file("tests/testthat/test-calculateOL.r")
# test_dir("tests/testthat")

test_that("calculateOL", {

  # Load test data.
  database <- getDb("ESX 17 Hill")
  kitInfo <- getKit("ESX17")

  # Load test data.
  database2 <- getDb("SGM Norway")
  kitInfo2 <- getKit("SGMPlus")
  
  # TEST 01 -------------------------------------------------------------------
  # Test excluding virtual alleles, actual frequencies.
  
  # Analyse dataframe.
  res <- calculateOL(kit=kitInfo, db=database,
                     virtual=FALSE, limit=FALSE, debug=FALSE)

  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))

  # Check that expected columns exist.  
  expect_true(any(grepl("Kit", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Database", names(res))))
  expect_true(any(grepl("Risk", names(res))))
  expect_true(any(grepl("Total", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Kit)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Database)))
  expect_false(any(is.na(res$Risk)))
  expect_false(any(is.na(res$Total)))
  
  # Check result.
  expect_that(ncol(res), equals(5))
  expect_that(nrow(res), equals(17))
  
  # Check result.
  expect_that(length(unique(res$Total)), equals(1))
  expect_that(round(unique(res$Total),4), equals(0.0609))
  
  
  # TEST 02 -------------------------------------------------------------------
  # Test including virtual alleles, actual frequencies.
  
  # Analyse dataframe.
  res <- calculateOL(kit=kitInfo, db=database,
                     virtual=TRUE, limit=FALSE, debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Kit", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Database", names(res))))
  expect_true(any(grepl("Risk", names(res))))
  expect_true(any(grepl("Total", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Kit)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Database)))
  expect_false(any(is.na(res$Risk)))
  expect_false(any(is.na(res$Total)))
  
  # Check result.
  expect_that(ncol(res), equals(5))
  expect_that(nrow(res), equals(17))
  
  # Check result.
  expect_that(length(unique(res$Total)), equals(1))
  expect_that(round(unique(res$Total),4), equals(0.0185))
  
  # TEST 03 -------------------------------------------------------------------
  # Test excluding virtual alleles, using minimum frequency.
  
  # Analyse dataframe.
  res <- calculateOL(kit=kitInfo, db=database,
                     virtual=FALSE, limit=TRUE, debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Kit", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Database", names(res))))
  expect_true(any(grepl("Risk", names(res))))
  expect_true(any(grepl("Total", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Kit)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Database)))
  expect_false(any(is.na(res$Risk)))
  expect_false(any(is.na(res$Total)))
  
  # Check result.
  expect_that(ncol(res), equals(5))
  expect_that(nrow(res), equals(17))
  
  # Check result.
  expect_that(length(unique(res$Total)), equals(1))
  expect_that(round(unique(res$Total),4), equals(0.1096))
  
  # TEST 04 -------------------------------------------------------------------
  # Test including virtual alleles, using minimum frequency.
  
  # Analyse dataframe.
  res <- calculateOL(kit=kitInfo, db=database,
                     virtual=TRUE, limit=TRUE, debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Kit", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Database", names(res))))
  expect_true(any(grepl("Risk", names(res))))
  expect_true(any(grepl("Total", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Kit)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Database)))
  expect_false(any(is.na(res$Risk)))
  expect_false(any(is.na(res$Total)))
  
  # Check result.
  expect_that(ncol(res), equals(5))
  expect_that(nrow(res), equals(17))
  
  # Check result.
  expect_that(length(unique(res$Total)), equals(1))
  expect_that(round(unique(res$Total),4), equals(0.0392))

  
  # TEST 05 -------------------------------------------------------------------
  # Test excluding virtual alleles, actual frequencies.
  
  # Analyse dataframe.
  res <- calculateOL(kit=kitInfo2, db=database2,
                     virtual=FALSE, limit=FALSE, debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Kit", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Database", names(res))))
  expect_true(any(grepl("Risk", names(res))))
  expect_true(any(grepl("Total", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Kit)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Database)))
  expect_false(any(is.na(res$Risk)))
  expect_false(any(is.na(res$Total)))
  
  # Check result.
  expect_that(ncol(res), equals(5))
  expect_that(nrow(res), equals(11))
  
  # Check result.
  expect_that(length(unique(res$Total)), equals(1))
  expect_that(round(unique(res$Total),4), equals(0))
  
  
  # TEST 06 -------------------------------------------------------------------
  # Test including virtual alleles, actual frequencies.
  
  # Analyse dataframe.
  res <- calculateOL(kit=kitInfo2, db=database2,
                     virtual=TRUE, limit=FALSE, debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Kit", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Database", names(res))))
  expect_true(any(grepl("Risk", names(res))))
  expect_true(any(grepl("Total", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Kit)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Database)))
  expect_false(any(is.na(res$Risk)))
  expect_false(any(is.na(res$Total)))
  
  # Check result.
  expect_that(ncol(res), equals(5))
  expect_that(nrow(res), equals(11))
  
  # Check result.
  expect_that(length(unique(res$Total)), equals(1))
  expect_that(round(unique(res$Total),4), equals(0))
  
  # TEST 07 -------------------------------------------------------------------
  # Test excluding virtual alleles, using minimum frequency.
  
  # Analyse dataframe.
  res <- calculateOL(kit=kitInfo2, db=database2,
                     virtual=FALSE, limit=TRUE, debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Kit", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Database", names(res))))
  expect_true(any(grepl("Risk", names(res))))
  expect_true(any(grepl("Total", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Kit)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Database)))
  expect_false(any(is.na(res$Risk)))
  expect_false(any(is.na(res$Total)))
  
  # Check result.
  expect_that(ncol(res), equals(5))
  expect_that(nrow(res), equals(11))
  
  # Check result.
  expect_that(length(unique(res$Total)), equals(1))
  expect_that(round(unique(res$Total),4), equals(0))
  
  # TEST 08 -------------------------------------------------------------------
  # Test including virtual alleles, using minimum frequency.
  
  # Analyse dataframe.
  res <- calculateOL(kit=kitInfo2, db=database2,
                     virtual=TRUE, limit=TRUE, debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Kit", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Database", names(res))))
  expect_true(any(grepl("Risk", names(res))))
  expect_true(any(grepl("Total", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Kit)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Database)))
  expect_false(any(is.na(res$Risk)))
  expect_false(any(is.na(res$Total)))
  
  # Check result.
  expect_that(ncol(res), equals(5))
  expect_that(nrow(res), equals(11))
  
  # Check result.
  expect_that(length(unique(res$Total)), equals(1))
  expect_that(round(unique(res$Total),4), equals(0))
  
  
})