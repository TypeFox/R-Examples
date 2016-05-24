context("trim")

################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG
# 28.04.2014: First tests for 'trim'.
# 
# test_dir("inst/tests/")
# test_file("tests/testthat/test-trim.r")
# test_dir("tests/testthat")

test_that("trim", {

  # Columns:
  allele <- c("10","","X","Y","12","OL")
  nastrcol <- c("NA","NA","NA","NA","NA","NA")
  nacol <- c(NA,NA,NA,NA,NA,NA)
  emptycol <- c("","","","","","")

  # Samples:
  samples1 <- c("01-Positive Control", "02-Negative Control", "03-Sample.1",
                "04-Sample.2", "05-Sample.3", "06-Allelic Ladder")
  
  samples2 <- c("C+", "C-", "C", "C1", "C2", "C3",
                "c+", "c-", "c", "c4", "c5", "c6")
  
  # Create a dataframe for testing:
  
  # Test set 1.
  df1 <- data.frame(Sample.Name=samples1, Allele=allele,
                    NAstr=nastrcol, NAcol=nacol, Empty=emptycol,
                    stringsAsFactors=FALSE)
  
  # Test set 2.
  df2 <- data.frame(Sample.Name=samples2, Allele=allele,
                    NAstr=nastrcol, NAcol=nacol, Empty=emptycol,
                    stringsAsFactors=FALSE)
  

  
  # TEST 01 -------------------------------------------------------------------
  # Test nothing trimmed or changed using df1.
  
  # Analyse dataframe.
  res <- trim(data=df1, samples=NULL, columns=NULL, word=FALSE,
              ignore.case=TRUE, invert.s=FALSE, invert.c=FALSE,
              rm.na.col=FALSE, rm.empty.col=FALSE, missing=NULL, debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))

  # Check dimensions.
  expect_true(ncol(res) == 5)
  expect_true(nrow(res) == 6)
  
  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res))
  expect_true("Allele" %in% names(res))
  expect_true("NAstr" %in% names(res))
  expect_true("NAcol" %in% names(res))
  expect_true("Empty" %in% names(res))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Allele)))
  expect_false(any(is.na(res$NAstr)))
  expect_true(any(is.na(res$NAcol)))
  expect_false(any(is.na(res$Empty)))
  
  # TEST 02 -------------------------------------------------------------------
  # Test remove empty columns.
  
  # Analyse dataframe.
  res <- trim(data=df1, samples=NULL, columns=NULL, word=FALSE,
              ignore.case=TRUE, invert.s=FALSE, invert.c=FALSE,
              rm.na.col=FALSE, rm.empty.col=TRUE, missing=NULL, debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check dimensions.
  expect_true(ncol(res) == 4)
  expect_true(nrow(res) == 6)
  
  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res))
  expect_true("Allele" %in% names(res))
  expect_true("NAstr" %in% names(res))
  expect_true("NAcol" %in% names(res))
  expect_false("Empty" %in% names(res))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Allele)))
  expect_false(any(is.na(res$NAstr)))
  expect_true(any(is.na(res$NAcol)))

  # TEST 03 -------------------------------------------------------------------
  # Test remove NA columns.
  
  # Analyse dataframe.
  res <- trim(data=df1, samples=NULL, columns=NULL, word=FALSE,
              ignore.case=TRUE, invert.s=FALSE, invert.c=FALSE,
              rm.na.col=TRUE, rm.empty.col=FALSE, missing=NULL, debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check dimensions.
  expect_true(ncol(res) == 4)
  expect_true(nrow(res) == 6)
  
  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res))
  expect_true("Allele" %in% names(res))
  expect_true("NAstr" %in% names(res))
  expect_false("NAcol" %in% names(res))
  expect_true("Empty" %in% names(res))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Allele)))
  expect_false(any(is.na(res$NAstr)))
  expect_false(any(is.na(res$Empty)))
  
  # TEST 04 -------------------------------------------------------------------
  # Test replace missing values with NA.
  
  # Analyse dataframe.
  res <- trim(data=df1, samples=NULL, columns=NULL, word=FALSE,
              ignore.case=TRUE, invert.s=FALSE, invert.c=FALSE,
              rm.na.col=FALSE, rm.empty.col=FALSE, missing=NA, debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check dimensions.
  expect_true(ncol(res) == 5)
  expect_true(nrow(res) == 6)
  
  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res))
  expect_true("Allele" %in% names(res))
  expect_true("NAstr" %in% names(res))
  expect_true("NAcol" %in% names(res))
  expect_true("Empty" %in% names(res))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_true(any(is.na(res$Allele)))
  expect_false(any(is.na(res$NAstr)))
  expect_true(any(is.na(res$NAcol)))
  expect_true(any(is.na(res$Empty)))
  
  # TEST 05 -------------------------------------------------------------------
  # Test replace missing values with string.
  
  # Analyse dataframe.
  res <- trim(data=df1, samples=NULL, columns=NULL, word=FALSE,
              ignore.case=TRUE, invert.s=FALSE, invert.c=FALSE,
              rm.na.col=FALSE, rm.empty.col=FALSE, missing="N/A", debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check dimensions.
  expect_true(ncol(res) == 5)
  expect_true(nrow(res) == 6)
  
  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res))
  expect_true("Allele" %in% names(res))
  expect_true("NAstr" %in% names(res))
  expect_true("NAcol" %in% names(res))
  expect_true("Empty" %in% names(res))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Allele)))
  expect_false(any(is.na(res$NAstr)))
  expect_true(any(is.na(res$NAcol)))
  expect_false(any(is.na(res$Empty)))
  
  # TEST 06 -------------------------------------------------------------------
  # Test to trim one samples and columns using df1.
  
  # Analyse dataframe.
  res <- trim(data=df1, samples=c("Sample.1"), columns=c("Sample.Name"), word=FALSE,
              ignore.case=TRUE, invert.s=FALSE, invert.c=FALSE,
              rm.na.col=FALSE, rm.empty.col=FALSE, missing=NA, debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check dimensions.
  expect_true(ncol(res) == 1)
  expect_true(nrow(res) == 1)
  
  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res))
  expect_false("Allele" %in% names(res))
  expect_false("NAstr" %in% names(res))
  expect_false("NAcol" %in% names(res))
  expect_false("Empty" %in% names(res))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))

  # TEST 07 -------------------------------------------------------------------
  # Test to trim using invert for sample and column with df1.
  
  # Analyse dataframe.
  res <- trim(data=df1, samples=c("Sample.1"), columns=c("Sample.Name"), word=FALSE,
              ignore.case=TRUE, invert.s=TRUE, invert.c=TRUE,
              rm.na.col=FALSE, rm.empty.col=FALSE, missing=NA, debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check dimensions.
  expect_true(ncol(res) == 4)
  expect_true(nrow(res) == 5)
  
  # Check that expected columns exist.  
  expect_false("Sample.Name" %in% names(res))
  expect_true("Allele" %in% names(res))
  expect_true("NAstr" %in% names(res))
  expect_true("NAcol" %in% names(res))
  expect_true("Empty" %in% names(res))
  
  # Check for NA's.
  expect_true(any(is.na(res$Allele)))
  expect_false(any(is.na(res$NAstr)))
  expect_true(any(is.na(res$NAcol)))
  expect_true(any(is.na(res$Empty)))
  
  # TEST 08 -------------------------------------------------------------------
  # Test to trim controls and columns using df1.
  
  # Analyse dataframe.
  res <- trim(data=df1, samples=c("pos|neg|ladder"), columns=c("Sample.Name|Allele"), word=FALSE,
              ignore.case=TRUE, invert.s=TRUE, invert.c=FALSE,
              rm.na.col=FALSE, rm.empty.col=FALSE, missing=NA, debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check dimensions.
  expect_true(ncol(res) == 2)
  expect_true(nrow(res) == 3)
  
  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res))
  expect_true("Allele" %in% names(res))
  expect_false("NAstr" %in% names(res))
  expect_false("NAcol" %in% names(res))
  expect_false("Empty" %in% names(res))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Allele)))
  
  # TEST 09 -------------------------------------------------------------------
  # Test to trim controls and columns with +/-, and respect case using df2.
  # ...also invert columns, remove NA/empty columns.
  
  # Analyse dataframe.
  res <- trim(data=df2, samples=c("c+|c-"), columns=c("NAstr"), word=FALSE,
              ignore.case=FALSE, invert.s=TRUE, invert.c=TRUE,
              rm.na.col=TRUE, rm.empty.col=TRUE, missing=NA, debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check dimensions.
  expect_true(ncol(res) == 2)
  expect_true(nrow(res) == 10)
  
  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res))
  expect_true("Allele" %in% names(res))
  expect_false("NAstr" %in% names(res))
  expect_false("NAcol" %in% names(res))
  expect_false("Empty" %in% names(res))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_true(any(is.na(res$Allele)))
  
  # TEST 10 -------------------------------------------------------------------
  # Test to trim controls and columns with +/-, and ignore case using df2.
  # ...also invert columns, remove NA/empty columns.
  
  # Analyse dataframe.
  res <- trim(data=df2, samples=c("c+|c-"), columns=c("NAstr"), word=FALSE,
              ignore.case=TRUE, invert.s=TRUE, invert.c=TRUE,
              rm.na.col=TRUE, rm.empty.col=TRUE, missing=NA, debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check dimensions.
  expect_true(ncol(res) == 2)
  expect_true(nrow(res) == 8)
  
  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res))
  expect_true("Allele" %in% names(res))
  expect_false("NAstr" %in% names(res))
  expect_false("NAcol" %in% names(res))
  expect_false("Empty" %in% names(res))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Allele)))
  
  # TEST 11 -------------------------------------------------------------------
  # Test word boundary using df2.
  
  # Analyse dataframe.
  res <- trim(data=df2, samples=c("c"), columns=c("NAstr"), word=TRUE,
              ignore.case=TRUE, invert.s=FALSE, invert.c=TRUE,
              rm.na.col=FALSE, rm.empty.col=FALSE, missing="*", debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check dimensions.
  expect_true(ncol(res) == 4)
  expect_true(nrow(res) == 6)
  
  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res))
  expect_true("Allele" %in% names(res))
  expect_false("NAstr" %in% names(res))
  expect_true("NAcol" %in% names(res))
  expect_true("Empty" %in% names(res))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Allele)))
  expect_true(any(is.na(res$NAcol)))
  expect_false(any(is.na(res$Empty)))
  
})