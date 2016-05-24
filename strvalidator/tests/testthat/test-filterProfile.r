context("filterProfile")

################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG
# 15.12.2015: Added test 10-12 to test new option 'exact'.
# 09.04.2015: Added test 09 to test new option 'invert'.
# 22.01.2014: Fixed test 2 in accordance with change in 'filterProfile'.
# 20.01.2014: All tests working.
# 
# require(testthat)
# test_dir("inst/tests/")
# test_file("tests/testthat/test-filterProfile.r")
# test_dir("tests/testthat")

test_that("filterProfile", {

  # Load test data.
  data(set1)
  data(ref1)
  data(set2)
  data(ref2)
  
  # Dataset with missing marker.
  set3 <- set2[set2$Marker!="TH01",]
  
  # Dataset with different case in sample name.
  set4 <- set2
  set4$Sample.Name[set4$Sample.Name=="SampleA02"] <- "sAmPLea02"
  
  # Dataset with dirty' data.
  set6 <- slim(set1, fix=c("Sample.Name", "Marker"), 
               stack=c("Allele", "Height"),
               keep.na=TRUE)
  ref6 <- slim(ref1, fix=c("Sample.Name", "Marker"),
               stack=c("Allele"),
               keep.na=TRUE)
  
  # Dataset with double notation (16/16 vs 16).
  set7 <- slim(set1, fix=c("Sample.Name", "Marker"), stack=c("Allele", "Height"))
  ref7 <- slim(ref1, fix=c("Sample.Name", "Marker"), stack=c("Allele"))
  sampleName7 <- ref7$Sample.Name
  marker7 <- ref7$Marker
  allele7 <- ref7$Allele
  # Duplicate homozygous alleles.
  allele7 <- c(allele7[1:19], allele7[19], allele7[20:33])
  marker7 <- c(marker7[1:19],marker7[19],marker7[20:33])
  sampleName7 <- c(sampleName7[1], sampleName7)
  # Make new dataframe.
  ref7 <- data.frame(Sample.Name=sampleName7, Marker=marker7, Allele=allele7, stringsAsFactors=FALSE)
  
  # Dataset without Sample.Name column.
  set8 <- slim(set1, fix=c("Sample.Name", "Marker"), stack=c("Allele", "Height"))
  ref8 <- slim(ref1, fix=c("Sample.Name", "Marker"), stack=c("Allele"))
  ref8 <- ref8[,c("Marker","Allele")]
  
  # Dataset testing exact matching.
  # 'PC1' is in 'NOPC15'
  set9 <- slim(set1, fix=c("Sample.Name", "Marker"), stack=c("Allele", "Height"))
  set9 <- trim(data = set9, samples = "PC1|ladder")
  set9[set9$Sample.Name == "Ladder", ]$Sample.Name <- "NOPC15"
  ref9 <- guessProfile(data = set9, ol.rm = TRUE)
  ref9$Height <- as.character(ref9$Height) # To enable comparison using 'all.equal'.
  
  # TEST 01 -------------------------------------------------------------------
  # Test that filtering of two 'clean' samples work.
  
  # Analyse dataframe.
  res <- filterProfile(data=set2, ref=ref2,
                       add.missing.loci=FALSE,
                       keep.na=FALSE,
                       ignore.case=FALSE,
                       debug=FALSE)

  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))

  # Check that expected columns exist.  
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  expect_true(any(grepl("Dye", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Allele)))
  expect_false(any(is.na(res$Height)))
  expect_false(any(is.na(res$Dye)))
  
  # Check result.
  expect_that(ncol(res), equals(5))
  expect_that(nrow(res), equals(31))

  # TEST 02 -------------------------------------------------------------------
  # Test that adding missing loci work.

  # Analyse dataframe.
  res <- filterProfile(data=set3, ref=ref2,
                       add.missing.loci=TRUE,
                       keep.na=TRUE,
                       ignore.case=FALSE,
                       debug=FALSE)
  
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_true(any(is.na(res$Allele)))
  expect_true(any(is.na(res$Height)))
  
  # Check result.
  expect_that(ncol(res), equals(5))
  expect_that(nrow(res), equals(30))

  # TEST 03 -------------------------------------------------------------------
  # Test that keep NA loci work.
  
  # Analyse dataframe.
  res <- filterProfile(data=set2, ref=ref2,
                       add.missing.loci=FALSE,
                       keep.na=TRUE,
                       ignore.case=FALSE,
                       debug=FALSE)
  
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_true(any(is.na(res$Allele)))
  expect_true(any(is.na(res$Height)))
  
  # Check result.
  expect_that(ncol(res), equals(5))
  expect_that(nrow(res), equals(32))

  # TEST 04 -------------------------------------------------------------------
  # Test that ignore case = FALSE work.

  # Analyse dataframe.
  res <- filterProfile(data=set2, ref=ref2,
                       add.missing.loci=FALSE,
                       keep.na=TRUE,
                       ignore.case=FALSE,
                       debug=FALSE)
  
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_true(any(is.na(res$Allele)))
  expect_true(any(is.na(res$Height)))
  
  # Check result.
  expect_that(ncol(res), equals(5))
  expect_that(nrow(res), equals(32))
  
  # TEST 05 -------------------------------------------------------------------
  # Test that 'fat' data gives error.
  
  # Test using 'fat' data.
  expect_that(filterProfile(data=set1, ref=ref2,
                       add.missing.loci=FALSE,
                       keep.na=TRUE,
                       ignore.case=FALSE,
                       debug=FALSE),
              throws_error())
  
  # Test using 'fat' ref.
  expect_that(filterProfile(data=set2, ref=ref1,
                            add.missing.loci=FALSE,
                            keep.na=TRUE,
                            ignore.case=FALSE,
                            debug=FALSE),
              throws_error())
  
  # Test using 'fat' data and ref.
  expect_that(filterProfile(data=set1, ref=ref1,
                            add.missing.loci=FALSE,
                            keep.na=TRUE,
                            ignore.case=FALSE,
                            debug=FALSE),
              throws_error())
  
  # TEST 06 -------------------------------------------------------------------
  # Test that filtering of 'dirty' samples work.

  # Analyse dataframe.
  res <- filterProfile(data=set6, ref=ref6,
                       add.missing.loci=FALSE,
                       keep.na=FALSE,
                       ignore.case=FALSE,
                       debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Allele)))
  expect_false(any(is.na(res$Height)))
  
  # Check result.
  expect_that(ncol(res), equals(4))
  expect_that(nrow(res), equals(264))
  
  # Check success of filtering.
  expect_true(all(unique(res$Allele[res$Marker=="AMEL"]) == c("X","Y")))
  expect_true(all(unique(res$Allele[res$Marker=="D3S1358"]) == c("17","18")))
  expect_true(all(unique(res$Allele[res$Marker=="TH01"]) == c("6","9.3")))
  expect_true(all(unique(res$Allele[res$Marker=="D21S11"]) == c("29","31.2")))
  expect_true(all(unique(res$Allele[res$Marker=="D18S51"]) == c("16","18")))
  expect_true(all(unique(res$Allele[res$Marker=="D10S1248"]) == c("13","15")))
  expect_true(all(unique(res$Allele[res$Marker=="D1S1656"]) == c("12","13")))
  expect_true(all(unique(res$Allele[res$Marker=="D2S1338"]) == c("22","25")))
  expect_true(all(unique(res$Allele[res$Marker=="D16S539"]) == c("9","13")))
  expect_true(all(unique(res$Allele[res$Marker=="D22S1045"]) == c("16")))
  expect_true(all(unique(res$Allele[res$Marker=="vWA"]) == c("16","19")))
  expect_true(all(unique(res$Allele[res$Marker=="D8S1179"]) == c("14","15")))
  expect_true(all(unique(res$Allele[res$Marker=="FGA"]) == c("20","23")))
  expect_true(all(unique(res$Allele[res$Marker=="D2S441"]) == c("10","14")))
  expect_true(all(unique(res$Allele[res$Marker=="D12S391"]) == c("18","23")))
  expect_true(all(unique(res$Allele[res$Marker=="D19S433"]) == c("13","14")))
  expect_true(all(unique(res$Allele[res$Marker=="SE33"]) == c("15","16")))

  # TEST 07 -------------------------------------------------------------------
  # Test that filtering using double notation work (16/16 vs 16).
  
  # Analyse dataframe.
  res <- filterProfile(data=set7, ref=ref7,
                       add.missing.loci=FALSE,
                       keep.na=FALSE,
                       ignore.case=FALSE,
                       debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Allele)))
  expect_false(any(is.na(res$Height)))
  
  # Check result.
  expect_that(ncol(res), equals(4))
  expect_that(nrow(res), equals(264))
  
  # Check success of filtering.
  expect_true(all(unique(res$Allele[res$Marker=="AMEL"]) == c("X","Y")))
  expect_true(all(unique(res$Allele[res$Marker=="D3S1358"]) == c("17","18")))
  expect_true(all(unique(res$Allele[res$Marker=="TH01"]) == c("6","9.3")))
  expect_true(all(unique(res$Allele[res$Marker=="D21S11"]) == c("29","31.2")))
  expect_true(all(unique(res$Allele[res$Marker=="D18S51"]) == c("16","18")))
  expect_true(all(unique(res$Allele[res$Marker=="D10S1248"]) == c("13","15")))
  expect_true(all(unique(res$Allele[res$Marker=="D1S1656"]) == c("12","13")))
  expect_true(all(unique(res$Allele[res$Marker=="D2S1338"]) == c("22","25")))
  expect_true(all(unique(res$Allele[res$Marker=="D16S539"]) == c("9","13")))
  expect_true(all(unique(res$Allele[res$Marker=="D22S1045"]) == c("16")))
  expect_true(all(unique(res$Allele[res$Marker=="vWA"]) == c("16","19")))
  expect_true(all(unique(res$Allele[res$Marker=="D8S1179"]) == c("14","15")))
  expect_true(all(unique(res$Allele[res$Marker=="FGA"]) == c("20","23")))
  expect_true(all(unique(res$Allele[res$Marker=="D2S441"]) == c("10","14")))
  expect_true(all(unique(res$Allele[res$Marker=="D12S391"]) == c("18","23")))
  expect_true(all(unique(res$Allele[res$Marker=="D19S433"]) == c("13","14")))
  expect_true(all(unique(res$Allele[res$Marker=="SE33"]) == c("15","16")))
  
  # TEST 08 -------------------------------------------------------------------
  # Test that filtering of 'dirty' samples without Sample.Name work.
  
  # Analyse dataframe.
  res <- filterProfile(data=set8, ref=ref8,
                       add.missing.loci=FALSE,
                       keep.na=FALSE,
                       ignore.case=FALSE,
                       debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Allele)))
  expect_false(any(is.na(res$Height)))
  
  # Check result.
  expect_that(ncol(res), equals(4))
  expect_that(nrow(res), equals(270))

  # Check that 'Ladder' is in result.
  expect_true("Ladder" %in% res$Sample.Name)
  
  # Check success of filtering.
  expect_true(all(unique(res$Allele[res$Marker=="AMEL"]) == c("X","Y")))
  expect_true(all(unique(res$Allele[res$Marker=="D3S1358"]) == c("17","18")))
  expect_true(all(unique(res$Allele[res$Marker=="TH01"]) == c("6","9.3")))
  expect_true(all(unique(res$Allele[res$Marker=="D21S11"]) == c("29","31.2")))
  expect_true(all(unique(res$Allele[res$Marker=="D18S51"]) == c("16","18")))
  expect_true(all(unique(res$Allele[res$Marker=="D10S1248"]) == c("13","15")))
  expect_true(all(unique(res$Allele[res$Marker=="D1S1656"]) == c("12","13")))
  expect_true(all(unique(res$Allele[res$Marker=="D2S1338"]) == c("22","25")))
  expect_true(all(unique(res$Allele[res$Marker=="D16S539"]) == c("9","13")))
  expect_true(all(unique(res$Allele[res$Marker=="D22S1045"]) == c("16")))
  expect_true(all(unique(res$Allele[res$Marker=="vWA"]) == c("16","19")))
  expect_true(all(unique(res$Allele[res$Marker=="D8S1179"]) == c("14","15")))
  expect_true(all(unique(res$Allele[res$Marker=="FGA"]) == c("20","23")))
  expect_true(all(unique(res$Allele[res$Marker=="D2S441"]) == c("10","14")))
  expect_true(all(unique(res$Allele[res$Marker=="D12S391"]) == c("18","23")))
  expect_true(all(unique(res$Allele[res$Marker=="D19S433"]) == c("13","14")))
  expect_true(all(unique(res$Allele[res$Marker=="SE33"]) == c("15","16")))
  
  # TEST 09 -------------------------------------------------------------------
  # Test that invert filtering of 'dirty' samples work.
  
  # Analyse dataframe.
  res <- filterProfile(data=set6, ref=ref6,
                       add.missing.loci=TRUE,
                       keep.na=TRUE,
                       ignore.case=FALSE,
                       invert=TRUE,
                       debug=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_true(any(is.na(res$Allele)))
  expect_true(any(is.na(res$Height)))
  
  # Check result.
  expect_that(ncol(res), equals(4))
  expect_that(nrow(res), equals(240))
  
  # Check success of filtering.
  expect_false(all(unique(res$Allele[res$Marker=="AMEL"]) %in% c("X","Y")))
  expect_false(all(unique(res$Allele[res$Marker=="D3S1358"]) %in% c("17","18")))
  expect_false(all(unique(res$Allele[res$Marker=="TH01"]) %in% c("6","9.3")))
  expect_false(all(unique(res$Allele[res$Marker=="D21S11"]) %in% c("29","31.2")))
  expect_false(all(unique(res$Allele[res$Marker=="D18S51"]) %in% c("16","18")))
  expect_false(all(unique(res$Allele[res$Marker=="D10S1248"]) %in% c("13","15")))
  expect_false(all(unique(res$Allele[res$Marker=="D1S1656"]) %in% c("12","13")))
  expect_false(all(unique(res$Allele[res$Marker=="D2S1338"]) %in% c("22","25")))
  expect_false(all(unique(res$Allele[res$Marker=="D16S539"]) %in% c("9","13")))
  expect_false(all(unique(res$Allele[res$Marker=="D22S1045"]) %in% c("16")))
  expect_false(all(unique(res$Allele[res$Marker=="vWA"]) %in% c("16","19")))
  expect_false(all(unique(res$Allele[res$Marker=="D8S1179"]) %in% c("14","15")))
  expect_false(all(unique(res$Allele[res$Marker=="FGA"]) %in% c("20","23")))
  expect_false(all(unique(res$Allele[res$Marker=="D2S441"]) %in% c("10","14")))
  expect_false(all(unique(res$Allele[res$Marker=="D12S391"]) %in% c("18","23")))
  expect_false(all(unique(res$Allele[res$Marker=="D19S433"]) %in% c("13","14")))
  expect_false(all(unique(res$Allele[res$Marker=="SE33"]) %in% c("15","16")))

  # TEST 10 -------------------------------------------------------------------
  # Test that filtering generates incorrect profile in 'fast' mode.
  
  # Analyse dataframe.
  res <- filterProfile(data = set9, ref = ref9,
                       add.missing.loci = FALSE,
                       keep.na = FALSE,
                       ignore.case = FALSE,
                       exact = FALSE,
                       debug = FALSE)
  
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))

  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Allele)))
  expect_false(any(is.na(res$Height)))

  # Check result.
  expect_that(ncol(res), equals(4))
  expect_that(nrow(res), equals(71))
  
  # Check if filtered profile equals reference.
  expect_false(is.logical(all.equal(res[res$Sample.Name == "NOPC15", ],
                                     ref9[ref9$Sample.Name == "NOPC15", ],
                                     check.attributes=FALSE)))
  
  # TEST 11 -------------------------------------------------------------------
  # Test that filtering generates correct profile in 'fast' mode.
  
  # Analyse dataframe.
  res <- filterProfile(data = set9, ref = ref9,
                       add.missing.loci = FALSE,
                       keep.na = FALSE,
                       ignore.case = FALSE,
                       exact = TRUE,
                       debug = FALSE)
  
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Allele)))
  expect_false(any(is.na(res$Height)))
  
  # Check result.
  expect_that(ncol(res), equals(4))
  expect_that(nrow(res), equals(66))
  
  # Check if filtered profile equals reference.
  expect_true(is.logical(all.equal(res[res$Sample.Name == "NOPC15", ],
                                   ref9[ref9$Sample.Name == "NOPC15", ],
                                   check.attributes=FALSE)))

  # TEST 12 -------------------------------------------------------------------
  # Test that filtering generates incorrect profile in 'slow' mode.
  
  # Analyse dataframe.
  res <- filterProfile(data = set9, ref = ref9,
                       add.missing.loci = TRUE,
                       keep.na = TRUE,
                       ignore.case = FALSE,
                       exact = FALSE,
                       debug = FALSE)
  
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Allele)))
  expect_false(any(is.na(res$Height)))
  
  # Check result.
  expect_that(ncol(res), equals(4))
  expect_that(nrow(res), equals(105))
  
  # Check if filtered profile equals reference.
  expect_false(is.logical(all.equal(res[res$Sample.Name == "NOPC15", ],
                         ref9[ref9$Sample.Name == "NOPC15", ],
                         check.attributes = FALSE)))

  # TEST 12 -------------------------------------------------------------------
  # Test that filtering generates correct profile in 'slow' mode.
  
  # Analyse dataframe.
  res <- filterProfile(data = set9, ref = ref9,
                       add.missing.loci = TRUE,
                       keep.na = TRUE,
                       ignore.case = FALSE,
                       exact = TRUE,
                       debug = FALSE)
  
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(any(grepl("Sample.Name", names(res))))
  expect_true(any(grepl("Marker", names(res))))
  expect_true(any(grepl("Allele", names(res))))
  expect_true(any(grepl("Height", names(res))))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Allele)))
  expect_false(any(is.na(res$Height)))
  
  # Check result.
  expect_that(ncol(res), equals(4))
  expect_that(nrow(res), equals(66))
  
  # Check if filtered profile equals reference.
  expect_true(is.logical(all.equal(res[res$Sample.Name == "NOPC15", ],
                                   ref9[ref9$Sample.Name == "NOPC15", ],
                                   check.attributes = FALSE)))

})