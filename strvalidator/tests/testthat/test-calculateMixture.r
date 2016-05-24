context("calculateMixture")

################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG
# 31.07.2014: First version.
# 
# library(testthat)
# test_dir("inst/tests/")
# test_file("tests/testthat/test-calculateMixture.r")
# test_dir("tests/testthat")

test_that("calculateMixture", {

  # Create reference profiles and mixtures.
  # Created using: dump("major", file="")
  major <-
    structure(list(Sample.Name = c("major", "major", "major", "major", 
                                   "major", "major", "major", "major",
                                   "major", "major", "major", "major",
                                   "major", "major", "major", "major",
                                   "major", "major", "major", "major",
                                   "major", "major", "major", "major",
                                   "major", "major", "major"),
                   Marker = c("AMEL", "D3S1358", "D3S1358", "TH01",
                              "D21S11", "D21S11", "D18S51", "D18S51",
                              "D10S1248", "D10S1248", "D1S1656", "D1S1656",
                              "D2S1338", "D16S539", "D16S539", "D22S1045",
                              "vWA", "D8S1179", "D8S1179", "FGA", "D2S441",
                              "D2S441", "D12S391", "D12S391", "D19S433",
                              "SE33", "SE33"),
                   Allele = c("X", "14", "15", "8", "25.2", "28", "12",
                              "15", "14", "15", "16", "16.1", "20", "11",
                              "12", "11", "17", "11", "13", "21", "11.3",
                              "14", "21", "23", "13", "18", "20")),
              .Names = c("Sample.Name", "Marker", "Allele"),
              row.names = c(NA, 27L), class = "data.frame")
  
  minor <-
    structure(list(Sample.Name = c("minor", "minor", "minor", "minor", 
                                   "minor", "minor", "minor", "minor",
                                   "minor", "minor", "minor", "minor",
                                   "minor", "minor", "minor", "minor",
                                   "minor", "minor", "minor", "minor",
                                   "minor", "minor", "minor", "minor",
                                   "minor", "minor", "minor", "minor",
                                   "minor"),
                   Marker = c("AMEL", "AMEL", "D3S1358", "D3S1358",
                              "TH01", "TH01", "D21S11", "D21S11",
                              "D18S51", "D18S51", "D10S1248", "D1S1656",
                              "D1S1656", "D2S1338", "D2S1338", "D16S539",
                              "D16S539", "D22S1045", "D22S1045", "vWA",
                              "D8S1179", "D8S1179", "FGA", "D2S441",
                              "D12S391", "D12S391", "D19S433", "SE33", "SE33"),
                   Allele = c("X", "Y", "16", "18", "8", "9.3",
                              "25.2", "28", "15", "17", "15", "13",
                              "16", "19", "20", "9", "10", "16", "17",
                              "19", "8", "10", "25", "14", "19", "22",
                              "13", "17", "28.2")),
              .Names = c("Sample.Name", "Marker", "Allele"),
              class = "data.frame",
              row.names = c(NA, 29L))
  
  mixture <- 
    structure(list(Sample.Name = c("major_minor_1", "major_minor_1",
                                   "major_minor_1", "major_minor_1",
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1",
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_1", "major_minor_1",
                                   "major_minor_1", "major_minor_1", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2", "major_minor_2", 
                                   "major_minor_2"),
                   Marker = c("AMEL", "AMEL", "D3S1358", "D3S1358", 
                              "D3S1358", "D3S1358", "TH01", "TH01",
                              "D21S11", "D21S11", "D18S51", "D18S51",
                              "D18S51", "D10S1248", "D10S1248", "D1S1656",
                              "D1S1656", "D1S1656", "D2S1338", "D2S1338",
                              "D16S539", "D16S539", "D16S539", "D16S539",
                              "D22S1045", "D22S1045", "D22S1045", "vWA",
                              "vWA", "D8S1179", "D8S1179", "D8S1179",
                              "D8S1179", "FGA", "FGA", "D2S441",
                              "D2S441", "D12S391", "D12S391", "D12S391",
                              "D12S391", "D19S433", "SE33", "SE33",
                              "SE33", "SE33", "AMEL", "AMEL",
                              "D3S1358", "D3S1358", "D3S1358", "D3S1358",
                              "TH01", "D21S11", "D21S11", "D21S11",
                              "D18S51", "D18S51", "D10S1248", "D10S1248",
                              "D1S1656", "D1S1656", "D1S1656","D1S1656",
                              "D2S1338", "D2S1338", "D16S539", "D16S539",
                              "D16S539", "D16S539", "D22S1045", "D22S1045",
                              "D22S1045", "vWA", "vWA", "D8S1179",
                              "D8S1179", "D8S1179", "FGA", "FGA",
                              "D2S441", "D2S441", "D12S391", "D12S391",
                              "D12S391", "D12S391", "D19S433", "D19S433",
                              "SE33", "SE33", "SE33"),
                   Allele = c("X", "Y", "14", "15", "16", "18", "8", "9.3",
                              "25.2", "28", "12", "15", "17", "14", "15",
                              "13", "16", "16.1", "19", "20", "9", "10",
                              "11", "12", "11", "16", "17", "17", "19", "8",
                              "10", "11", "13", "21", "25", "11.3", "14",
                              "19", "21", "22", "23", "13", "17", "18", "20",
                              "28.2", "OL", "X", "14", "15", "16", "18", "8",
                              "25.2", "28", "29", "12", "15", "14", "15",
                              "13", "16", "16.1", "OL", "19", "20", "9", 
                              "10", "11", "12", "11", "16", "17", "17", "19",
                              "10", "11", "13", "21", "22", "11.3", "14",
                              "19", "21", "22", "23", "13", "14", "17", "18",
                              "20"),
                   Height = c("7533", "1503", "4500", "4200", "1000", "1100",
                              "9300", "1200", "4100", "3800", "1700", "3300",
                              "1500", "3000", "7000", "1300", "5600", "3800",
                              "1200", "8500", "1600", "1600", "3500", "3400",
                              "8600", "1200", "900", "8200", "1900", "1200",
                              "1400", "5000", "4600", "7800", "2100", "6000",
                              "8000", "900", "3400", "1200", "3800", "12000",
                              "1000", "4000", "4100", "1100", "215", "7533",
                              "4500", "4200", "1000", "1100", "9300", "4100",
                              "3800", "340", "1700", "3300", "3000", "7000", 
                              "1300", "5600", "3800", "1662", "1200", "8500",
                              "1600", "1600", "3500", "3400", "8600", "1200",
                              "900", "8200", "1900", "1400", "5000", "4600",
                              "7800", "250", "6000", "8000", "900", "3400", 
                              "1200", "3800", "12000", "100", "1000", "4000",
                              "4100"),
                   Style = c("AA:AB", "AA:AB", "AB:CD", "AB:CD", "AB:CD",
                               "AB:CD", "AA:AB", "AA:AB", "AB:AB", "AB:AB",
                               "AB:AC", "AB:AC", "AB:AC", "AB:AA", "AB:AA",
                               "AB:AC", "AB:AC", "AB:AC", "AB:AA", "AB:AA",
                               "AB:CD", "AB:CD", "AB:CD", "AB:CD", "AA:BC",
                               "AA:BC", "AA:BC", "AA:BB", "AA:BB", "AB:CD",
                               "AB:CD", "AB:CD", "AB:CD", "AA:BB", "AA:BB",
                               "AA:AB", "AA:AB", "AB:CD", "AB:CD", "AB:CD",
                               "AB:CD", "AA:AA", "AB:CD", "AB:CD", "AB:CD",
                               "AB:CD", "OL", "AA:A!B", "AB:CD", "AB:CD", 
                               "AB:CD", "AB:CD", "AA:A!B", "AB:AB", "AB:AB",
                               "Dropin", "AB:A!C", "AB:A!C", "AB:AA",
                               "AB:AA", "AB:AC", "AB:AC", "AB:AC", "OL", 
                               "AB:AA", "AB:AA", "AB:CD", "AB:CD", "AB:CD",
                               "AB:CD", "AA:BC", "AA:BC", "AA:BC", "AA:BB",
                               "AA:BB", "AB:C!D", "AB:C!D", "AB:C!D", 
                               "AA:!B!B", "Dropin", "AA:AB", "AA:AB", "AB:CD",
                               "AB:CD", "AB:CD", "AB:CD", "AA:AA", "Dropin",
                               "AB:C!D", "AB:C!D", "AB:C!D"),
                   Major.minor = c("MAJOR/MAJOR/minor",
                                   "minor", "MAJOR", "MAJOR", "minor",
                                   "minor", "MAJOR/MAJOR/minor", "minor",
                                   "MAJOR/minor", "MAJOR/minor", "MAJOR",
                                   "MAJOR/minor", "minor", "MAJOR",
                                   "MAJOR/minor/minor", "minor", 
                                   "MAJOR/minor", "MAJOR", "minor", 
                                   "MAJOR/MAJOR/minor", "minor", "minor",
                                   "MAJOR", "MAJOR", "MAJOR/MAJOR", "minor",
                                   "minor", "MAJOR/MAJOR", "minor/minor", 
                                   "minor", "minor", "MAJOR", "MAJOR", 
                                   "MAJOR/MAJOR", "minor/minor", "MAJOR",
                                   "MAJOR/minor/minor", "minor", "MAJOR", 
                                   "minor", "MAJOR", "MAJOR/MAJOR/minor/minor",
                                   "minor", "MAJOR", "MAJOR", "minor", NA, 
                                   "MAJOR/MAJOR/minor", "MAJOR", "MAJOR", 
                                   "minor", "minor", "MAJOR/MAJOR/minor", 
                                   "MAJOR/minor", "MAJOR/minor", NA, "MAJOR",
                                   "MAJOR/minor", "MAJOR", "MAJOR/minor/minor",
                                   "minor", "MAJOR/minor", "MAJOR", NA, 
                                   "minor", "MAJOR/MAJOR/minor", "minor", 
                                   "minor", "MAJOR", "MAJOR", "MAJOR/MAJOR",
                                   "minor", "minor", "MAJOR/MAJOR", 
                                   "minor/minor", "minor", "MAJOR", "MAJOR",
                                   "MAJOR/MAJOR", NA, "MAJOR", 
                                   "MAJOR/minor/minor", "minor", "MAJOR",
                                   "minor", "MAJOR", "MAJOR/MAJOR/minor/minor",
                                   NA, "minor", "MAJOR", "MAJOR")),
              .Names = c("Sample.Name", "Marker", "Allele", "Height",
                         "Style", "Major.minor"),
              class = "data.frame", row.names = c(NA, 91L))
  
  # TEST 01 -------------------------------------------------------------------
  # Remove off-ladder alleles, ignore drop-outs.
  
  # Analyse dataframe.
  res <- calculateMixture(data=mixture, ref1=major, ref2=minor,
                          ol.rm=TRUE, ignore.dropout=TRUE)

  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_false(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$Style))
  expect_false(is.null(res$Mx))
  expect_false(is.null(res$Average))
  expect_false(is.null(res$Difference))
  expect_false(is.null(res$Observed))
  expect_false(is.null(res$Expected))
  expect_false(is.null(res$Profile))
  expect_false(is.null(res$Dropin))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Style)))
  expect_true(any(is.na(res$Mx)))
  expect_false(any(is.na(res$Average)))
  expect_true(any(is.na(res$Difference)))
  expect_false(any(is.na(res$Observed)))
  expect_false(any(is.na(res$Expected)))
  expect_false(any(is.na(res$Profile)))
  expect_false(any(is.na(res$Dropin)))
  
  # Check result: sample name.
  expect_that(unique(res$Sample.Name),
              equals(c("major_minor_1","major_minor_2")))
  
  # Check result: average Mx.
  expect_that(round(unique(res$Average),4),
              equals(c(0.2552,0.1603)))
  
  # Check result: profile.
  expect_that(round(unique(res$Profile),2),
              equals(c(100.00,68.42)))
  
  # Check result: Mx.
  expect_that(res$Mx[1], equals((2*1503)/(1503+7533)))
  expect_that(res$Mx[2], equals((1000+1100)/(1000+1100+4500+4200)))
  expect_that(res$Mx[3], equals((2*1200)/(1200+9300)))
  expect_that(res$Mx[4], equals(as.numeric(NA)))
  expect_that(res$Mx[5], equals((1500)/(1500+1700)))
  expect_that(res$Mx[6], equals((7000-3000)/(7000+3000)))
  expect_that(res$Mx[7], equals(1300/(1300+3800)))
  expect_that(res$Mx[8], equals((2*1200)/(1200+8500)))
  expect_that(res$Mx[9], equals((1600+1600)/(1600+1600+3500+3400)))
  expect_that(res$Mx[10], equals((1200+900)/(1200+900+8600)))
  expect_that(res$Mx[11], equals(1900/(1900+8200)))
  expect_that(res$Mx[12], equals((1200+1400)/(1200+1400+5000+4600)))
  expect_that(res$Mx[13], equals(2100/(2100+7800)))
  expect_that(res$Mx[14], equals((8000-6000)/(8000+6000)))
  expect_that(res$Mx[15], equals((900+1200)/(900+1200+3400+3800)))
  expect_that(res$Mx[16], equals(as.numeric(NA)))
  expect_that(res$Mx[17], equals((1000+1100)/(1000+1100+4000+4100)))
  
  expect_that(res$Mx[18], equals(0/7533))
  expect_that(res$Mx[19], equals((1000+1100)/(1000+1100+4500+4200)))
  expect_that(res$Mx[20], equals(0/9300))
  expect_that(res$Mx[21], equals(as.numeric(NA)))
  expect_that(res$Mx[22], equals(0/1700))
  expect_that(res$Mx[23], equals((7000-3000)/(7000+3000)))
  expect_that(res$Mx[24], equals(1300/(1300+3800)))
  expect_that(res$Mx[25], equals((2*1200)/(1200+8500)))
  expect_that(res$Mx[26], equals((1600+1600)/(1600+1600+3500+3400)))
  expect_that(res$Mx[27], equals((1200+900)/(1200+900+8600)))
  expect_that(res$Mx[28], equals(1900/(1900+8200)))
  expect_that(res$Mx[29], equals(1400/(1400+5000+4600)))
  expect_that(res$Mx[30], equals(0/7800))
  expect_that(res$Mx[31], equals((8000-6000)/(8000+6000)))
  expect_that(res$Mx[32], equals((900+1200)/(900+1200+3400+3800)))
  expect_that(res$Mx[33], equals(as.numeric(NA)))
  expect_that(res$Mx[34], equals(1000/(1000+4000+4100)))

  # Check result: Style.
  expect_that(res$Style[1], equals("AB:AA"))
  expect_that(res$Style[2], equals("AB:CD"))
  expect_that(res$Style[3], equals("AB:AA"))
  expect_that(res$Style[4], equals("AB:AB"))
  expect_that(res$Style[5], equals("AB:AC"))
  expect_that(res$Style[6], equals("AA:AB"))
  expect_that(res$Style[7], equals("AB:AC"))
  expect_that(res$Style[8], equals("AB:AA"))
  expect_that(res$Style[9], equals("AB:CD"))
  expect_that(res$Style[10], equals("AB:CC"))
  expect_that(res$Style[11], equals("AA:BB"))
  expect_that(res$Style[12], equals("AB:CD"))
  expect_that(res$Style[13], equals("AA:BB"))
  expect_that(res$Style[14], equals("AA:AB"))
  expect_that(res$Style[15], equals("AB:CD"))
  expect_that(res$Style[16], equals("AA:AA"))
  expect_that(res$Style[17], equals("AB:CD"))
  
  expect_that(res$Style[18], equals("AB:AA"))
  expect_that(res$Style[19], equals("AB:CD"))
  expect_that(res$Style[20], equals("AB:AA"))
  expect_that(res$Style[21], equals("AB:AB"))
  expect_that(res$Style[22], equals("AB:AC"))
  expect_that(res$Style[23], equals("AA:AB"))
  expect_that(res$Style[24], equals("AB:AC"))
  expect_that(res$Style[25], equals("AB:AA"))
  expect_that(res$Style[26], equals("AB:CD"))
  expect_that(res$Style[27], equals("AB:CC"))
  expect_that(res$Style[28], equals("AA:BB"))
  expect_that(res$Style[29], equals("AB:CD"))
  expect_that(res$Style[30], equals("AA:BB"))
  expect_that(res$Style[31], equals("AA:AB"))
  expect_that(res$Style[32], equals("AB:CD"))
  expect_that(res$Style[33], equals("AA:AA"))
  expect_that(res$Style[34], equals("AB:CD"))
  
  # Check result: Observed.
  expect_that(res$Observed[res$Sample.Name=="major_minor_1"],
              equals(c(1,2,1,0,1,0,1,1,2,2,1,2,1,0,2,0,2)))
  expect_that(res$Observed[res$Sample.Name=="major_minor_2"],
              equals(c(0,2,0,0,0,0,1,1,2,2,1,1,0,0,2,0,1)))

  # Check result: Expected.
  expect_that(res$Expected[res$Sample.Name=="major_minor_1"],
              equals(c(1,2,1,0,1,0,1,1,2,2,1,2,1,0,2,0,2)))
  expect_that(res$Expected[res$Sample.Name=="major_minor_2"],
              equals(c(1,2,1,0,1,0,1,1,2,2,1,2,1,0,2,0,2)))
  
  # Check result: Dropin.
  expect_that(res$Dropin[res$Sample.Name=="major_minor_1"],
              equals(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)))
  expect_that(res$Dropin[res$Sample.Name=="major_minor_2"],
              equals(c(0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0)))
  
  
  # TEST 02 -------------------------------------------------------------------
  # Count OL as dropin, ignore drop-outs.
  
  # Analyse dataframe.
  res <- calculateMixture(data=mixture, ref1=major, ref2=minor,
                          ol.rm=FALSE, ignore.dropout=TRUE)

  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_false(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$Style))
  expect_false(is.null(res$Mx))
  expect_false(is.null(res$Average))
  expect_false(is.null(res$Difference))
  expect_false(is.null(res$Observed))
  expect_false(is.null(res$Expected))
  expect_false(is.null(res$Profile))
  expect_false(is.null(res$Dropin))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Style)))
  expect_true(any(is.na(res$Mx)))
  expect_false(any(is.na(res$Average)))
  expect_true(any(is.na(res$Difference)))
  expect_false(any(is.na(res$Observed)))
  expect_false(any(is.na(res$Expected)))
  expect_false(any(is.na(res$Profile)))
  expect_false(any(is.na(res$Dropin)))
  
  # Check result: sample name.
  expect_that(unique(res$Sample.Name),
              equals(c("major_minor_1","major_minor_2")))
  
  # Check result: average Mx.
  expect_that(round(unique(res$Average),4),
              equals(c(0.2552,0.1603)))
  
  # Check result: profile.
  expect_that(round(unique(res$Profile),2),
              equals(c(100.00,68.42)))
  
  # Check result: Mx.
  expect_that(res$Mx[1], equals((2*1503)/(1503+7533)))
  expect_that(res$Mx[2], equals((1000+1100)/(1000+1100+4500+4200)))
  expect_that(res$Mx[3], equals((2*1200)/(1200+9300)))
  expect_that(res$Mx[4], equals(as.numeric(NA)))
  expect_that(res$Mx[5], equals((1500)/(1500+1700)))
  expect_that(res$Mx[6], equals((7000-3000)/(7000+3000)))
  expect_that(res$Mx[7], equals(1300/(1300+3800)))
  expect_that(res$Mx[8], equals((2*1200)/(1200+8500)))
  expect_that(res$Mx[9], equals((1600+1600)/(1600+1600+3500+3400)))
  expect_that(res$Mx[10], equals((1200+900)/(1200+900+8600)))
  expect_that(res$Mx[11], equals(1900/(1900+8200)))
  expect_that(res$Mx[12], equals((1200+1400)/(1200+1400+5000+4600)))
  expect_that(res$Mx[13], equals(2100/(2100+7800)))
  expect_that(res$Mx[14], equals((8000-6000)/(8000+6000)))
  expect_that(res$Mx[15], equals((900+1200)/(900+1200+3400+3800)))
  expect_that(res$Mx[16], equals(as.numeric(NA)))
  expect_that(res$Mx[17], equals((1000+1100)/(1000+1100+4000+4100)))
  
  expect_that(res$Mx[18], equals(0/7533))
  expect_that(res$Mx[19], equals((1000+1100)/(1000+1100+4500+4200)))
  expect_that(res$Mx[20], equals(0/9300))
  expect_that(res$Mx[21], equals(as.numeric(NA)))
  expect_that(res$Mx[22], equals(0/1700))
  expect_that(res$Mx[23], equals((7000-3000)/(7000+3000)))
  expect_that(res$Mx[24], equals(1300/(1300+3800)))
  expect_that(res$Mx[25], equals((2*1200)/(1200+8500)))
  expect_that(res$Mx[26], equals((1600+1600)/(1600+1600+3500+3400)))
  expect_that(res$Mx[27], equals((1200+900)/(1200+900+8600)))
  expect_that(res$Mx[28], equals(1900/(1900+8200)))
  expect_that(res$Mx[29], equals(1400/(1400+5000+4600)))
  expect_that(res$Mx[30], equals(0/7800))
  expect_that(res$Mx[31], equals((8000-6000)/(8000+6000)))
  expect_that(res$Mx[32], equals((900+1200)/(900+1200+3400+3800)))
  expect_that(res$Mx[33], equals(as.numeric(NA)))
  expect_that(res$Mx[34], equals(1000/(1000+4000+4100)))

  # Check result: Style.
  expect_that(res$Style[1], equals("AB:AA"))
  expect_that(res$Style[2], equals("AB:CD"))
  expect_that(res$Style[3], equals("AB:AA"))
  expect_that(res$Style[4], equals("AB:AB"))
  expect_that(res$Style[5], equals("AB:AC"))
  expect_that(res$Style[6], equals("AA:AB"))
  expect_that(res$Style[7], equals("AB:AC"))
  expect_that(res$Style[8], equals("AB:AA"))
  expect_that(res$Style[9], equals("AB:CD"))
  expect_that(res$Style[10], equals("AB:CC"))
  expect_that(res$Style[11], equals("AA:BB"))
  expect_that(res$Style[12], equals("AB:CD"))
  expect_that(res$Style[13], equals("AA:BB"))
  expect_that(res$Style[14], equals("AA:AB"))
  expect_that(res$Style[15], equals("AB:CD"))
  expect_that(res$Style[16], equals("AA:AA"))
  expect_that(res$Style[17], equals("AB:CD"))
  
  expect_that(res$Style[18], equals("AB:AA"))
  expect_that(res$Style[19], equals("AB:CD"))
  expect_that(res$Style[20], equals("AB:AA"))
  expect_that(res$Style[21], equals("AB:AB"))
  expect_that(res$Style[22], equals("AB:AC"))
  expect_that(res$Style[23], equals("AA:AB"))
  expect_that(res$Style[24], equals("AB:AC"))
  expect_that(res$Style[25], equals("AB:AA"))
  expect_that(res$Style[26], equals("AB:CD"))
  expect_that(res$Style[27], equals("AB:CC"))
  expect_that(res$Style[28], equals("AA:BB"))
  expect_that(res$Style[29], equals("AB:CD"))
  expect_that(res$Style[30], equals("AA:BB"))
  expect_that(res$Style[31], equals("AA:AB"))
  expect_that(res$Style[32], equals("AB:CD"))
  expect_that(res$Style[33], equals("AA:AA"))
  expect_that(res$Style[34], equals("AB:CD"))
  
  # Check result: Observed.
  expect_that(res$Observed[res$Sample.Name=="major_minor_1"],
              equals(c(1,2,1,0,1,0,1,1,2,2,1,2,1,0,2,0,2)))
  expect_that(res$Observed[res$Sample.Name=="major_minor_2"],
              equals(c(0,2,0,0,0,0,1,1,2,2,1,1,0,0,2,0,1)))
  
  # Check result: Expected.
  expect_that(res$Expected[res$Sample.Name=="major_minor_1"],
              equals(c(1,2,1,0,1,0,1,1,2,2,1,2,1,0,2,0,2)))
  expect_that(res$Expected[res$Sample.Name=="major_minor_2"],
              equals(c(1,2,1,0,1,0,1,1,2,2,1,2,1,0,2,0,2)))
  
  # Check result: Dropin.
  expect_that(res$Dropin[res$Sample.Name=="major_minor_1"],
              equals(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)))
  expect_that(res$Dropin[res$Sample.Name=="major_minor_2"],
              equals(c(1,0,0,1,0,0,1,0,0,0,0,0,1,0,0,1,0)))
  
  # TEST 03 -------------------------------------------------------------------
  # Remove off-ladder alleles, do not use drop-outs.
  
  # Analyse dataframe.
  res <- calculateMixture(data=mixture, ref1=major, ref2=minor,
                          ol.rm=TRUE, ignore.dropout=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_false(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$Style))
  expect_false(is.null(res$Mx))
  expect_false(is.null(res$Average))
  expect_false(is.null(res$Difference))
  expect_false(is.null(res$Observed))
  expect_false(is.null(res$Expected))
  expect_false(is.null(res$Profile))
  expect_false(is.null(res$Dropin))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Style)))
  expect_true(any(is.na(res$Mx)))
  expect_false(any(is.na(res$Average)))
  expect_true(any(is.na(res$Difference)))
  expect_false(any(is.na(res$Observed)))
  expect_false(any(is.na(res$Expected)))
  expect_false(any(is.na(res$Profile)))
  expect_false(any(is.na(res$Dropin)))
  
  # Check result: sample name.
  expect_that(unique(res$Sample.Name),
              equals(c("major_minor_1","major_minor_2")))
  
  # Check result: average Mx.
  expect_that(round(unique(res$Average),4),
              equals(c(0.2552,0.2407)))
  
  # Check result: profile.
  expect_that(round(unique(res$Profile),2),
              equals(c(100.00,68.42)))
  
  # Check result: Mx.
  expect_that(res$Mx[1], equals((2*1503)/(1503+7533)))
  expect_that(res$Mx[2], equals((1000+1100)/(1000+1100+4500+4200)))
  expect_that(res$Mx[3], equals((2*1200)/(1200+9300)))
  expect_that(res$Mx[4], equals(as.numeric(NA)))
  expect_that(res$Mx[5], equals((1500)/(1500+1700)))
  expect_that(res$Mx[6], equals((7000-3000)/(7000+3000)))
  expect_that(res$Mx[7], equals(1300/(1300+3800)))
  expect_that(res$Mx[8], equals((2*1200)/(1200+8500)))
  expect_that(res$Mx[9], equals((1600+1600)/(1600+1600+3500+3400)))
  expect_that(res$Mx[10], equals((1200+900)/(1200+900+8600)))
  expect_that(res$Mx[11], equals(1900/(1900+8200)))
  expect_that(res$Mx[12], equals((1200+1400)/(1200+1400+5000+4600)))
  expect_that(res$Mx[13], equals(2100/(2100+7800)))
  expect_that(res$Mx[14], equals((8000-6000)/(8000+6000)))
  expect_that(res$Mx[15], equals((900+1200)/(900+1200+3400+3800)))
  expect_that(res$Mx[16], equals(as.numeric(NA)))
  expect_that(res$Mx[17], equals((1000+1100)/(1000+1100+4000+4100)))
  
  expect_that(res$Mx[18], equals(as.numeric(NA)))
  expect_that(res$Mx[19], equals((1000+1100)/(1000+1100+4500+4200)))
  expect_that(res$Mx[20], equals(as.numeric(NA)))
  expect_that(res$Mx[21], equals(as.numeric(NA)))
  expect_that(res$Mx[22], equals(as.numeric(NA)))
  expect_that(res$Mx[23], equals((7000-3000)/(7000+3000)))
  expect_that(res$Mx[24], equals(1300/(1300+3800)))
  expect_that(res$Mx[25], equals((2*1200)/(1200+8500)))
  expect_that(res$Mx[26], equals((1600+1600)/(1600+1600+3500+3400)))
  expect_that(res$Mx[27], equals((1200+900)/(1200+900+8600)))
  expect_that(res$Mx[28], equals(1900/(1900+8200)))
  expect_that(res$Mx[29], equals(as.numeric(NA)))
  expect_that(res$Mx[30], equals(as.numeric(NA)))
  expect_that(res$Mx[31], equals((8000-6000)/(8000+6000)))
  expect_that(res$Mx[32], equals((900+1200)/(900+1200+3400+3800)))
  expect_that(res$Mx[33], equals(as.numeric(NA)))
  expect_that(res$Mx[34], equals(as.numeric(NA)))
  
  # Check result: Style.
  expect_that(res$Style[1], equals("AB:AA"))
  expect_that(res$Style[2], equals("AB:CD"))
  expect_that(res$Style[3], equals("AB:AA"))
  expect_that(res$Style[4], equals("AB:AB"))
  expect_that(res$Style[5], equals("AB:AC"))
  expect_that(res$Style[6], equals("AA:AB"))
  expect_that(res$Style[7], equals("AB:AC"))
  expect_that(res$Style[8], equals("AB:AA"))
  expect_that(res$Style[9], equals("AB:CD"))
  expect_that(res$Style[10], equals("AB:CC"))
  expect_that(res$Style[11], equals("AA:BB"))
  expect_that(res$Style[12], equals("AB:CD"))
  expect_that(res$Style[13], equals("AA:BB"))
  expect_that(res$Style[14], equals("AA:AB"))
  expect_that(res$Style[15], equals("AB:CD"))
  expect_that(res$Style[16], equals("AA:AA"))
  expect_that(res$Style[17], equals("AB:CD"))
  
  expect_that(res$Style[18], equals("Dropout"))
  expect_that(res$Style[19], equals("AB:CD"))
  expect_that(res$Style[20], equals("Dropout"))
  expect_that(res$Style[21], equals("AB:AB"))
  expect_that(res$Style[22], equals("Dropout"))
  expect_that(res$Style[23], equals("AA:AB"))
  expect_that(res$Style[24], equals("AB:AC"))
  expect_that(res$Style[25], equals("AB:AA"))
  expect_that(res$Style[26], equals("AB:CD"))
  expect_that(res$Style[27], equals("AB:CC"))
  expect_that(res$Style[28], equals("AA:BB"))
  expect_that(res$Style[29], equals("Dropout"))
  expect_that(res$Style[30], equals("Dropout"))
  expect_that(res$Style[31], equals("AA:AB"))
  expect_that(res$Style[32], equals("AB:CD"))
  expect_that(res$Style[33], equals("AA:AA"))
  expect_that(res$Style[34], equals("Dropout"))
  
  # Check result: Observed.
  expect_that(res$Observed[res$Sample.Name=="major_minor_1"],
              equals(c(1,2,1,0,1,0,1,1,2,2,1,2,1,0,2,0,2)))
  expect_that(res$Observed[res$Sample.Name=="major_minor_2"],
              equals(c(0,2,0,0,0,0,1,1,2,2,1,1,0,0,2,0,1)))
  
  # Check result: Expected.
  expect_that(res$Expected[res$Sample.Name=="major_minor_1"],
              equals(c(1,2,1,0,1,0,1,1,2,2,1,2,1,0,2,0,2)))
  expect_that(res$Expected[res$Sample.Name=="major_minor_2"],
              equals(c(1,2,1,0,1,0,1,1,2,2,1,2,1,0,2,0,2)))
  
  # Check result: Dropin.
  expect_that(res$Dropin[res$Sample.Name=="major_minor_1"],
              equals(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)))
  expect_that(res$Dropin[res$Sample.Name=="major_minor_2"],
              equals(c(0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0)))
  
  
  # TEST 04 -------------------------------------------------------------------
  # Count OL as dropin, do not use drop-outs.
  
  # Analyse dataframe.
  res <- calculateMixture(data=mixture, ref1=major, ref2=minor,
                          ol.rm=FALSE, ignore.dropout=FALSE)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_false(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$Style))
  expect_false(is.null(res$Mx))
  expect_false(is.null(res$Average))
  expect_false(is.null(res$Difference))
  expect_false(is.null(res$Observed))
  expect_false(is.null(res$Expected))
  expect_false(is.null(res$Profile))
  expect_false(is.null(res$Dropin))
  
  # Check for NA's.
  expect_false(any(is.na(res$Sample.Name)))
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Style)))
  expect_true(any(is.na(res$Mx)))
  expect_false(any(is.na(res$Average)))
  expect_true(any(is.na(res$Difference)))
  expect_false(any(is.na(res$Observed)))
  expect_false(any(is.na(res$Expected)))
  expect_false(any(is.na(res$Profile)))
  expect_false(any(is.na(res$Dropin)))
  
  # Check result: sample name.
  expect_that(unique(res$Sample.Name),
              equals(c("major_minor_1","major_minor_2")))
  
  # Check result: average Mx.
  expect_that(round(unique(res$Average),4),
              equals(c(0.2552,0.2407)))
  
  # Check result: profile.
  expect_that(round(unique(res$Profile),2),
              equals(c(100.00,68.42)))
  
  # Check result: Mx.
  expect_that(res$Mx[1], equals((2*1503)/(1503+7533)))
  expect_that(res$Mx[2], equals((1000+1100)/(1000+1100+4500+4200)))
  expect_that(res$Mx[3], equals((2*1200)/(1200+9300)))
  expect_that(res$Mx[4], equals(as.numeric(NA)))
  expect_that(res$Mx[5], equals((1500)/(1500+1700)))
  expect_that(res$Mx[6], equals((7000-3000)/(7000+3000)))
  expect_that(res$Mx[7], equals(1300/(1300+3800)))
  expect_that(res$Mx[8], equals((2*1200)/(1200+8500)))
  expect_that(res$Mx[9], equals((1600+1600)/(1600+1600+3500+3400)))
  expect_that(res$Mx[10], equals((1200+900)/(1200+900+8600)))
  expect_that(res$Mx[11], equals(1900/(1900+8200)))
  expect_that(res$Mx[12], equals((1200+1400)/(1200+1400+5000+4600)))
  expect_that(res$Mx[13], equals(2100/(2100+7800)))
  expect_that(res$Mx[14], equals((8000-6000)/(8000+6000)))
  expect_that(res$Mx[15], equals((900+1200)/(900+1200+3400+3800)))
  expect_that(res$Mx[16], equals(as.numeric(NA)))
  expect_that(res$Mx[17], equals((1000+1100)/(1000+1100+4000+4100)))
  
  expect_that(res$Mx[18], equals(as.numeric(NA)))
  expect_that(res$Mx[19], equals((1000+1100)/(1000+1100+4500+4200)))
  expect_that(res$Mx[20], equals(as.numeric(NA)))
  expect_that(res$Mx[21], equals(as.numeric(NA)))
  expect_that(res$Mx[22], equals(as.numeric(NA)))
  expect_that(res$Mx[23], equals((7000-3000)/(7000+3000)))
  expect_that(res$Mx[24], equals(1300/(1300+3800)))
  expect_that(res$Mx[25], equals((2*1200)/(1200+8500)))
  expect_that(res$Mx[26], equals((1600+1600)/(1600+1600+3500+3400)))
  expect_that(res$Mx[27], equals((1200+900)/(1200+900+8600)))
  expect_that(res$Mx[28], equals(1900/(1900+8200)))
  expect_that(res$Mx[29], equals(as.numeric(NA)))
  expect_that(res$Mx[30], equals(as.numeric(NA)))
  expect_that(res$Mx[31], equals((8000-6000)/(8000+6000)))
  expect_that(res$Mx[32], equals((900+1200)/(900+1200+3400+3800)))
  expect_that(res$Mx[33], equals(as.numeric(NA)))
  expect_that(res$Mx[34], equals(as.numeric(NA)))
  
  # Check result: Style.
  expect_that(res$Style[1], equals("AB:AA"))
  expect_that(res$Style[2], equals("AB:CD"))
  expect_that(res$Style[3], equals("AB:AA"))
  expect_that(res$Style[4], equals("AB:AB"))
  expect_that(res$Style[5], equals("AB:AC"))
  expect_that(res$Style[6], equals("AA:AB"))
  expect_that(res$Style[7], equals("AB:AC"))
  expect_that(res$Style[8], equals("AB:AA"))
  expect_that(res$Style[9], equals("AB:CD"))
  expect_that(res$Style[10], equals("AB:CC"))
  expect_that(res$Style[11], equals("AA:BB"))
  expect_that(res$Style[12], equals("AB:CD"))
  expect_that(res$Style[13], equals("AA:BB"))
  expect_that(res$Style[14], equals("AA:AB"))
  expect_that(res$Style[15], equals("AB:CD"))
  expect_that(res$Style[16], equals("AA:AA"))
  expect_that(res$Style[17], equals("AB:CD"))
  
  expect_that(res$Style[18], equals("Dropout"))
  expect_that(res$Style[19], equals("AB:CD"))
  expect_that(res$Style[20], equals("Dropout"))
  expect_that(res$Style[21], equals("AB:AB"))
  expect_that(res$Style[22], equals("Dropout"))
  expect_that(res$Style[23], equals("AA:AB"))
  expect_that(res$Style[24], equals("AB:AC"))
  expect_that(res$Style[25], equals("AB:AA"))
  expect_that(res$Style[26], equals("AB:CD"))
  expect_that(res$Style[27], equals("AB:CC"))
  expect_that(res$Style[28], equals("AA:BB"))
  expect_that(res$Style[29], equals("Dropout"))
  expect_that(res$Style[30], equals("Dropout"))
  expect_that(res$Style[31], equals("AA:AB"))
  expect_that(res$Style[32], equals("AB:CD"))
  expect_that(res$Style[33], equals("AA:AA"))
  expect_that(res$Style[34], equals("Dropout"))
  
  # Check result: Observed.
  expect_that(res$Observed[res$Sample.Name=="major_minor_1"],
              equals(c(1,2,1,0,1,0,1,1,2,2,1,2,1,0,2,0,2)))
  expect_that(res$Observed[res$Sample.Name=="major_minor_2"],
              equals(c(0,2,0,0,0,0,1,1,2,2,1,1,0,0,2,0,1)))
  
  # Check result: Expected.
  expect_that(res$Expected[res$Sample.Name=="major_minor_1"],
              equals(c(1,2,1,0,1,0,1,1,2,2,1,2,1,0,2,0,2)))
  expect_that(res$Expected[res$Sample.Name=="major_minor_2"],
              equals(c(1,2,1,0,1,0,1,1,2,2,1,2,1,0,2,0,2)))
  
  # Check result: Dropin.
  expect_that(res$Dropin[res$Sample.Name=="major_minor_1"],
              equals(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)))
  expect_that(res$Dropin[res$Sample.Name=="major_minor_2"],
              equals(c(1,0,0,1,0,0,1,0,0,0,0,0,1,0,0,1,0)))
  
})