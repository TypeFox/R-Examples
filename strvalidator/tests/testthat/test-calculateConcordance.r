context("calculateConcordance")

################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG
# 28.04.2014: First tests for 'calculateConcordance'.
# 
# test_dir("inst/tests/")
# test_file("tests/testthat/test-calculateConcordance.r")
# test_dir("tests/testthat")

test_that("calculateConcordance", {

  # Create a dataframe for testing:
  # Load data.
  data(set2)
  
  # Set 1.
  x1 <- set2[set2$Marker!="TH01",]
  x1$Allele[11] <- "31.3"

  # Set 2.
  x2 <- set2[set2$Marker!="FGA" & set2$Marker!="vWA",]
  x2$Allele[2] <- ""
  
  # Set 3.
  x3 <- rbind(set2,set2)
  x3$Sample.Name[33:48] <- "SampleA03"
  x3$Sample.Name[49:64] <- "SampleA04"
  x3$Allele[6] <- "OL"
  x3$Allele[32] <- "15"

  # Set 4.
  x4 <- rbind(x2,x3[x3$Sample.Name=="SampleA03",])
  x4$Allele[41] <- "16"
  
  # Kit/dataset name vector.
  kitVector <- c("KitA", "KitB","KitC", "KitD")
  
  # List with dataset.
  dataList <- list(x1, x2, x3, x4)
  
  # TEST 01 -------------------------------------------------------------------
  # Test all differences with default settings.
  
  # Analyse dataframe.
  resList <- calculateConcordance(data=dataList, debug=FALSE)
  
  # Extract result tables.
  res1 <- resList[[1]]
  res2 <- resList[[2]]
  
  # Check return class.  
  expect_that(class(resList), matches(class(list())))
  expect_that(class(res1), matches(class(data.frame())))
  expect_that(class(res2), matches(class(data.frame())))
  
  # Check dimensions.
  expect_true(ncol(res1) == 6)
  expect_true(nrow(res1) == 8)
  
  expect_true(ncol(res2) == 6)
  expect_true(nrow(res2) == 6)
  
  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res1))
  expect_true("Marker" %in% names(res1))
  expect_true("Kit.1" %in% names(res1))
  expect_true("Kit.2" %in% names(res1))
  expect_true("Kit.3" %in% names(res1))
  expect_true("Kit.4" %in% names(res1))

  expect_true("Kits" %in% names(res2))
  expect_true("Samples" %in% names(res2))
  expect_true("Loci" %in% names(res2))
  expect_true("Alleles" %in% names(res2))
  expect_true("Discordances" %in% names(res2))
  expect_true("Concordance" %in% names(res2))
  
  # Check result.
  expect_true(res1$Sample.Name[1]=="SampleA01")
  expect_true(res1$Sample.Name[2]=="SampleA01")
  expect_true(res1$Sample.Name[3]=="SampleA01")
  expect_true(res1$Sample.Name[4]=="SampleA01")
  expect_true(res1$Sample.Name[5]=="SampleA01")
  expect_true(res1$Sample.Name[6]=="SampleA02")
  expect_true(res1$Sample.Name[7]=="SampleA02")
  expect_true(res1$Sample.Name[8]=="SampleA03")
  expect_true(res1$Marker[1]=="D3S1358")
  expect_true(res1$Marker[2]=="vWA")
  expect_true(res1$Marker[3]=="D2S1338")
  expect_true(res1$Marker[4]=="D18S51")
  expect_true(res1$Marker[5]=="FGA")
  expect_true(res1$Marker[6]=="vWA")
  expect_true(res1$Marker[7]=="FGA")
  expect_true(res1$Marker[8]=="D19S433")
  expect_true(res1$Kit.1[1]=="15,18")
  expect_true(res1$Kit.1[2]=="14")
  expect_true(res1$Kit.1[3]=="19")
  expect_true(res1$Kit.1[4]=="31.3,16")
  expect_true(res1$Kit.1[5]=="25")
  expect_true(res1$Kit.1[6]=="14")
  expect_true(res1$Kit.1[7]=="NA")
  expect_true(res1$Kit.1[8]=="NO SAMPLE")
  expect_true(res1$Kit.2[1]=="15,")
  expect_true(res1$Kit.2[2]=="NO MARKER")
  expect_true(res1$Kit.2[3]=="19")
  expect_true(res1$Kit.2[4]=="31.2,16")
  expect_true(res1$Kit.2[5]=="NO MARKER")
  expect_true(res1$Kit.2[6]=="NO MARKER")
  expect_true(res1$Kit.2[7]=="NO MARKER")
  expect_true(res1$Kit.2[8]=="NO SAMPLE")
  expect_true(res1$Kit.3[1]=="15,18")
  expect_true(res1$Kit.3[2]=="14")
  expect_true(res1$Kit.3[3]=="OL")
  expect_true(res1$Kit.3[4]=="31.2,16")
  expect_true(res1$Kit.3[5]=="25")
  expect_true(res1$Kit.3[6]=="14")
  expect_true(res1$Kit.3[7]=="15")
  expect_true(res1$Kit.3[8]=="15")
  expect_true(res1$Kit.4[1]=="15,")
  expect_true(res1$Kit.4[2]=="")
  expect_true(res1$Kit.4[3]=="19")
  expect_true(res1$Kit.4[4]=="31.2,16")
  expect_true(res1$Kit.4[5]=="")
  expect_true(res1$Kit.4[6]=="")
  expect_true(res1$Kit.4[7]=="")
  expect_true(res1$Kit.4[8]=="16")

  expect_true(res2$Kits[1]=="Kit.1 vs. Kit.2")
  expect_true(res2$Kits[2]=="Kit.1 vs. Kit.3")
  expect_true(res2$Kits[3]=="Kit.1 vs. Kit.4")
  expect_true(res2$Kits[4]=="Kit.2 vs. Kit.3")
  expect_true(res2$Kits[5]=="Kit.2 vs. Kit.4")
  expect_true(res2$Kits[6]=="Kit.3 vs. Kit.4")
  expect_true(res2$Samples[1]==2)
  expect_true(res2$Samples[2]==2)
  expect_true(res2$Samples[3]==2)
  expect_true(res2$Samples[4]==2)
  expect_true(res2$Samples[5]==2)
  expect_true(res2$Samples[6]==3)
  expect_true(res2$Loci[1]==8)
  expect_true(res2$Loci[2]==10)
  expect_true(res2$Loci[3]==10)
  expect_true(res2$Loci[4]==9)
  expect_true(res2$Loci[5]==9)
  expect_true(res2$Loci[6]==11)
  expect_true(res2$Alleles[1]==32)
  expect_true(res2$Alleles[2]==40)
  expect_true(res2$Alleles[3]==40)
  expect_true(res2$Alleles[4]==36)
  expect_true(res2$Alleles[5]==36)
  expect_true(res2$Alleles[6]==66)
  expect_true(res2$Discordances[1]==2)
  expect_true(res2$Discordances[2]==3)
  expect_true(res2$Discordances[3]==6)
  expect_true(res2$Discordances[4]==2)
  expect_true(res2$Discordances[5]==0)
  expect_true(res2$Discordances[6]==7)
  expect_true(res2$Concordance[1]==100*(32-2)/32)
  expect_true(res2$Concordance[2]==100*(40-3)/40)
  expect_true(res2$Concordance[3]==100*(40-6)/40)
  expect_true(res2$Concordance[4]==100*(36-2)/36)
  expect_true(res2$Concordance[5]==100*(36-0)/36)
  expect_true(res2$Concordance[6]==100*(66-7)/66)
  
  # TEST 02 -------------------------------------------------------------------
  # Test all differences with custom settings.
  
  # Analyse dataframe.
  resList <- calculateConcordance(data=dataList, kit.name=kitVector, 
                                  no.marker="M", no.sample="S", delimeter="|",
                                  debug=FALSE)
  
  # Extract result tables.
  res1 <- resList[[1]]
  res2 <- resList[[2]]
  
  # Check return class.  
  expect_that(class(resList), matches(class(list())))
  expect_that(class(res1), matches(class(data.frame())))
  expect_that(class(res2), matches(class(data.frame())))
  
  # Check dimensions.
  expect_true(ncol(res1) == 6)
  expect_true(nrow(res1) == 8)
  
  expect_true(ncol(res2) == 6)
  expect_true(nrow(res2) == 6)
  
  # Check that expected columns exist.  
  expect_true("Sample.Name" %in% names(res1))
  expect_true("Marker" %in% names(res1))
  expect_true("KitA" %in% names(res1))
  expect_true("KitB" %in% names(res1))
  expect_true("KitC" %in% names(res1))
  expect_true("KitD" %in% names(res1))
  
  expect_true("Kits" %in% names(res2))
  expect_true("Samples" %in% names(res2))
  expect_true("Loci" %in% names(res2))
  expect_true("Alleles" %in% names(res2))
  expect_true("Discordances" %in% names(res2))
  expect_true("Concordance" %in% names(res2))
  
  # Check result.
  expect_true(res1$Sample.Name[1]=="SampleA01")
  expect_true(res1$Sample.Name[2]=="SampleA01")
  expect_true(res1$Sample.Name[3]=="SampleA01")
  expect_true(res1$Sample.Name[4]=="SampleA01")
  expect_true(res1$Sample.Name[5]=="SampleA01")
  expect_true(res1$Sample.Name[6]=="SampleA02")
  expect_true(res1$Sample.Name[7]=="SampleA02")
  expect_true(res1$Sample.Name[8]=="SampleA03")
  expect_true(res1$Marker[1]=="D3S1358")
  expect_true(res1$Marker[2]=="vWA")
  expect_true(res1$Marker[3]=="D2S1338")
  expect_true(res1$Marker[4]=="D18S51")
  expect_true(res1$Marker[5]=="FGA")
  expect_true(res1$Marker[6]=="vWA")
  expect_true(res1$Marker[7]=="FGA")
  expect_true(res1$Marker[8]=="D19S433")
  expect_true(res1$KitA[1]=="15|18")
  expect_true(res1$KitA[2]=="14")
  expect_true(res1$KitA[3]=="19")
  expect_true(res1$KitA[4]=="31.3|16")
  expect_true(res1$KitA[5]=="25")
  expect_true(res1$KitA[6]=="14")
  expect_true(res1$KitA[7]=="NA")
  expect_true(res1$KitA[8]=="S")
  expect_true(res1$KitB[1]=="15|")
  expect_true(res1$KitB[2]=="M")
  expect_true(res1$KitB[3]=="19")
  expect_true(res1$KitB[4]=="31.2|16")
  expect_true(res1$KitB[5]=="M")
  expect_true(res1$KitB[6]=="M")
  expect_true(res1$KitB[7]=="M")
  expect_true(res1$KitB[8]=="S")
  expect_true(res1$KitC[1]=="15|18")
  expect_true(res1$KitC[2]=="14")
  expect_true(res1$KitC[3]=="OL")
  expect_true(res1$KitC[4]=="31.2|16")
  expect_true(res1$KitC[5]=="25")
  expect_true(res1$KitC[6]=="14")
  expect_true(res1$KitC[7]=="15")
  expect_true(res1$KitC[8]=="15")
  expect_true(res1$KitD[1]=="15|")
  expect_true(res1$KitD[2]=="")
  expect_true(res1$KitD[3]=="19")
  expect_true(res1$KitD[4]=="31.2|16")
  expect_true(res1$KitD[5]=="")
  expect_true(res1$KitD[6]=="")
  expect_true(res1$KitD[7]=="")
  expect_true(res1$KitD[8]=="16")
  
  expect_true(res2$Kits[1]=="KitA vs. KitB")
  expect_true(res2$Kits[2]=="KitA vs. KitC")
  expect_true(res2$Kits[3]=="KitA vs. KitD")
  expect_true(res2$Kits[4]=="KitB vs. KitC")
  expect_true(res2$Kits[5]=="KitB vs. KitD")
  expect_true(res2$Kits[6]=="KitC vs. KitD")
  expect_true(res2$Samples[1]==2)
  expect_true(res2$Samples[2]==2)
  expect_true(res2$Samples[3]==2)
  expect_true(res2$Samples[4]==2)
  expect_true(res2$Samples[5]==2)
  expect_true(res2$Samples[6]==3)
  expect_true(res2$Loci[1]==8)
  expect_true(res2$Loci[2]==10)
  expect_true(res2$Loci[3]==10)
  expect_true(res2$Loci[4]==9)
  expect_true(res2$Loci[5]==9)
  expect_true(res2$Loci[6]==11)
  expect_true(res2$Alleles[1]==32)
  expect_true(res2$Alleles[2]==40)
  expect_true(res2$Alleles[3]==40)
  expect_true(res2$Alleles[4]==36)
  expect_true(res2$Alleles[5]==36)
  expect_true(res2$Alleles[6]==66)
  expect_true(res2$Discordances[1]==2)
  expect_true(res2$Discordances[2]==3)
  expect_true(res2$Discordances[3]==6)
  expect_true(res2$Discordances[4]==2)
  expect_true(res2$Discordances[5]==0)
  expect_true(res2$Discordances[6]==7)
  expect_true(res2$Concordance[1]==100*(32-2)/32)
  expect_true(res2$Concordance[2]==100*(40-3)/40)
  expect_true(res2$Concordance[3]==100*(40-6)/40)
  expect_true(res2$Concordance[4]==100*(36-2)/36)
  expect_true(res2$Concordance[5]==100*(36-0)/36)
  expect_true(res2$Concordance[6]==100*(66-7)/66)
  
  # TEST 03 -------------------------------------------------------------------
  # Test with no discordance.
  
  # List with dataset.
  dataList <- list(x1, x1)

  # Analyse dataframe.
  resList <- calculateConcordance(data=dataList, debug=FALSE)
  
  # Extract result tables.
  res1 <- resList[[1]]
  res2 <- resList[[2]]
  
  # Check return class.  
  expect_that(class(resList), matches(class(list())))
  expect_that(class(res1), matches(class(data.frame())))
  expect_that(class(res2), matches(class(data.frame())))
  
  # Check dimensions.
  expect_true(ncol(res1) == 1)
  expect_true(nrow(res1) == 1)
  
  expect_true(ncol(res2) == 6)
  expect_true(nrow(res2) == 1)
  
  # Check that expected columns exist.  
  expect_true("Kits" %in% names(res2))
  expect_true("Samples" %in% names(res2))
  expect_true("Loci" %in% names(res2))
  expect_true("Alleles" %in% names(res2))
  expect_true("Discordances" %in% names(res2))
  expect_true("Concordance" %in% names(res2))
  
  # Check result.
  expect_true(res1== "NO DISCORDANCE")
  
  expect_true(res2$Kits[1]=="Kit.1 vs. Kit.2")
  expect_true(res2$Samples[1]==2)
  expect_true(res2$Loci[1]==10)
  expect_true(res2$Alleles[1]==40)
  expect_true(res2$Discordances[1]==0)
  expect_true(res2$Concordance[1]==100)
  
})