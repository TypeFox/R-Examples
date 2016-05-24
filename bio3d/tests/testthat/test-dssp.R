context("Testing dssp()")

test_that("SSE assignment still works", {
  skip_on_cran()

  if(!check.utility('dssp')) {
     skip('Need DSSP installed to run this test')
  }

  ## Simple test with PDB ID 1HEL
  invisible(capture.output(pdb <- read.pdb("3ERJ")))
  sse <- dssp(pdb)
  
  ## helices
  sse.stored   <- c(17,  37,  58, 101,  19,  37, 58, 101)
  expect_that(as.numeric(sse$helix$start), equals(sse.stored))

  sse.stored   <- c(18,  9, 14,  8, 16, 10, 13,  7)
  expect_that(as.numeric(sse$helix$length), equals(sse.stored))

  ## sheet
  sse.stored   <- c(3, 50, 75, 93,  3, 50, 75, 93 )
  expect_that(as.numeric(sse$sheet$start), equals(sse.stored))

  sse.stored   <- c(8, 6, 4, 8, 8, 6, 4, 8)
  expect_that(as.numeric(sse$sheet$length), equals(sse.stored))

  
  ## With RESNO=FALSE
  sse <- dssp(pdb, resno=FALSE)
  
  ## helices
  sse.stored   <- c(16,  36,  57, 100, 134, 152, 173, 216)
  expect_that(as.numeric(sse$helix$start), equals(sse.stored))

  sse.stored <- c(rep("A", 4), rep("B", 4))
  expect_that(as.character(sse$helix$chain), equals(sse.stored))


  ## sheet
  sse.stored   <- c(2,  49,  74,  92, 118, 165, 190, 208)
  expect_that(as.numeric(sse$sheet$start), equals(sse.stored))

  sse.stored <- c(rep("A", 4), rep("B", 4))
  expect_that(as.character(sse$sheet$chain), equals(sse.stored))



  ## With FULL=TRUE
  sse <- dssp(pdb, full=TRUE)
  
  expect_that(sum(as.numeric(sse$hbonds[,"BP1"]), na.rm=T), equals(2127))
  expect_that(sum(as.numeric(sse$hbonds[,"BP2"]), na.rm=T), equals(1355))

  expect_that(sum(as.numeric(sse$hbonds[,"NH-O.1"]), na.rm=T), equals(12017))
  expect_that(sum(as.numeric(sse$hbonds[,"E1"]), na.rm=T), equals(-315.8))
  expect_that(sum(as.numeric(sse$hbonds[,"O-HN.1"]), na.rm=T), equals(13347))
  expect_that(sum(as.numeric(sse$hbonds[,"E2"]), na.rm=T), equals(-313.4))
  expect_that(sum(as.numeric(sse$hbonds[,"NH-O.2"]), na.rm=T), equals(11859))
  expect_that(sum(as.numeric(sse$hbonds[,"E3"]), na.rm=T), equals(-51.2))
  expect_that(sum(as.numeric(sse$hbonds[,"O-HN.2"]), na.rm=T), equals(13076))
  expect_that(sum(as.numeric(sse$hbonds[,"E4"]), na.rm=T), equals(-53.6))
  
  expect_that(length(which(sse$hbonds[,"Chain1"]=="A")), equals(112))
  expect_that(length(which(sse$hbonds[,"Chain1"]=="B")), equals(105))
  expect_that(length(which(sse$hbonds[,"Chain2"]=="A")), equals(115))
  expect_that(length(which(sse$hbonds[,"Chain3"]=="B")), equals(108))
  expect_that(length(which(sse$hbonds[,"Chain4"]=="A")), equals(116))
  
}
          )
