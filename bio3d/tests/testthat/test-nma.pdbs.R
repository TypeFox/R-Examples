context("Testing nma.pdbs()")

test_that("eNMA works", {

  skip_on_cran()

  "mysign" <- function(a,b) {
    if(all(sign(a)==sign(b)))
      return(1)
    else
      return(-1)
  }

  attach(transducin)
  inds <- unlist(lapply(c("1TND_A", "1TAG", "1AS0", "1AS2"), grep, pdbs$id))
  pdbs <- trim.pdbs(pdbs, row.inds=inds)
  gaps <- gap.inspect(pdbs$xyz)


  ## Calc modes
  invisible(capture.output(modes <- nma.pdbs(pdbs, fit=TRUE, rm.gaps=TRUE, ncore=1)))

  ## check dimensions
  expect_that(dim(modes$U), equals(c(939, 933, 4)))
  expect_that(dim(modes$L), equals(c(4, 933)))
  expect_that(dim(modes$fluctuations), equals(c(4, 313)))
  
  ## structure 1- mode1:
  U1 <- c(-0.046380187,  0.007816393,  0.078933552, -0.041820111, 0.009819448, 0.047598954)
  nowU1 <- head(modes$U.subspace[,1,1], n=6)
  expect_that(nowU1 * mysign(U1, nowU1), equals(U1, tolerance=1e-6))

  ## structure 1- mode2:
  U2 <- c(0.007819148, -0.005717370, -0.051933927, 0.002265788, -0.013386437, -0.039409712)
  nowU2 <- head(modes$U.subspace[,2,1], n=6)
  expect_that(nowU2 * mysign(U2, nowU2), equals(U2, tolerance=1e-6))

  ## structure 4- mode3:
  U3 <- c(-0.13697744, -0.05401178,  0.09615700, -0.11026525, -0.01510157, 0.05796091)
  nowU3 <- head(modes$U.subspace[,3,4], n=6)
  expect_that(nowU3 * mysign(U3, nowU3), equals(U3, tolerance=1e-6))

  ## structure 4-mode1 - tail:
  U1 <- c(0.009820174, 0.004566910, -0.055544781, 0.009938013, 0.006707436, -0.068154672)
  nowU1 <- tail(modes$U.subspace[,1,4], n=6)
  expect_that(nowU1 * mysign(U1, nowU1), equals(U1, tolerance=1e-6))


  ## Fluctuations:
  f1 <- c(0.3657288, 0.2196504, 0.1449100, 0.1217517, 0.1130416, 0.1007862)
  f2 <- c(0.5295304, 0.3004094, 0.2011062, 0.1401276, 0.1205508, 0.1011349)
  f4 <- c(0.6488982, 0.3079024, 0.2166070, 0.1504672, 0.1282899, 0.1029634)

  expect_that(modes$fluctuations[1,1:6], equals(f1, tolerance=1e-6))
  expect_that(modes$fluctuations[2,1:6], equals(f2, tolerance=1e-6))
  expect_that(modes$fluctuations[4,1:6], equals(f4, tolerance=1e-6))
  
  ## Orthognal
  expect_that(as.numeric(modes$U.subspace[,1,1] %*% modes$U.subspace[,1,1]),
              equals(1, tolerance=1e-6))
  expect_that(as.numeric(modes$U.subspace[,1,1] %*% modes$U.subspace[,2,1]),
              equals(0, tolerance=1e-6))

  ## RMSIP
  rmsips <- c(1.0000, 0.9174, 0.9441, 0.9251)
  expect_that(as.vector(modes$rmsip[1,]), equals(rmsips, tolerance=1e-6))

  
  ## Multicore (same arguments as above!)
  invisible(capture.output(mmc <- nma.pdbs(pdbs, fit=TRUE, rm.gaps=TRUE, ncore=NULL)))
  expect_that(mmc$fluctuations, equals(modes$fluctuations, tolerance=1e-6))
  expect_that(mmc$U.subspace, equals(modes$U.subspace, tolerance=1e-6))


  
  ## Calc modes with rm.gaps=FALSE
  invisible(capture.output(modes <- nma.pdbs(pdbs, fit=TRUE, rm.gaps=FALSE, ncore=NULL)))

  ## structure 1-mode1 - tail:
  U1 <- c(0.04397369, -0.01400912, -0.02123377, 0.08124317, -0.02660929, -0.02619898)
  nowU1 <- tail(modes$U.subspace[,1,1], n=6)
  expect_that(nowU1 * mysign(U1, nowU1), equals(U1, tolerance=1e-6))

  ## structure 1-mode1 - tail:
  U1 <- c(0.001300939, -0.053317847,  0.011292943, 0.003307034, -0.071684004, NA)
  nowU1 <- modes$U.subspace[938:943,1,2]
  U1[is.na(U1)] <- 0
  nowU1[is.na(nowU1)] <- 0
  expect_that(nowU1 * mysign(nowU1, U1), equals(U1, tolerance=1e-6))

  ## fluctuations
  na.expected <- c(3, 4, 1258, 1259, 1262, 1263, 1266, 1267, 1268, 1270, 1271, 1272, 1274, 1275, 1276,
                   1278, 1279, 1280, 1282, 1283, 1284, 1286, 1287, 1288, 1290, 1291, 1292)
                   
  expect_that(which(is.na(modes$fluctuations)), equals(na.expected))
  
 
  f1 <- c(0.59967448, 0.34438649, 0.20382435, 0.13350449)
  f4 <- c(0.3335200, 0.4255609, 0.5941589, rep(NA, 7))
          
  expect_that(modes$fluctuations[1,1:4], equals(f1, tolerance=1e-6))
  expect_that(tail(modes$fluctuations[4,], n=10), equals(f4, tolerance=1e-6))


   ## Calc modes with mass=FALSE and temp=NULL
  invisible(capture.output(modes <- nma.pdbs(pdbs, mass=FALSE, temp=NULL, ncore=NULL)))

  ## structure 1- mode1:
  U1 <- c(-0.04043330,  0.00730273,  0.07000757, -0.04520831,  0.01130271,  0.05337233)
  nowU1 <- head(modes$U.subspace[,1,1], n=6)
  expect_that(nowU1 * mysign(U1, nowU1), equals(U1, tolerance=1e-6))
  
  ## structure 1- mode2:
  U2 <- c(-0.002813312,  0.005765808,  0.039807147,  0.001667637,  0.014587327, 0.038682763)
  nowU2 <- head(modes$U.subspace[,2,1], n=6)
  expect_that(nowU2 * mysign(U2, nowU2), equals(U2, tolerance=1e-6))

  ## structure 5- mode3:
  U3 <- c(0.11324262, 0.04220159, -0.07597465,  0.10267147,  0.01483591, -0.05261675)
  nowU3 <- head(modes$U.subspace[,3,4], n=6)
  expect_that(nowU3 * mysign(U3, nowU3), equals(U3, tolerance=1e-6))


  
  ## Calc modes with mass=FALSE and temp=NULL and ff="anm"
  invisible(capture.output(modes <- nma.pdbs(pdbs, mass=FALSE, temp=NULL, ff="anm", ncore=NULL)))

  ## structure 3- mode10:
  U1 <- c(0.03630660,  0.03078575, -0.02376714,  0.01906218,  0.01110582, -0.01361602)
  nowU1 <- head(modes$U.subspace[,10,3], n=6)
  expect_that(nowU1 * mysign(U1, nowU1), equals(U1, tolerance=1e-6))
  
  ## structure 4- mode1:
  U1 <- c(-0.04113844,  0.01096919,  0.07368620, -0.04250786,  0.01320282,  0.05550216)
  nowU1 <- head(modes$U.subspace[,1,4], n=6)
  expect_that(nowU1 * mysign(U1, nowU1), equals(U1, tolerance=1e-6))

  f1 <- c(0.3630744, 0.2768045, 0.1996179, 0.1766148)
  f2 <- c(0.5231570, 0.3519813, 0.2331503, 0.2003372)
  expect_that(modes$fluctuations[1,1:4], equals(f1, tolerance=1e-6))
  expect_that(modes$fluctuations[2,1:4], equals(f2, tolerance=1e-6))

  detach(transducin)
})
          

  
