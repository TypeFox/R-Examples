context("Testing nma()")


test_that("NMA", {

  "mysign" <- function(a,b) {
    if(all(sign(a)==sign(b)))
      return(1)
    else
      return(-1)
  }

  ## Simple test with PDB ID 1HEL
  file <- system.file("examples/1hel.pdb",package="bio3d")
  invisible(capture.output(pdb <- read.pdb(file)))
  
  ## Calculate modes with default arguments
  invisible(capture.output(modes <- nma(pdb, ff='calpha',
                                        mass=TRUE, temp=300.0)))

  ## Check first eigenvector
  U7 <- c(-0.05471209, -0.054333625, 0.001052514,
          -0.041171891, -0.049232935, -0.001588035)
  nowU7 <- head(modes$U[,7])
  expect_that(nowU7 * mysign(U7, nowU7), equals(U7, tolerance=1e-6))
  
  ## Check second eigenvector
  U8 <- c(0.064185522, 0.027349834, -0.024359816,
          0.011493963, 0.029426825, -0.014397686)
  nowU8 <- head(modes$U[,8])
  expect_that(nowU8 * mysign(U8, nowU8), equals(U8, tolerance=1e-6))

  ## Check Mode vector
  mode7 <- c(-0.092579348, -0.091938941, 0.001780978,
             -0.079838481, -0.095470057, -0.003079439)
  nowMode7 <- head(modes$modes[,7])
  expect_that(nowMode7 * mysign(mode7, nowMode7), equals(mode7, tolerance=1e-6))

  ## Check eigenvalues
  eival <- c(0.013383, 0.013933, 0.022355, 0.025518, 0.029944, 0.033954)
  nowEival <- modes$L[7:12]
  expect_that(nowEival, equals(eival, tolerance=1e-6))

  ## Check frequencies
  freqs <- c(0.018411826, 0.018786352, 0.023796192,
             0.025423975, 0.027540704, 0.029326863)
  nowFreqs <- modes$frequencies[7:12]
  expect_that(nowFreqs, equals(freqs, tolerance=1e-6))

  ## Dimensions
  expect_that(dim(modes$U), equals(c(387, 387)))
  expect_that(dim(modes$modes), equals(c(387, 387)))
  expect_that(length(modes$L), equals(387))
  expect_that(length(modes$frequencies), equals(387))
  expect_that(length(modes$mass), equals(129))
  expect_that(modes$natoms, equals(129))
  expect_that(modes$temp, equals(300))

  ## Orthognals
  expect_that(as.numeric(modes$U[,7] %*% modes$U[,7]),
              equals(1, tolerance=1e-6))
  expect_that(as.numeric(modes$U[,7] %*% modes$U[,8]),
              equals(0, tolerance=1e-6))
  expect_that(all((round(c(modes$U[,7] %*% modes$U),6)==0)[-7]),
              equals(TRUE))
  
  expect_that(all(round(c(modes$L[1:6]), 6)==0), equals(TRUE))
  

  ###################################################################
  #
  # Test with ouput from MMTK
  #
  ###################################################################
  
  
  "calpha.mmtk" <- function(r, ...) {
    ## MMTK Units: kJ / mol / nm^2
    a <- 128; b <- 8.6 * 10^5; c <- 2.39 * 10^5;
    ifelse( r<4.0,
           b*(r/10) - c,
           a*(r/10)^(-6) )
  }

  ## Vibrational Modes
  invisible(capture.output(modes <- nma(pdb, pfc.fun=calpha.mmtk,
                                        mmtk=TRUE, addter=FALSE)))

  ## Mode vector 7 (mmtk: modes[6])
  mmtk7 <- c(0.009399498664059314, 0.009162216956173577, -0.00018940255982217028,
             0.008013487647313355, 0.009521462401750403, 0.000300410055738782,
             0.006725323170414416, 0.00613075811499374, 0.0007167801244317134,
             0.003911230038334056, 0.0031036402193391484, 0.00011224732577516142,
             5.015756380626851e-05, 0.00122307913030356, 0.0005064454471294721,
             0.0014876013838084666, -0.003968053761191632, 0.00020389385408319644)
  nowMmtk7 <- head(modes$modes[,7], n=18)
  expect_that(nowMmtk7 * mysign(mmtk7, nowMmtk7), equals(mmtk7, tolerance=1e-6))

  ## Raw mode vector 7 (mmtk: modes.rawMode(6))
  mmtk7 <- c(0.05535176862194779, 0.053954464078107375, -0.0011153538121951093,
             0.04133856415826275, 0.049117637874840574, 0.0015497023155839553,
             0.04227251082807122, 0.03853532867244931, 0.00450537391343214)
  nowMmtk7 <- head(modes$U[,7], n=9)
  expect_that(nowMmtk7 * mysign(mmtk7, nowMmtk7), equals(mmtk7, tolerance=1e-6)) 
 
  ## Frequencies
  mmtkFreqs <- c(0.18417800523842359, 0.18804324107310424, 0.23820080688206749,
                 0.25592672017449125, 0.2798133442063071, 0.29367413814307064)
  nowMmtkFreqs <- modes$frequencies[7:12]
  expect_that(nowMmtkFreqs, equals(mmtkFreqs, tolerance=1e-6))

  ## Fluctuations
  mmtk.flucts <- c(0.00195060853392, 0.00113764918589, 0.00167187530508,
                   0.00175346604072, 0.00151209078542, 0.00130098648001,
                   0.00133495588156, 0.00107978100112, 0.000924829566202,
                   0.00109689698409)
  nowFlucts <- modes$fluctuations[1:10]
  expect_that(nowFlucts, equals(mmtk.flucts, tolerance=1e-6))
  
  
  ## Energetic Modes (mass=FALSE)
  invisible(capture.output(modes <- nma(pdb, pfc.fun=calpha.mmtk, mass=FALSE,
                                        mmtk=TRUE, addter=FALSE)))
  
  mmtk7 <- c(0.010350805923938345, 0.009267077807430083, -3.701643999426641e-05,
             0.008268033266170226, 0.009606710315232818, 0.0003705525203545053,
             0.006767227535558591, 0.005694101052352917, 0.001077079483122824)
  nowMmtk7 <- head(modes$modes[,7], n=9) 
  expect_that(nowMmtk7 * mysign(mmtk7, nowMmtk7), equals(mmtk7, tolerance=1e-6))

  mmtk7 <- c(0.05478481030376396, 0.04904884735365174, -0.00019592084498809212,
             0.04376109816000157, 0.05084645641420892, 0.001961262696318724,
             0.03581762420652206, 0.030137773647408828, 0.005700773021792685)
  nowMmtk7 <- head(modes$U[,7], n=9)
  expect_that(nowMmtk7 * mysign(mmtk7, nowMmtk7), equals(mmtk7, tolerance=1e-6))

   ## Fluctuations
  mmtk.flucts <- c(0.00195600136572, 0.00114595965451, 0.00168855332538,
                   0.00175330685712, 0.00152428233485, 0.00130978174806,
                   0.00134381308059, 0.00108408194319, 0.000924316154921,
                   0.0010985505357)
  nowFlucts <- modes$fluctuations[1:10]
  expect_that(nowFlucts, equals(mmtk.flucts, tolerance=1e-6))


  
  ## Energetic Modes (mass=FALSE, temp=NULL)
  invisible(capture.output(modes <- nma(pdb, pfc.fun=calpha.mmtk, mass=FALSE, temp=NULL)))
 
  mmtk7 <- c(0.05478481030376396, 0.04904884735365174, -0.00019592084498809212,
             0.04376109816000157, 0.05084645641420892, 0.001961262696318724,
             0.03581762420652206, 0.030137773647408828, 0.005700773021792685)
  nowMmtk7 <- head(modes$modes[,7], n=9)
  expect_that(nowMmtk7 * mysign(mmtk7, nowMmtk7), equals(mmtk7, tolerance=1e-6))

  mmtk7 <- c(0.05478481030376396, 0.04904884735365174, -0.00019592084498809212,
             0.04376109816000157, 0.05084645641420892, 0.001961262696318724,
             0.03581762420652206, 0.030137773647408828, 0.005700773021792685)
  nowMmtk7 <- head(modes$U[,7], n=9)
  expect_that(nowMmtk7 * mysign(mmtk7, nowMmtk7), equals(mmtk7, tolerance=1e-6))
  
  ## Fluctuations
  mmtk.flucts <- c(0.000784175524735, 0.000459423765829, 0.000676953612192,
                   0.000702913785645, 0.000611096147848, 0.000525101264027,
                   0.00053874467886, 0.000434616530213, 0.00037056523503,
                   0.000440417096776)
  nowFlucts <- modes$fluctuations[1:10]
  expect_that(nowFlucts, equals(mmtk.flucts, tolerance=1e-6))

  
  
  ## ANM eigenvectors
  invisible(capture.output(modes <- nma(pdb, ff='anm', mass=FALSE, temp=NULL, cutoff=15)))
  
  anm7 <- c(0.041345308400364066, 0.03345000499525146, 0.008604839963113613,
            0.03755854024944313, 0.036973377719312125, 0.008638534251932818,
            0.033187347539802, 0.022779436981185324, 0.004702511702428035)
  nowAnm7 <- head(modes$modes[,7], n=9)
  expect_that(nowAnm7 * mysign(anm7, nowAnm7), equals(anm7, tolerance=1e-6))


  ## ANM eigenvalues check
  eivalsANM <- c(0.84962016107196869, 1.0327718030407862, 1.3724207202555807,
                 1.7545168246132175, 1.9606866740784614, 2.2429459260702607)
  nowEivalANM <- modes$L[7:12]
  expect_that(nowEivalANM, equals(eivalsANM, tolerance=1e-6))

  
  ###################################################################
  #
  # Test mass custom stuff
  #
  ###################################################################

  mc <- list(ALA=500, SER=1000)
  invisible(capture.output(modes <- nma(pdb, mass.custom=mc)))

  mass.expected <- c(500.000, 500.000, 500.000, 131.196, 129.180, 157.194)
  expect_that(modes$mass[9:14], equals(mass.expected, tolerance=1e-6))

  sum.expected <- 28564.36
  expect_that(sum(modes$mass), equals(sum.expected, tolerance=1e-6))

  modes.expected <- c(-0.128550854, -0.069409382,  0.011821391,
                      -0.056729257, -0.076231424, 0.004736013)
  nowMode7 <- modes$modes[1:6, 7]
  expect_that(nowMode7 * mysign(nowMode7, modes.expected), equals(modes.expected, tolerance=1e-6))
  
  L.expected <- c(0.007375, 0.009036, 0.013007, 0.015084, 0.020111, 0.022608)
  expect_that(modes$L[7:12], equals(L.expected, tolerance=1e-6))
  

  ###################################################################
  #
  # Test build.hessian
  #
  ###################################################################

  sele <- atom.select(pdb, chain="A", elety="CA")
  xyz <- pdb$xyz[sele$xyz]

  i <- 2; j <- 5;
  hessian <- build.hessian(xyz, pfc.fun=calpha.mmtk, fc.weights=NULL)

  subhess <- matrix(c(-79.568546,  80.648662, -19.298073,
                      80.648662, -81.743440,  19.560037,
                      -19.298073, 19.560037,  -4.680438),
                    ncol=3, byrow=TRUE)
  
  expect_that(hessian[atom2xyz(j), atom2xyz(i)],
              equals(subhess, tolerance=1e-6))
  
  ## Force constant weighting
  weight <- 0.5; 
  fc.mat <- matrix(1, nrow=length(xyz)/3, ncol=length(xyz)/3)
  fc.mat[i, j] <- weight; fc.mat[j, i] <- weight;
  hessian2 <- build.hessian(xyz, pfc.fun=calpha.mmtk, fc.weights=fc.mat)
  
  expect_that(hessian[atom2xyz(j), atom2xyz(i)] * weight,
              equals(hessian2[atom2xyz(j), atom2xyz(i)], tolerance=1e-6))
  
   
}
          )
