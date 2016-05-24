context("Testing basic PDB structure operation")

test_that("read.pdb() reads a normal pdb file", {

  ## Simple test with PDB ID 1HEL
  file <- system.file("examples/1dpx.pdb",package="bio3d")
  invisible(capture.output(pdb <- read.pdb(file)))
  
  expect_is(pdb$atom, "data.frame")
  expect_true(inherits(pdb, "pdb"))
  expect_true(inherits(pdb$xyz, "xyz"))

  expect_equal(nrow(pdb$atom), 1177)
  expect_equal(sum(pdb$calpha), 129)

  expect_equal(sum(pdb$atom$resid=="HOH"), 177)
  expect_equal(sum(pdb$atom$resid=="CL"), 2)
  expect_that(sum(pdb$xyz), equals(44657.12, tolerance=1e-6))

  expect_equal(sum(pdb$atom$type=="ATOM"), 998)
  expect_equal(sum(pdb$atom$type=="HETATM"), 179)

  expect_equal(pdb$remark$biomat$num, 1)
  expect_equal(pdb$remark$biomat$chain[[1]], "A")
  true_mat <- matrix(c(1.0, 0.0, 0.0, 0.0,
                       0.0, 1.0, 0.0, 0.0,
                       0.0, 0.0, 1.0, 0.0), nrow=3, byrow=TRUE)
  expect_equivalent(pdb$remark$biomat$mat[[1]][[1]], true_mat)

  invisible(capture.output(spdb <- read.pdb(file, ATOM.only=TRUE)))
  expect_equal(spdb$atom, pdb$atom)
  expect_equal(spdb$xyz, pdb$xyz)
  expect_equal(spdb$calpha, pdb$calpha)
  expect_true(is.null(spdb$helix)) 
  expect_true(is.null(spdb$sheet)) 
  expect_true(is.null(spdb$seqres))
  expect_true(is.null(spdb$remark))
  
})


test_that("read.pdb() reads and stores data properly", {
  skip_on_cran()
  
  datdir <- tempdir()
  invisible(capture.output(get.pdb(c("3DRC", "1P3Q", "1SVK", "1L2Y"), path=datdir,
                                   overwrite = FALSE, verbose = FALSE)))
  
   # "3DRC" example PDB has a CA calcium ion and a CA containing ligand. 
   expect_error(read.pdb("nothing"))
   invisible(capture.output(pdb <- read.pdb(file.path(datdir, "3DRC.pdb"))))

   expect_equal(nrow(pdb$atom), 2954)
   expect_equal(sum(pdb$calpha), 318)
#   expect_equivalent(aa321(pdb$seqres), pdbseq(pdb))
   expect_equal(pdb$xyz[1:6], c(24.317, 59.447, 4.079, 25.000, 58.475, 4.908), tolerance=1e-6)
   expect_equal(length(pdb$helix$start), 8)
   expect_equal(length(pdb$sheet$start), 16)
    
   # "1SVK" example PDB has alternate location indicator 
   invisible(capture.output(pdb <- read.pdb(file.path(datdir, "1SVK.pdb"))))
   expect_equal(sum(pdb$calpha), 313)
   expect_equal(sum(pdb$atom$resno==47), 6)
   expect_equal(sum(pdb$atom$resid=="GDP"), 28)
    
   # multi-model structure
   invisible(capture.output(pdb <- read.pdb(file.path(datdir, "1L2Y.pdb"), multi=TRUE)))
   expect_equal(dim(pdb$xyz), c(38, 912))
   expect_equal(pdb$xyz[20, 1:6], c(-8.559, 6.374, -1.226, -7.539, 6.170, -0.168), tolerance=1e-6)

   # one atom
   cat("ATOM      1  N   SER Q 398      48.435  21.981  -6.393  1.00 56.10           N\n",
      file=file.path(datdir, "t1a.pdb"))
   pdb <- read.pdb(file.path(datdir, "t1a.pdb"))
   expect_is(pdb$atom, "data.frame") 



  ### write.pdb()
   invisible(capture.output(pdb <- read.pdb(file.path(datdir, "3DRC.pdb"))))
   write.pdb(pdb, file=file.path(datdir, "t1.pdb"))
   invisible(capture.output(pdb1 <- read.pdb(file.path(datdir, "t1.pdb"))))
   expect_identical(pdb$atom, pdb1$atom)
   expect_identical(pdb$xyz, pdb1$xyz)
   expect_identical(pdb$calpha, pdb1$calpha)
 
   # multi-model structure
   invisible(capture.output(pdb <- read.pdb(file.path(datdir, "1L2Y.pdb"), multi=TRUE)))
   write.pdb(pdb, file=file.path(datdir, "t2.pdb"))
   invisible(capture.output(pdb2 <- read.pdb(file.path(datdir, "t2.pdb"), multi=TRUE)))
   # SSE and SEQRES missing in write.pdb()
   pdb[c("seqres", "helix", "sheet", "call")] <- NULL
   pdb2[c("seqres", "helix", "sheet", "call")] <- NULL
   expect_identical(pdb, pdb2)
  


  ### trim.pdb() 
   invisible(capture.output(pdb <- read.pdb(file.path(datdir, "1P3Q.pdb"))))
   pdb1 <- trim.pdb(pdb, inds = atom.select(pdb, "calpha", verbose=FALSE))
   expect_is(pdb1, "pdb")
   expect_equal(nrow(pdb1$atom), 228)
   expect_equal(sum(pdb1$calpha), 228)
   expect_equivalent(pdb1$helix$start, pdb$helix$start)
   expect_equivalent(sort(pdb1$sheet$end), sort(pdb$sheet$end))

   pdb2 <- trim.pdb(pdb, inds = atom.select(pdb, "protein", chain="U", verbose=FALSE))
   expect_equal(nrow(pdb2$atom), 593)
   expect_equal(sum(pdb2$calpha), 74)
   expect_equivalent(pdb2$helix, list(start=c(22, 37, 56), end=c(35, 39, 60), 
       chain=rep("U",3), type=c("1", "5", "5")))
   expect_equivalent(pdb2$sheet, list(start=c(12,2,66,41,48), end=c(16,7,71,45,49), 
       chain=rep("U",5), sense=c("0","-1","1","-1","-1")))
})


