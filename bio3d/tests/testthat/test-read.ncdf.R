context("Testing basic operation with NetCDF trajectory")

test_that("read.ncdf() and write.ncdf() works properly", {
   skip_on_cran()
   
   ##- Prepare files
   trjfile <- tempfile()
   file <- system.file("examples/hivp.dcd", package="bio3d")
   invisible(capture.output(trj0 <- read.dcd(file)))
   time0 <- sort(round(runif(nrow(trj0), 0, 1000), digit=3))
   cell0 <- matrix(rep(runif(6, 0, 100), nrow(trj0)), ncol=6, byrow=TRUE)
   rownames(trj0) <- time0

   ##- Write
   out <- try(write.ncdf(trj0, trjfile, cell = cell0))
   expect_false(inherits(out, "try-error"))

   ##- Read
   trj <- read.ncdf(trjfile, headonly = TRUE, verbose = FALSE)
   expect_output(str(trj), "frames: int 351" )
   expect_output(str(trj), "atoms : int 198" )

   trj <- read.ncdf(trjfile, cell = TRUE, verbose = FALSE)
   ##expect_equal(trj, as.data.frame(cell0, stringsAsFactors=FALSE), tolerance = 1e-6)
   expect_equal(trj, cell0, tolerance = 1e-6)

   trj <- read.ncdf(trjfile, verbose = FALSE, time = TRUE)
   expect_equal(as.numeric(rownames(trj)), time0, tolerance = 1e-6)
   expect_equivalent(trj, trj0)

   pdb <- read.pdb(system.file("examples/hivp.pdb", package="bio3d"))
   inds <- atom.select(pdb, chain="A", verbose=FALSE)
   trj <- read.ncdf(trjfile, verbose = FALSE, first=10, last=20, stride=2,
      at.sel = inds)
   expect_equivalent(trj, as.xyz(trj0[seq(10, 20, 2), inds$xyz]))

   # multiple files
   files <- rep(trjfile, 4)
   txt <- capture.output(trj <- read.ncdf(files, headonly = TRUE))
   expect_output(txt, "Frames: 1404")
   expect_output(txt, "Atoms: 198")
})

