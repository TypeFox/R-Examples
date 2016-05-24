context("Testing the utility function clean.pdb()")

test_that("clean.pdb() does nothing for 'clean' pdb by default", {

  file <- system.file("examples/1dpx.pdb", package="bio3d")
  invisible(capture.output(pdb <- read.pdb(file)))

  invisible(capture.output(npdb <- clean.pdb(pdb)))
 
  expect_true(is.null(npdb$log))

  npdb$call <- NULL
  npdb$log <- NULL
  pdb$call <- NULL
  expect_equal(pdb, npdb)

})


test_that("clean.pdb() does renumbering properly", {

  skip_on_cran()

  invisible(capture.output(pdb <- read.pdb("1tag")))
  invisible(capture.output(npdb <- clean.pdb(pdb, force.renumber = TRUE)))

  resno <- npdb$atom[npdb$calpha, "resno"]
  expect_equal(resno[1:10], 1:10)

  # A PDB with 'insert' residues: Should do automatic renumbering
  invisible(capture.output(pdb <- read.pdb("1a7l")))
  invisible(capture.output(npdb <- clean.pdb(pdb)))
  
  resno <- npdb$atom[npdb$calpha, "resno"]
  expect_equal(resno[1:10], 1:10)

  # Renumbering for each chain
  invisible(capture.output(npdb <- clean.pdb(pdb, consecutive = FALSE)))

  resno <- npdb$atom[npdb$calpha, "resno"]
  chain <- npdb$atom[npdb$calpha, "chain"]
  expect_equal(resno[chain=="B"][1:10], 1:10)
 
  # Is SSE update correct?
  invisible(capture.output(ss <- pdb2sse(pdb)))
  invisible(capture.output(nss <- pdb2sse(npdb)))

  expect_equal(as.character(ss), as.character(nss))
 
} )


test_that("clean.pdb() relabels chains properly (fix.chain = TRUE)", {

  file <- system.file("examples/1dpx.pdb", package="bio3d")
  invisible(capture.output(pdb <- read.pdb(file)))

  # remove chain ID
  pdb$atom[, "chain"] <- as.character(NA)
  pdb$helix$chain <- "" 
  pdb$sheet$chain <- "" 

  invisible(capture.output(npdb <- clean.pdb(pdb, fix.chain = TRUE)))

  expect_equal(npdb$atom[npdb$calpha, "chain"], rep("A", sum(npdb$calpha)))

  # A case with wrong chain labels but consecutive residue numbering
  file <- system.file("examples/hivp.pdb", package="bio3d")
  invisible(capture.output(pdb0 <- read.pdb(file)))

  # Manually renumbering all residues
  pdb <- pdb0
  pdb$atom[, "resno"] <- vec2resno(1:sum(pdb$calpha), 
                  paste(pdb$atom[, "resno"], pdb$atom[, "chain"], sep = "_") ) 

  # Label both chains as "A" 
  pdb$atom[, "chain"] <- "A"
    
  invisible(capture.output(npdb <- clean.pdb(pdb, consecutive = FALSE, force.renumber = TRUE, fix.chain = TRUE)))
  
  pdb0$call <- NULL
  npdb$call <- NULL
  npdb$log <- NULL
  expect_equal(pdb0, npdb)
 
} )


