context("Testing atom.select function")

test_that("atom.select() gets correct selections with various options", {

  # Use the curated test-purpose pdb file
  file <- system.file("examples/test.pdb", package="bio3d")
#  invisible(capture.output(pdb <- read.pdb(file, rm.alt = FALSE)))
  invisible(capture.output(pdb <- read.pdb(file)))

  # Select all atoms: omit 1 ALT
  capture.output(all.inds <- atom.select(pdb, "all"))
  expect_equal(length(all.inds$atom), 175)

  # Select chain A: return everything in chain A
  capture.output(a.inds <- atom.select(pdb, chain = "A"))
  expect_equal(a.inds$xyz, atom2xyz(a.inds$atom))  # self-consistent
  expect_equal(length(a.inds$atom), 98)
  expect_equal(a.inds$atom[1], 1)
  expect_equal(a.inds$atom[98], 170)

  # Select protein: omit 1 ALT
  #                 omit 2 'unknown' AA
  capture.output(pro.inds <- atom.select(pdb, "protein")$atom)
  expect_equal(length(pro.inds), 103)

  # Select C-alpha: omit 2 'unknown' AA
  #                 omit 1 with missing 'CA'
  capture.output(ca.inds <- atom.select(pdb, "calpha")$atom)
  expect_equal(length(ca.inds), 7)
  expect_equal(ca.inds[c(1, 5, 7)], c(3, 53, 110))

  capture.output(ca2.inds <- atom.select(pdb, elety = "CA")$atom)
  expect_equal(length(ca2.inds), 10) # include calcium

  # Select 'N': return the number of all amino acids
  capture.output(cb.inds <- atom.select(pdb, elety = "N")$atom)
  expect_equal(length(cb.inds), 10)
  
  # Select 'unknown' AA
  capture.output(unk.inds <- atom.select(pdb, resid = "TES")$atom)
  expect_equal(length(unk.inds), 26)
  expect_equal(unk.inds[c(3, 13)], c(10, 94))

  # Select 'ATOM' record: omit 1 'ALT'
  capture.output(ATOM.inds <- atom.select(pdb, type = "ATOM")$atom)
  expect_equal(length(ATOM.inds), 129)
 
  # Select first 5 CA atoms by 'eleno'
  capture.output(ca5.inds <- atom.select(pdb, eleno = c(3, 10, 20, 27, 42))$atom)
  expect_equal(length(ca5.inds), 5)
  expect_true(all(pdb$atom[ca5.inds, "elety"] == "CA"))

  # Select hydrogen/ligand/water
  capture.output(h.inds <- atom.select(pdb, "h")$atom)
  expect_equal(length(h.inds), 77)
  
  capture.output(lig.inds <- atom.select(pdb, "ligand")$atom)
  expect_equal(length(lig.inds), 67) # include "TES"

  capture.output(wat.inds <- atom.select(pdb, "water")$atom)
  expect_equal(length(wat.inds), 5)

  capture.output(ion.inds <- atom.select(pdb, resid = "CA")$atom)
  expect_equal(ion.inds, 170)

  capture.output(gdp.inds <- atom.select(pdb, resid = "GDP")$atom)
  expect_equal(length(gdp.inds), 40)
  
  # More string test
  capture.output(bb.inds <- atom.select(pdb, "backbone")$atom)
  capture.output(bb2.inds <- atom.select(pdb, "back")$atom)
  expect_equal(bb.inds, bb2.inds)
  expect_equal(length(bb.inds), 31)
  
  capture.output(cb2.inds <- atom.select(pdb, "cbeta")$atom)
  expect_equal(length(cb2.inds), 36)
  
  capture.output(npro.inds <- atom.select(pdb, "notprotein")$atom)
  capture.output(npro2.inds <- atom.select(pdb, "protein", inverse = TRUE)$atom)
  expect_equal(npro.inds, npro2.inds)
  expect_equal(length(intersect(pro.inds, npro.inds)), 0)
  expect_equal(length(pro.inds) + length(npro.inds), nrow(pdb$atom)) # omit ALT
  
  capture.output(nwat.inds <- atom.select(pdb, "notwater")$atom)
  capture.output(nwat2.inds <- atom.select(pdb, "water", inverse = TRUE)$atom)
  expect_equal(nwat.inds, nwat2.inds)
  expect_equal(length(intersect(wat.inds, nwat.inds)), 0)
  expect_equal(length(wat.inds) + length(nwat.inds), nrow(pdb$atom)) # omit ALT
   
  capture.output(noh.inds <- atom.select(pdb, "noh")$atom)
  capture.output(noh2.inds <- atom.select(pdb, "h", inverse = TRUE)$atom)
  expect_equal(noh.inds, noh2.inds)
  expect_equal(length(intersect(h.inds, noh.inds)), 0)
  expect_equal(length(h.inds) + length(noh.inds), nrow(pdb$atom)) # omit ALT
  
  # Test on combination of select
  capture.output(comb1.inds <- atom.select(pdb, chain = "B", resno = c(1,4)))
  expect_equal(length(comb1.inds$atom), 33) # omit ALT
  
  capture.output(comb2.inds <- atom.select(pdb, resid = "GDP", elety = "PA") )
  expect_equal(comb2.inds$atom, 135)

  capture.output(comb3.inds <- atom.select(pdb, "noh", resid = "TES", chain = "A") )
  expect_equal(comb3.inds$atom, c(8, 10, 12, 16, 17)) 
  
  capture.output(comb4.inds <- atom.select(pdb, chain = "B", resid = "GDP", operator = "OR") )
  expect_equal(length(comb4.inds$atom), 117)
})
