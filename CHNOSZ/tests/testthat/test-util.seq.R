context("util.seq")

test_that("count.aa() warns about unrecognized amino acids and performs substring operations", {
  expect_message(count.aa("ABCDEFGHIJ"), "count.aa: unrecognized letter\\(s\\) in protein sequence: B J")
  myseq <- "AAAAAGGGGG"
  expect_equal(count.aa(myseq, stop=5)[, "G"], 0)
  expect_equal(count.aa(myseq, start=6)[, "A"], 0)
  expect_equal(count.aa(myseq, start=5, stop=6)[, c("A", "G")], c(1, 1), check.attributes=FALSE)
})

test_that("nucleobase sequences can be processed with count.aa(), nucleic.formula() and nucleic.complement()", {
  expect_message(dna <- count.aa("ABCDEFGHIJ", type="DNA"), "count.aa: unrecognized letter\\(s\\) in DNA sequence: B D E F H I J")
  expect_equal(as.numeric(dna), c(1, 1, 1, 0))
  expect_equal(nucleic.formula(dna), "C14H15N13O2")
  # nucleobases can be in any order
  expect_equal(nucleic.formula(dna[, 4:1, drop=FALSE]), "C14H15N13O2")
  # ACG -> UGC (RNA complement)
  expect_equal(nucleic.formula(nucleic.complement(dna, "RNA")), "C13H14N10O4")
})
