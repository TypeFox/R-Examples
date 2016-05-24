context("Testing core.find function")


test_that("core.find() works properly", {
  skip_on_cran()

  attach(transducin)
  inds <- unlist(lapply(c("1TND_A", "1TAG", "1AS0", "1AS2"), grep, pdbs$id))
  pdbs <- trim.pdbs(pdbs, row.inds=inds)
 
  invisible(capture.output(core <- core.find(pdbs, ncore=1)))
  resnos.1 <- c(202, 206, 209, 205, 203, 201)
  resnos.2 <- c(332, 334, 335, 336, 337, 340)
  
  expect_equal(length(core$resno), 313)
  expect_equal(resnos.1, as.numeric(core$resno[1:6]))
  expect_equal(resnos.2, as.numeric(tail(core$resno)))
  
  xyz <- c(16, 17,  18,  19,  20,  21,  25,  26,  27, 34)
  expect_equal(xyz, core$xyz[1:10])
  expect_equal(sum(core$xyz), 234006)
  
  ## Check multicore 
  invisible(capture.output(core.mc <- core.find(pdbs, ncore=NULL)))
  expect_identical(core, core.mc)

  detach(transducin)
})
