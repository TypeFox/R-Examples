# Copyright (c) 2016 Nuno Fachada
# Distributed under the MIT License (http://opensource.org/licenses/MIT)

library(micompr)
context("grpoutputs")

test_that("grpoutputs constructs the expected objects", {

  # Outputs names
  outputs <- c("PopSheep", "PopWolf", "QtyGrass",
               "EnSheep", "EnWolf", "EnGrass",
               "All")

  # Determine location of extdata files
  dir_nl_ok <- system.file("extdata", "nl_ok", package = "micompr")
  dir_jex_ok <- system.file("extdata", "j_ex_ok", package = "micompr")
  dir_jex_noshuff <- system.file("extdata", "j_ex_noshuff", package = "micompr")
  dir_jex_diff <- system.file("extdata", "j_ex_diff", package = "micompr")
  dir_na <- system.file("extdata", "testdata", "NA", package = "micompr")
  files <- glob2rx("stats400v1*.tsv")
  filesA_na <- glob2rx("stats400v1*n20A.tsv")
  filesB_na <- glob2rx("stats400v1*n20B.tsv")

  # Instantiate several grpoutputs objects
  go_ok <- grpoutputs(outputs, c(dir_nl_ok, dir_jex_ok),
                      c(files, files),
                      lvls = c("NLOK", "JEXOK"), concat = T)
  go_noshuff <- grpoutputs(outputs, c(dir_nl_ok, dir_jex_noshuff),
                           c(files, files),
                           lvls = c("NLOK", "JEXNOSHUF"), concat = T)
  go_diff <- grpoutputs(outputs, c(dir_nl_ok, dir_jex_diff),
                        c(files, files),
                        lvls = c("NLOK", "JEXDIFF"), concat = T)
  go_tri <- grpoutputs(6,
                       c(dir_nl_ok, dir_jex_noshuff, dir_jex_diff),
                       c(files, files, files))
  go_1out <- grpoutputs("OnlyOne", c(dir_nl_ok, dir_jex_ok),
                        c(files, files), concat = F)
  go_1lvl <- grpoutputs(3, dir_nl_ok, files)
  go_diflencatT <- grpoutputs(7, dir_na, c(filesA_na, filesB_na), concat = T)
  go_diflencatF <- grpoutputs(6, dir_na, c(filesA_na, filesB_na), concat = F)

  #### Start testing ####

  # Test if objects have the correct class
  expect_is(go_ok, "grpoutputs")
  expect_is(go_noshuff, "grpoutputs")
  expect_is(go_diff, "grpoutputs")
  expect_is(go_tri, "grpoutputs")
  expect_is(go_1out, "grpoutputs")
  expect_is(go_1lvl, "grpoutputs")
  expect_is(go_diflencatT, "grpoutputs")
  expect_is(go_diflencatF, "grpoutputs")

  # Test if output names are as expected
  expect_equal(names(go_ok$data), outputs)
  expect_equal(names(go_noshuff$data), outputs)
  expect_equal(names(go_diff$data), outputs)
  expect_equal(names(go_tri$data),
               c("out1", "out2", "out3", "out4", "out5", "out6"))
  expect_equal(names(go_1out$data), "OnlyOne")
  expect_equal(names(go_1lvl$data), c("out1", "out2", "out3"))
  expect_equal(names(go_diflencatT$data),
               c("out1", "out2", "out3", "out4", "out5", "out6", "out7"))
  expect_equal(names(go_diflencatF$data),
               c("out1", "out2", "out3", "out4", "out5", "out6"))

  # Test if length of concatenated output is equal to the sum of normal
  # outputs
  expect_equal(sum(sapply(go_ok$data[1:6], function(x) dim(x)[2])),
               dim(go_ok$data[[7]])[2])
  expect_equal(sum(sapply(go_noshuff$data[1:6], function(x) dim(x)[2])),
               dim(go_noshuff$data[[7]])[2])
  expect_equal(sum(sapply(go_diff$data[1:6], function(x) dim(x)[2])),
               dim(go_diff$data[[7]])[2])
  expect_equal(sum(sapply(go_diflencatT$data[1:6], function(x) dim(x)[2])),
               dim(go_diflencatT$data[[7]])[2])

  # Test if groups are as expected
  expect_equal(go_ok$groupsize, c(10, 10))
  expect_equal(go_noshuff$groupsize, c(10, 10))
  expect_equal(go_diff$groupsize, c(10, 10))
  expect_equal(go_tri$groupsize, c(10, 10, 10))
  expect_equal(go_1out$groupsize, c(10, 10))
  expect_equal(go_1lvl$groupsize, 10)
  expect_equal(go_diflencatT$groupsize, c(3, 3))
  expect_equal(go_diflencatF$groupsize, c(3, 3))

  # Test if levels are as expected
  expect_equal(go_ok$lvls, c("NLOK", "JEXOK"))
  expect_equal(go_noshuff$lvls, c("NLOK", "JEXNOSHUF"))
  expect_equal(go_diff$lvls, c("NLOK", "JEXDIFF"))
  expect_equal(go_tri$lvls, c(1, 2, 3))
  expect_equal(go_1out$lvls, c(1, 2))
  expect_equal(go_1lvl$lvls, 1)
  expect_equal(go_diflencatT$lvls, c(1, 2))
  expect_equal(go_diflencatF$lvls, c(1, 2))

  # Test if concat is as expected
  expect_true(go_ok$concat)
  expect_true(go_noshuff$concat)
  expect_true(go_diff$concat)
  expect_false(go_tri$concat)
  expect_false(go_1out$concat)
  expect_false(go_1lvl$concat)
  expect_true(go_diflencatT$concat)
  expect_false(go_diflencatF$concat)

})

test_that("grpoutputs throws errors when improperly invoked", {

  # OS-specific file separator
  fs <- .Platform$file.sep

  #### Start testing ####

  # Should throw error because number of levels != number of files
  expect_error(
    grpoutputs(4, c("dir1", "dir2"), glob2rx("*.tsv"), lvls = c("A", "B")),
    "Number of file sets is not the same as the given number of factor levels.",
    fixed = TRUE
  )

  # Should throw error because no files were found
  expect_error(
    grpoutputs(4, "some_fake_folder",
               c(glob2rx("fake_files*.csv"), glob2rx("also_fakes*.csv"))),
    paste("No files were found: some_fake_folder", fs,
          glob2rx("fake_files*.csv"),
          sep = ""),
    fixed = TRUE
  )

  # Should throw error because number of specified outputs is larger than
  # outputs available in file
  expect_error(
    grpoutputs(7, system.file("extdata", "nl_ok", package = "micompr"),
               "stats400v1r1.tsv", lvls = "just_the_one", concat = F),
    paste("Specified number of outputs is larger than the number ",
          "of outputs in file '",
          system.file("extdata", "nl_ok", package = "micompr"), fs,
          "stats400v1r1.tsv'.", sep = ""),
    fixed = TRUE
  )

  # Should throw error because outputs in files have different lengths
  expect_error(
    grpoutputs(4, c(system.file("extdata", "nl_ok", package = "micompr"),
                    system.file("extdata", "testdata", "n50",
                                package = "micompr")),
               c("stats400v1r1.tsv", "stats400v1r1n50.tsv")),
    paste("Length of outputs in file '",
          system.file("extdata", "testdata", "n50", package = "micompr"),
          fs, "stats400v1r1n50.tsv",
          "' does not match the length of outputs in file '",
          system.file("extdata", "nl_ok", package = "micompr"),
          fs, "stats400v1r1.tsv",
          "'.", sep = ""),
    fixed = TRUE
  )

  # Should expect error because at least 3 outputs are required when requesting
  # output concatenation
  expect_error(
    grpoutputs(2, c(system.file("extdata", "nl_ok", package = "micompr"),
                    system.file("extdata", "j_ex_ok", package = "micompr")),
               c(glob2rx("stats400v1*.tsv"), glob2rx("stats400v1*.tsv")),
               concat = T),
    paste("A minimum of 3 outputs must be specified in order to use ",
          "output concatenation.", sep = ""),
    fixed = TRUE
  )

  # Should expect error because groups have outputs of different lengths
  expect_error(
    grpoutputs(6, c(system.file("extdata", "nl_ok", package = "micompr"),
                    system.file("extdata", "testdata", "NA",
                                package = "micompr")),
               c(glob2rx("stats400v1*.tsv"), glob2rx("stats400v1*A.tsv"))),
    paste("Length of outputs in file '",
          system.file("extdata", "testdata", "NA", package = "micompr"), fs,
          "stats400v1r[0-9]+n20A.tsv' ",
          "does not match the length of outputs in file '",
          system.file("extdata", "nl_ok", package = "micompr"), fs,
          "stats400v1r[0-9]+.tsv'.", sep = "")
  )

})
