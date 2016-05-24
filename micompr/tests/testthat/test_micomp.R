# Copyright (c) 2016 Nuno Fachada
# Distributed under the MIT License (http://opensource.org/licenses/MIT)

library(micompr)
context("micomp")

test_that("micomp constructs the expected objects", {

  # Output names
  outputs <- c("PopSheep", "PopWolf", "QtyGrass",
               "EnSheep", "EnWolf", "EnGrass",
               "All")

  # Minimum percentage of variance to be explained
  ve_npcs <- c(0.5, 2, 0.9, 4)

  # Which of these are variances?
  idxvar <- which(ve_npcs < 1)
  # And which are PCs?
  idxpcs <- which(ve_npcs > 1)

  # Determine location of extdata files
  dir_nl_ok <- system.file("extdata", "nl_ok", package = "micompr")
  dir_jex_ok <- system.file("extdata", "j_ex_ok", package = "micompr")
  dir_jex_noshuff <- system.file("extdata", "j_ex_noshuff", package = "micompr")
  dir_jex_diff <- system.file("extdata", "j_ex_diff", package = "micompr")
  dir_na <- system.file("extdata", "testdata", "NA", package = "micompr")
  files <- glob2rx("stats400v1*.tsv")
  filesA_na <- glob2rx("stats400v1*n20A.tsv")
  filesB_na <- glob2rx("stats400v1*n20B.tsv")

  # 1 - Build a micomp object using data from extdata files

  # 1a - Use files containing package datasets, three comparisons
  mic1a <- micomp(outputs, ve_npcs,
                  list(
                    list(name = "NLOKvsJEXOK",
                         folders = c(dir_nl_ok, dir_jex_ok),
                         files = c(files, files),
                         lvls = c("NLOK", "JEXOK")),
                    list(name = "NLOKvsJEXNOSHUFF",
                         folders = c(dir_nl_ok, dir_jex_noshuff),
                         files = c(files, files),
                         lvls = c("NLOK", "JEXNOSHUFF")),
                    list(name = "NLOKvsJEXDIFF",
                         folders = c(dir_nl_ok, dir_jex_diff),
                         files = c(files, files),
                         lvls = c("NLOK", "JEXDIFF"))),
                  concat = T)

  # 1b - Use files containing test dataset, one comparison, just five outputs
  # (unnamed), no concatenation, unnamed levels
  mic1b <- micomp(5, ve_npcs,
                  list(
                    list(name = "testVLOdata",
                         folders = dir_na,
                         files = c(filesA_na, filesB_na))),
                  mnv_test = "Roy")

  # 2 - Use package datasets (i.e. grpoutputs objects) directly
  mic2 <- micomp(outputs, ve_npcs,
                 list(
                   list(name = "NLOKvsJEXOK", grpout = pphpc_ok),
                   list(name = "NLOKvsJEXNOSHUFF", grpout = pphpc_noshuff),
                   list(name = "NLOKvsJEXDIFF", grpout = pphpc_diff)),
                 concat = T,
                 lim_npcs = T,
                 mnv_test = "Wilks")

  # 3 - Use manually inserted data, unnamed outputs, no concatenation
  mic3 <- micomp(6, ve_npcs,
                 list(
                   list(name = "NLOKvsJEXOK",
                        grpout = list(data = pphpc_ok$data,
                                      obs_lvls = pphpc_ok$obs_lvls)),
                   list(name = "NLOKvsJEXNOSHUFF",
                        grpout = list(data = pphpc_noshuff$data,
                                      obs_lvls = pphpc_noshuff$obs_lvls)),

                   list(name = "NLOKvsJEXDIFF",
                        grpout = list(data = pphpc_diff$data,
                                      obs_lvls = pphpc_diff$obs_lvls))),
                 concat = F,
                 mnv_test = "Hotelling-Lawley")

  ##### Start testing #####

  # Check object dimensions
  expect_equal(dim(mic1a), c(7, 3))
  expect_equal(dim(mic1b), c(5, 1))
  expect_equal(dim(mic2), c(7, 3))
  expect_equal(dim(mic3), c(6, 3))

  # Check object row names
  expect_equal(rownames(mic1a), outputs)
  expect_equal(rownames(mic1b), c("out1", "out2", "out3", "out4", "out5"))
  expect_equal(rownames(mic2), outputs)
  expect_equal(rownames(mic3),
               c("out1", "out2", "out3", "out4", "out5", "out6"))

  # Check object column names
  expect_equal(colnames(mic1a),
               c("NLOKvsJEXOK", "NLOKvsJEXNOSHUFF", "NLOKvsJEXDIFF"))
  expect_equal(colnames(mic1b),
               "testVLOdata")
  expect_equal(colnames(mic2),
               c("NLOKvsJEXOK", "NLOKvsJEXNOSHUFF", "NLOKvsJEXDIFF"))
  expect_equal(colnames(mic3),
               c("NLOKvsJEXOK", "NLOKvsJEXNOSHUFF", "NLOKvsJEXDIFF"))

  # Check properties of sub-objects
  for (i in 1:dim(mic1a)[1]) {
    for (j in 1:dim(mic1a)[2]) {

      # Get current subobject
      sobj <- mic1a[[i, j]]

      # Is subobject a cmpoutput object?
      expect_is(sobj, "cmpoutput")

      # Check that the number of PCs which explain the specified minimum
      # percentage of variance has the expected value

      # Case 1 - Variance to explain was specified
      expect_equal(sobj$npcs[idxvar],
                   sapply(ve_npcs[idxvar],
                          function(mv, ve) match(T, cumsum(ve) > mv),
                          sobj$varexp))
      # Case 2 - Number of PCs was directly specified
      expect_equal(sobj$ve[idxpcs],
                   sapply(ve_npcs[idxpcs],
                          function(npcs, ve) sum(ve[1:npcs]),
                          sobj$varexp))

    }
  }

})

test_that("micomp throws errors when improperly invoked", {

  # Don't specify files in the first list
  expect_error(
    micomp(7, 0.6,
           list(list(name = "A", folders = c("dir1", "dir2")),
                list(name = "B", files = c("file1", "file2")))),
    "Invalid 'comps' parameter",
    fixed = TRUE
  )

  # Don't specify a comparison name in the second list
  expect_error(
    micomp(c("o1", "o2", "o3"), 0.75,
           list(list(name = "aName", grpout = pphpc_ok),
                list(grpout = pphpc_noshuff))),
    "Invalid 'comps' parameter",
    fixed = TRUE
  )

  # Don't specify observation levels in the third list. This will provoke an
  # error in cmpoutput.
  expect_error(
    micomp(6, 0.5,
           list(
             list(name = "NLOKvsJEXOK",
                  grpout = list(data = pphpc_ok$data,
                                obs_lvls = pphpc_ok$obs_lvls)),
             list(name = "NLOKvsJEXNOSHUFF",
                  grpout = list(data = pphpc_noshuff$data,
                                obs_lvls = pphpc_noshuff$obs_lvls)),

             list(name = "NLOKvsJEXDIFF",
                  grpout = list(data = pphpc_diff$data)))),
    "Number of observations in 'data' and 'obs_lvls' does not match.",
    fixed = TRUE

  )

})

test_that("micomp assumptions have the correct properties", {

  # Minimum percentage of variance to be explained
  minvar <- 0.8

  # Determine location of extdata files
  dir_nl_ok <- system.file("extdata", "nl_ok", package = "micompr")
  dir_jex_ok <- system.file("extdata", "j_ex_ok", package = "micompr")
  dir_jex_noshuff <- system.file("extdata", "j_ex_noshuff", package = "micompr")
  dir_jex_diff <- system.file("extdata", "j_ex_diff", package = "micompr")
  dir_na <- system.file("extdata", "testdata", "NA", package = "micompr")
  files <- glob2rx("stats400v1*.tsv")
  filesA_na <- glob2rx("stats400v1*n20A.tsv")
  filesB_na <- glob2rx("stats400v1*n20B.tsv")

  ##### Create micomp objects #####

  # 1 - Build a micomp object using data from extdata files

  # 1a - Use files containing package datasets, three comparisons
  mic1a <- micomp(7, minvar,
                  list(
                    list(name = "NLOKvsJEXOK",
                         folders = c(dir_nl_ok, dir_jex_ok),
                         files = c(files, files),
                         lvls = c("NLOK", "JEXOK")),
                    list(name = "NLOKvsJEXNOSHUFF",
                         folders = c(dir_nl_ok, dir_jex_noshuff),
                         files = c(files, files),
                         lvls = c("NLOK", "JEXNOSHUFF")),
                    list(name = "NLOKvsJEXDIFF",
                         folders = c(dir_nl_ok, dir_jex_diff),
                         files = c(files, files),
                         lvls = c("NLOK", "JEXDIFF"))),
                  concat = T)

  # 1b - Use files containing test dataset, one comparison, just five outputs
  # (unnamed), no concatenation, unnamed levels
  mic1b <- micomp(5, minvar,
                  list(
                    list(name = "testVLOdata",
                         folders = dir_na,
                         files = c(filesA_na, filesB_na))))

  # 2 - Use package datasets (i.e. grpoutputs objects) directly
  mic2 <- micomp(7, minvar,
                 list(
                   list(name = "NLOKvsJEXOK", grpout = pphpc_ok),
                   list(name = "NLOKvsJEXNOSHUFF", grpout = pphpc_noshuff),
                   list(name = "NLOKvsJEXDIFF", grpout = pphpc_diff)),
                 concat = T)

  # 3 - Use manually inserted data, unnamed outputs, no concatenation
  mic3 <- micomp(6, minvar,
                 list(
                   list(name = "NLOKvsJEXOK",
                        grpout = list(data = pphpc_ok$data,
                                      obs_lvls = pphpc_ok$obs_lvls)),
                   list(name = "NLOKvsJEXNOSHUFF",
                        grpout = list(data = pphpc_noshuff$data,
                                      obs_lvls = pphpc_noshuff$obs_lvls)),

                   list(name = "NLOKvsJEXDIFF",
                        grpout = list(data = pphpc_diff$data,
                                      obs_lvls = pphpc_diff$obs_lvls))),
                 concat = F)

  ##### Create an assumptions_micomp object for each micomp object #####

  oldw <- getOption("warn")
  options(warn = -1)

  am1a <- assumptions(mic1a)
  am1b <- assumptions(mic1b)
  am2 <- assumptions(mic2)
  am3 <- assumptions(mic3)

  options(warn = oldw)

  ##### Start testing #####

  # Check that the objects are of the correct type
  expect_is(am1a, "assumptions_micomp")
  expect_is(am1b, "assumptions_micomp")
  expect_is(am2, "assumptions_micomp")
  expect_is(am3, "assumptions_micomp")

  # Check that assumptions objects have the same dimensions as the respective
  # micomp objects
  expect_equal(dim(am1a), dim(mic1a))
  expect_equal(dim(am1b), dim(mic1b))
  expect_equal(dim(am2), dim(mic2))
  expect_equal(dim(am3), dim(mic3))

  # Check that assumptions objects have the same row names as the respective
  # micomp objects
  expect_equal(rownames(am1a), rownames(mic1a))
  expect_equal(rownames(am1b), rownames(mic1b))
  expect_equal(rownames(am2), rownames(mic2))
  expect_equal(rownames(am3), rownames(mic3))

  # Check that assumptions objects have the same column names as the respective
  # micomp objects
  expect_equal(colnames(am1a), colnames(mic1a))
  expect_equal(colnames(am1b), colnames(mic1b))
  expect_equal(colnames(am2), colnames(mic2))
  expect_equal(colnames(am3), colnames(mic3))

  # Check properties of sub-objects
  for (a in list(am1a, am1b, am2, am3)) {
    for (i in 1:(dim(a)[1])) {
      for (j in 1:(dim(a)[2])) {

        # Get current subobject
        sobj <- a[[i, j]]

        # Is subobject a assumptions_cmpoutput object?
        expect_is(sobj, "assumptions_cmpoutput")
      }
    }
  }

  #### Test the summary function
  sam1a <- summary(am1a)
  sam1b <- summary(am1b)
  sam2 <- summary(am2)
  sam3 <- summary(am3)

  expect_equal(names(sam1a),
               c("NLOKvsJEXOK", "NLOKvsJEXNOSHUFF", "NLOKvsJEXDIFF"))
  expect_equal(names(sam1b), "testVLOdata")
  expect_equal(names(sam2),
               c("NLOKvsJEXOK", "NLOKvsJEXNOSHUFF", "NLOKvsJEXDIFF"))
  expect_equal(names(sam3),
               c("NLOKvsJEXOK", "NLOKvsJEXNOSHUFF", "NLOKvsJEXDIFF"))
})