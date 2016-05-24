context("Testing pdbsplit()")


test_that("pdbsplit works", {
  skip_on_cran()
  
  path <- tempdir()
  invisible(capture.output(rawfiles <- get.pdb("3R1C", path=path)))
  invisible(capture.output(files <- pdbsplit(rawfiles, ids=NULL, path=path)))
  
  expected <- c('3R1C_A.pdb', '3R1C_B.pdb', '3R1C_C.pdb', '3R1C_D.pdb',
                '3R1C_E.pdb', '3R1C_F.pdb', '3R1C_G.pdb', '3R1C_H.pdb',
                '3R1C_I.pdb', '3R1C_J.pdb', '3R1C_K.pdb', '3R1C_L.pdb',
                '3R1C_M.pdb', '3R1C_N.pdb', '3R1C_O.pdb', '3R1C_P.pdb',
                '3R1C_Q.pdb', '3R1C_R.pdb', '3R1C_S.pdb', '3R1C_Y.pdb',
                '3R1C_T.pdb', '3R1C_U.pdb', '3R1C_W.pdb', '3R1C_X.pdb',
                '3R1C_V.pdb', '3R1C_Z.pdb', '3R1C_a.pdb', '3R1C_b.pdb',
                '3R1C_c.pdb', '3R1C_d.pdb', '3R1C_e.pdb', '3R1C_f.pdb',
                '3R1C_g.pdb', '3R1C_h.pdb', '3R1C_i.pdb', '3R1C_j.pdb')
  expect_that(expected, equals(basename(files)))
  
  ids <- c('3R1C')
  invisible(capture.output(files <- pdbsplit(rawfiles, ids=ids , path=path)))
  expect_that(expected, equals(basename(files)))
  
  ids <- c('3R1C_e', '3R1C_E')
  invisible(capture.output(files <- pdbsplit(rawfiles, ids=ids , path=path)))
  expected <- c("3R1C_e.pdb", "3R1C_E.pdb")
  expect_that(expected, equals(basename(files)))

  ids <- c('3R1C_XX')
  invisible(capture.output(files <- pdbsplit(rawfiles, ids=ids , path=path)))
  expected <- NULL
  expect_that(expected, equals(files))

  ## multi=TRUE
  invisible(capture.output(rawfiles <- get.pdb("1UD7", path=path)))
  invisible(capture.output(files <- pdbsplit(rawfiles, ids=NULL, path=path, multi=TRUE)))
  expected <- c('1UD7_A.01.pdb', '1UD7_A.02.pdb', '1UD7_A.03.pdb', '1UD7_A.04.pdb',
                '1UD7_A.05.pdb', '1UD7_A.06.pdb', '1UD7_A.07.pdb', '1UD7_A.08.pdb',
                '1UD7_A.09.pdb', '1UD7_A.10.pdb', '1UD7_A.11.pdb', '1UD7_A.12.pdb',
                '1UD7_A.13.pdb', '1UD7_A.14.pdb', '1UD7_A.15.pdb', '1UD7_A.16.pdb',
                '1UD7_A.17.pdb', '1UD7_A.18.pdb', '1UD7_A.19.pdb', '1UD7_A.20.pdb')
  expect_that(expected, equals(basename(files)))
  
  
  pdb <- read.pdb(files[1])
  expect_that(nrow(pdb$atom), equals(1230))


  ## non standard amino acids:
  invisible(capture.output(rawfiles <- get.pdb("1cdk", path=path)))
  invisible(capture.output(files <- pdbsplit(rawfiles, path=path)))
  invisible(capture.output(pdb <- read.pdb(files[1])))
  invisible(capture.output(inds <- atom.select(pdb, resno=197)$atom))

  inds.expected <- c(1568, 1569, 1570, 1571, 1572,
                     1573, 1574, 1575, 1576, 1577, 1578)
  expect_that(inds, equals(inds.expected))
  
}
          )
