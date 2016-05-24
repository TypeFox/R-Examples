library(grpSLOPE)

context("orthogonalizeGroups()")

A <- as.matrix(read.table("./test_data/gaussianMC_test_mat.txt"))

# add an additional group of size 1 to A
A <- cbind(A, 1:nrow(A))

group  <- c(10, "A", 7,10, 6, 9, 9, 2, 5, "0", 6, 2, 6, "0", 2, "0",
            2, "0", 5, 9, 5, "0", 3,10, "A", 5, 3, 7, 6, 9, 8)
grp.id <- getGroupID(group)

test_that("for each group i returns P, Q, R such that A_i[ , P] = Q %*% R", {
  grp.ortho <- orthogonalizeGroups(A, grp.id)
  for (i in 1:length(grp.id)) {
    Ai <- as.matrix(A[ , grp.id[[i]]])
    P  <- grp.ortho[[i]]$P
    Q  <- grp.ortho[[i]]$Q
    R  <- grp.ortho[[i]]$R
    expect_equal(as.matrix(Ai[ , P]), Q %*% R)
  }
})


context("getGroupID()")

test_that("returns a list of vectors with coefficient for each group", {
  expect_equal(grp.id[["10"]], which(group==10)) 
  expect_equal(grp.id[["A"]], which(group=="A")) 
  expect_equal(grp.id[["7"]], which(group==7)) 
  expect_equal(grp.id[["6"]], which(group==6))
  expect_equal(grp.id[["9"]], which(group==9))
  expect_equal(grp.id[["2"]], which(group==2))
  expect_equal(grp.id[["5"]], which(group==5))
  expect_equal(grp.id[["0"]], which(group=="0"))
  expect_equal(grp.id[["3"]], which(group==3))
  expect_equal(grp.id[["8"]], which(group==8))
})
