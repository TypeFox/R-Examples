context("Score matrix functions")

testneurons <- readRDS('testdata/testneurons.rds')

dense_smat <- nblast(testneurons, testneurons)

test_that("can generate non square score matrix from square score matrix", {
  nn=colnames(dense_smat)

  expect_equal(sub_score_mat(nn, nn, scoremat = dense_smat), dense_smat)
  expect_equal(q234t12<-sub_score_mat(query = nn[2:4], target = nn[1:2],
                             scoremat = dense_smat),
               dense_smat[nn[1:2],nn[2:4]])
  expect_equal(q231t21<-sub_score_mat(query = nn[c(2,3,1)], target = nn[2:1],
                             scoremat = dense_smat),
               dense_smat[nn[2:1],nn[c(2,3,1)]])

  # check that we get equivalent answers for all normalisation types
  for(normtype in c("raw", "normalised", "mean")){
    for(distarg in c(FALSE,TRUE)){
      # raw scores are always similarity scores
      if(normtype=='raw' && distarg) next
      q=nn[c(2,4,3)]
      t=nn[2:1]
      baseline=sub_score_mat(scoremat = dense_smat, normalisation = normtype,
                             distance = distarg)[t, q]
      expect_equivalent(sub_score_mat(query = q, target = t, scoremat=dense_smat,
                                 normalisation = normtype, distance = distarg),
                   baseline, info=paste('norm:',normtype,', distance:', distarg))
    }
  }
})


context("Sparse score matrix functions")

sparse_smat <- sparse_score_mat(names(testneurons)[2:4], dense_smat)

test_that("conversion to sparse matrix representation is correct", {
  expect_equal(0, sparse_smat[1, 5])
  expect_equal(diagonal(dense_smat), diagonal(sparse_smat))
})

test_that("dimnames for sparse matrix are correct", {
  expect_equal(dimnames(sparse_smat), dimnames(dense_smat))
})

test_that("subsetting sparse matrix by characters works", {
  expect_equal(sparse_smat["5HT1bMARCM-F000001_seg001", "5HT1bMARCM-F000002_seg001"], sparse_smat[1, 2])
  expect_equal(sparse_smat[, "5HT1bMARCM-F000002_seg001"], sparse_smat[, 2])
  # spam's dropping behaviour is suboptimal...
  expect_equal(sparse_smat["5HT1bMARCM-F000001_seg001", ], sparse_smat[1, , drop=TRUE])
})
