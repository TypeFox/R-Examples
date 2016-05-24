

context("value substitution")

A <- matrix(1:12,nrow=4)
b <- rep(1,4)

subst_value(A,b,1,1)

test_that("value substitution",{

  expect_equivalent(
    subst_value(A,b,1,0,remove_columns=TRUE)
    , list(A=A[,2:3],b=b)
  )
  expect_equivalent(
    subst_value(A,b,c(1,3),c(0,1),remove_columns=TRUE)
    ,list(A=A[,2,drop=FALSE],b=b-9:12)
  )
  colnames(A) <- paste0("x",seq_len(ncol(A)))
  expect_equivalent(
    subst_value(A,b,"x1",0,remove_columns=TRUE)
    , list(A=A[,2:3],b=b)
  )
  
})
