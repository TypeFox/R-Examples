context("Finding Feasible Matrices")

set.seed(123787)


test_that("Basic Checks",{
  skip_on_cran()


  p <- matrix(c(1,0,0,1),nrow=2)
  expect_error(findFeasibleMatrix(c(0,0,0),c(0,0),p))  # wrong number of rows
  expect_error(findFeasibleMatrix(c(0,0),c(0,1,2),p))  # wrong number of columns
  expect_error(findFeasibleMatrix(c(0,-1),c(0,0),p))  # negative row
  expect_error(findFeasibleMatrix(c(0,0),c(0,-2),p))  # negative column
  expect_error(findFeasibleMatrix(c(0,1),c(2,0),p),
               "Sums of r and c differ")  # row and column sum do not match

  #Violation of existence condition.
  expect_error(findFeasibleMatrix(c(0,1),c(1,0),matrix(nrow=2,ncol=2,0))) ## all p=0
  expect_error(findFeasibleMatrix(c(1,2),c(2,1),p),
               "row\\(s\\) 2 and column\\(s\\) 1")
  p2 <- diag(x=1,nrow=3,ncol=4)
  expect_error(findFeasibleMatrix(c(1,2,3),c(3,2,1,0),p2))
  expect_error(findFeasibleMatrix(c(1,2,3),c(3,2,1),matrix(c(1,1,0,1,1,1,0,0,1),byrow=TRUE,nrow=3)),
               "row\\(s\\) 3 and column\\(s\\) 1 2")
  expect_error(findFeasibleMatrix(c(1,2,3),c(3,2,1),matrix(c(1,1,0,0,1,1,0,1,1),byrow=TRUE,nrow=3)),
               "row\\(s\\) 2 3 and column\\(s\\) 1")
  expect_error(findFeasibleMatrix(c(1,2,4),c(1,2,4),matrix(c(0,1,1,1,0,1,1,1,0),byrow=TRUE,nrow=3)),
               "row\\(s\\) 3 and column\\(s\\) 3")


  expect_error(findFeasibleMatrix(c(0,1),c(2,0),matrix(nrow=2,ncol=2,0)))  # row and column sums different


  expect_equal(findFeasibleMatrix(c(0,0),c(0,0),p),matrix(c(0,0,0,0),nrow=2))
  expect_equal(findFeasibleMatrix(c(1,1),c(1,1),p),matrix(c(1,0,0,1),nrow=2))



})

test_that("Find feasible matrix with desired average degree",{
    set.seed(12180980)
    Lorig <- matrix(rbinom(200,1,0.5)*abs(rcauchy(200)),nrow=20,ncol=10)
    p <- pmin(1+Lorig,0.5)
    L <- findFeasibleMatrix_targetmean(rowSums(Lorig),colSums(Lorig),p)
    expect_true(all(L>=0))
    expect_equal(rowSums(L),rowSums(Lorig))
    expect_equal(colSums(L),colSums(Lorig))
    expect_warning(findFeasibleMatrix_targetmean(rowSums(Lorig),colSums(Lorig),p,targetmean=0.01))
})

test_that("Checks of Finding a Feasible Matrix for Randomly Generated Matrices",{
  skip_on_cran()

  for (alpha in seq(0,1,by=0.25)){
      for (n in 1:30){
          for (nr in c(5,n)){
              M <- matrix(nrow=nr,ncol=n,rexp(nr*n)*(runif(nr*n)>alpha))
              r <- rowSums(M)
              c <- colSums(M)
              Mnew <- findFeasibleMatrix(r=r,c=c,p=(M>0)*0.5)
              expect_equal(rowSums(Mnew),r)
              expect_equal(colSums(Mnew),c)
              expect_equal(all(!((M==0)&(Mnew>0))),TRUE) ## check that all 0s are present
          }
      }
  }


})

