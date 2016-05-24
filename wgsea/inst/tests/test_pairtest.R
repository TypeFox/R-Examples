library(snpStats)
data(for.exercise, package="snpStats")
test.data <- snps.10[,11:20]

## pairtest.R

context("pairtest.R")

## setup

data.case<-test.data[subject.support$cc==1, ]
data.control<-test.data[subject.support$cc==0, ]
## create a test matrix of phenotypes
total.cc<-sum(nrow(data.case),nrow(data.control))
known.pheno<-matrix(as.integer(0),total.cc,1)
known.pheno[1:nrow(data.case),1]<-1
p.1df<-pairtest(data.case,data.control,quiet=TRUE)
p.1df.perm<-pairtest(data.case,data.control,5,quiet=TRUE)
p.1df.pheno<-pairtest(data.case,data.control,pheno.perm=known.pheno,quiet=TRUE)

test_that("our inputs are valid",{
  expect_equal(colnames(data.control),colnames(data.case))
  expect_that(data.case, is_a("SnpMatrix"))
  expect_that(data.control, is_a("SnpMatrix"))
})


test_that("pairtest throws error when either data.case or data.control param is not a snpMatrix object",{
  expect_that(pairtest(1,data.control),throws_error())
  expect_that(pairtest(data.case,1),throws_error())
})

test_that("pairtest throws error when data.case and data.control snpMatrix objects have different SNPs",{
  expect_that(pairtest(data.case[,1:5],data.control),throws_error())
})


test_that("pairtest return when n.perm=0",{
  ## this is probably enough
  expect_equal(names(p.1df),colnames(data.case))
  ## however we can check that none >1 and <0
  expect_that(unique(p.1df<0),is_false())
  expect_that(unique(p.1df>1),is_false())
})

test_that("pairtest returns when n.perm>0",{
  expect_that(p.1df.perm,is_a("matrix"))
  expect_equal(nrow(p.1df.perm),ncol(data.case))
  expect_equal(ncol(p.1df.perm),5)
})


test_that("pairtest throws error when pheno.perm does not match total sample number",{
  expect_that(pairtest(data.case,data.control,pheno.perm=known.pheno[-5,]),throws_error())
})

test_that("pairtest pheno.perm parameter returns correct values",{
  expect_that(p.1df.pheno,is_a("matrix"))
  names(p.1df)<-NULL
  expect_identical(p.1df,p.1df.pheno[,1])
})

