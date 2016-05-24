context("wilcoxon.R internal methods")

test.mu<- 0.5

test_that("mu.wilcox computes correctly",{
  expect_equal(mu.wilcox(1,1),test.mu)
})


test_that("sd.wilcox computes correctly",{
  expect_equal(sd.wilcox(1,1),sd.wilcox(1,1))
})

p<-1:20
weights<-p/length(p)
snps.in<-seq(2,20,2)

test_that("calc.wilcoxon computes correctly",{
  expect_equal(calc.wilcoxon(p,snps.in,length(snps.in)),55)
  expect_equal(calc.wilcoxon(p,snps.in,length(snps.in),w=1),calc.wilcoxon(p,snps.in,length(snps.in)))
  expect_equal(calc.wilcoxon(p,snps.in,length(snps.in),w=2),165)
})

context("wilcoxon public method")

#in the case of no weight then just a wrapper for calc.wilcoxon


test_that("wilcoxon computes correctly in the context of vector input",{
  expect_equal(wilcoxon(p,snps.in),55)
  ##in turn expect
  expect_equal(wilcoxon(p,snps.in),calc.wilcoxon(p,snps.in,length(snps.in)))
  ## test propensity score stuff
  weights<-p/length(p)
  ## this should throw an error for being too sparse
  expect_that(wilcoxon(p,snps.in,weights,0.1),throws_error(regexp="sparse"))
  #this should be OK
  w<-wilcoxon(p,snps.in,weights,0.5)
  ##should be equivalent to wilcoxon without weighting
  expect_equal(w,wilcoxon(p,snps.in))
})
	

test_that("wilcoxon computes correctly in the context of matrix input",{
  p.matrix<-matrix(1:200,nrow=20,ncol=10)
  p.matrix<-p.matrix[sample.int(nrow(p.matrix)),]
  out<-wilcoxon(p.matrix,snps.in)
  expect_is(out,"numeric")
  expect_equal(length(out),ncol(p.matrix))
  #these should be equivalent
  expect_equal(unique(out),wilcoxon(p.matrix[,1],snps.in))
})
