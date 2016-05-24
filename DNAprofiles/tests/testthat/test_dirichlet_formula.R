library(DNAprofiles)
context("Dirichlet formula (pr.next.alleles)")

test_that(desc = "Dirichlet formula gives correct results",{    
  data(freqsNLngm)
  fr0 <- as.vector(freqsNLngm$FGA)
  
  # when theta is 0, sampling should be independent of what is seen
  n <- length(fr0)
  seen <- matrix(sample.int(n = n,size = n^2,replace = TRUE),nrow = n) # bogus
  expect_true(all.equal(pr.next.allele(i = seq(n),seen = seen,f = fr0),fr0))
  
  theta <- 0.3
  expect_true(all.equal(pr.next.allele((c(1,2)),seen=matrix(c(1,1,1,2,2,1),nrow=2,byrow=TRUE),f=c(1/4,3/4),theta=theta),
            c((3*theta+(1-theta)*(1/4))/(1+2*theta),(2*theta+(1-theta)*(3/4))/(1+2*theta))))
  
  # check pr.next.alleles with theta=0
  seen <- matrix(sample.int(n = n,size = n^2,replace = TRUE),nrow = n) # bogus
  ij <- matrix(rep(seq(n),3),ncol=3)
  
  expect_true(all.equal(pr.next.alleles(ij = ij,seen = seen,fr = fr0),fr0^3))  
  
  # check pr.next.alleles with pr.next.allele
  ij <-  matrix(sample.int(n = n,size=2*n,replace = TRUE),ncol=2)
  theta <-0.342
  r1 <- pr.next.allele(i = ij[,2],seen = seen,fr = fr0,theta = theta)*
    pr.next.allele(i = ij[,1],seen = cbind(ij[,2],seen),fr = fr0,theta = theta)
  r2 <- pr.next.alleles(ij = ij,seen = seen,fr = fr0,theta = theta)
  expect_true(all.equal(r1,r2))  
})


