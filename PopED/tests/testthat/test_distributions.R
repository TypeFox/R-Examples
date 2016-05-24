context("Distributions")

test_that("pargen works", {
  
  set.seed(1234)
  
  source("examples_fcn_doc/warfarin_basic.R")
  source("examples_fcn_doc/examples_pargen.R")
  
  # with ln-normal distributions
  #   colMeans(pars.ln)
  #   var(pars.ln)
  mean.diff.ln <- (colMeans(pars.ln) - bpop_vals_ed_ln[,2])/bpop_vals_ed_ln[,2]*100
  var.diff.ln <- (diag(var(pars.ln)) - bpop_vals_ed_ln[,3])/bpop_vals_ed_ln[,3]*100
  
  expect_that(all(mean.diff.ln<20),is_true())
  expect_that(all(var.diff.ln[1:3]<50),is_true())
  
  set.seed(1234)
  pars.ln.1 <- pargen(par=bpop_vals_ed_ln,
                    user_dist_pointer=NULL,
                    sample_size=1,
                    bLHS=1,
                    sample_number=NULL,
                    poped.db)
  
  set.seed(1234)
  pars.ln.2 <- pargen(par=bpop_vals_ed_ln,
                      user_dist_pointer=NULL,
                      sample_size=1,
                      bLHS=1,
                      sample_number=NULL,
                      poped.db)
  
  expect_true(all(pars.ln.1==pars.ln.2))
  
  p.vals <- apply(pars.ln[,1:3],2,function(x) shapiro.test(log(x))[["p.value"]])
  expect_that(all(p.vals > 0.05), is_true())
  
  # Adding 10% Uncertainty to fixed effects normal-distribution (not Favail)
  # with normal distributions
  #   # Looks ok
  #   colMeans(pars.n)
  #   var(pars.n)
  #   
  mean.diff.n <- (colMeans(pars.n) - bpop_vals_ed_n[,2])/bpop_vals_ed_n[,2]*100
  var.diff.n <- (diag(var(pars.n)) - bpop_vals_ed_n[,3])/bpop_vals_ed_n[,3]*100
  
  expect_that(all(mean.diff.n<20),is_true())
  expect_that(all(var.diff.n[1:3]<50),is_true())
  
  p.vals <- apply(pars.n[,1:3],2,function(x) shapiro.test(x)[["p.value"]])
  expect_that(all(p.vals > 0.05), is_true())
  
  # Adding 10% Uncertainty to fixed effects uniform-distribution (not Favail)
  mean.diff.u <- (colMeans(pars.u) - bpop_vals_ed_u[,2])/bpop_vals_ed_u[,2]*100
  range.diff.u <- mean.diff.u
  for(i in 1:4){
    range.diff.u[i] <- (diff(range(pars.u[,i])) - bpop_vals_ed_u[i,3])/bpop_vals_ed_u[i,3]*100
  }
  
  expect_that(all(mean.diff.u<20),is_true())
  expect_that(all(range.diff.u[1:3]<50),is_true())

  # Adding user defined distributions
  mean.diff.ud <- (colMeans(pars.ud) - bpop_vals_ed_ud[,2])/bpop_vals_ed_ud[,2]*100
  sd.diff.ud <- (sqrt(diag(var(pars.ud))) - bpop_vals_ed_ud[,3])/bpop_vals_ed_ud[,3]*100
  
  expect_that(all(mean.diff.ud<20),is_true())
  expect_that(all(sd.diff.ud[1:3]<50),is_true())
  
    
})

test_that("log_prior_pdf works", {
  
  source("examples_fcn_doc/warfarin_optimize.R")
  source("examples_fcn_doc/examples_log_prior_pdf.R")
  
  # then compute the log density from log_prior_pdf
  alpha <- bpop_vals_ed_n[bpop_vals_ed_n[,1]!=0,2]
  lp.1 <- log_prior_pdf(alpha=alpha, bpopdescr=bpop_vals_ed_n, ddescr=poped.db$parameters$d)
  
  # compute the log density from dnorm
  par <- bpop_vals_ed_n[bpop_vals_ed_n[,1]==1,]
  lp.2 <- log(prod(dnorm(par[,2],par[,2],sqrt(par[,3]))))
  
  test_that(lp.1,equals(lp.2))
  
  # then compute the log density 
  alpha <- bpop_vals_ed_ln[bpop_vals_ed_ln[,1]!=0,2]
  lp.3 <- log_prior_pdf(alpha=alpha, bpopdescr=bpop_vals_ed_ln, ddescr=poped.db$parameters$d)
  
  # compute the log density from dlnorm
  par <- bpop_vals_ed_ln[bpop_vals_ed_ln[,1]==4,]
  lp.4 <- log(prod(dlnorm(par[,2],log(par[,2]^2/sqrt(par[,3]+par[,2]^2)),sqrt(log(par[,3]/par[,2]^2+1)))))
  
  test_that(lp.3,equals(lp.4))
  
  # then compute the log density 
  alpha <- bpop_vals_ed_ln[bpop_vals_ed_u[,1]!=0,2]
  lp.5 <- log_prior_pdf(alpha=alpha, bpopdescr=bpop_vals_ed_u, ddescr=poped.db$parameters$d)
  
  # compute the log density from dlnorm
  par <- bpop_vals_ed_u[bpop_vals_ed_u[,1]==2,]
  lp.6 <- log(prod(dunif(par[,2],min=par[,2]-par[,3]/2,max=par[,2]+par[,3]/2)))
  
  test_that(lp.5,equals(lp.6))
  
  
})
