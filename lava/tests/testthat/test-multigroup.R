context("Multiple Group")

test_that("Multiple group I", {
  m <- lvm(y~x)
  set.seed(1)
  d <- sim(m,100)
  ## Just a stratified analysis
  e <- estimate(list("Group A"=m,"Group B"=m),list(d,d))
  expect_equivalent(coef(e)[c(1,3)],coef(lm(y~x,d)))
  expect_equivalent(coef(e)[c(2,5)],coef(lm(y~x,d)))
})

test_that("Multiple group II", {
  m <- baptize(lvm(y~x))
  set.seed(1)
  d <- sim(m,100)
  ## Just a standard linear regression (single group)
  e <- estimate(list(m,m),list(d,d))
  expect_identical(coef(e,level=2)[[1]],coef(e,level=2)[[2]])
  expect_equivalent(coef(e,level=2)[[1]][1:2,1],coef(lm(y~x,cbind(d,d)))) 
})

context("Missing data")

test_that("Missing data analysis", {
  ## Random intercept model
  m <- lvm(c(y1,y2,y3)~x+u); latent(m) <- ~u
  set.seed(1)
  ## Missing on first two outcomes
  d <- makemissing(sim(m,200),p=0.3,cols=c("y1","y2"))  
  e <- estimate(m,d,missing=TRUE)
  expect_true("lvm.missing"%in%class(e))
  expect_true(sum(unlist(lapply(e$estimate$model$data,nrow)))==200)
  ## Convergence:
  g <- gof(e)
  expect_true(mean(score(e))<1e-3)
  expect_true(g$rankV==length(pars(e)))
})

test_that("Multiple group, missing data analysis", {
  m <- lvm(list(c(y1,y2,y3)~u,u~x)); latent(m) <- ~u
  m <- baptize(fixsome(m))
  regression(m,u~x) <- NA
  covariance(m,~u) <- NA
  set.seed(1)
  ## Missing on all outcomes
  d1 <- makemissing(sim(m,500),cols=c("y1","y2"),p=0.3)
  d2 <- makemissing(sim(m,500),cols=c("y1","y2"),p=0.3)
  e <- estimate(list(m,m),list(d1,d2),missing=TRUE)
  g <- gof(e)
  expect_true(g$n==1000)
  expect_true(mean(score(e))<1e-3)
  expect_true(g$rankV==length(pars(e)))
})


test_that("Multiple group, constraints", {
    m1 <- lvm(y ~ f(x,beta)+f(z,beta2))
    m2 <- lvm(y ~ f(x,psi) + z)
    ## And simulate data from them
    set.seed(1)
    d1 <- sim(m1,100)
    d2 <- sim(m2,100)
    ## Add 'non'-linear parameter constraint
    constrain(m2,psi ~ beta2) <- function(x) x
    ## Add parameter beta2 to model 2, now beta2 exists in both models
    parameter(m2) <- ~ beta2
    ee <- estimate(list(m1,m2),list(d1,d2))

    m <- lvm(y1 ~ x1 + beta2*z1)
    regression(m) <- y2 ~ beta2*x2 + z2
    d <- cbind(d1,d2); names(d) <- c(paste0(names(d1),1),paste0(names(d1),2))
    e <- estimate(m,d)

    b1 <- coef(e,2)["beta2",1]
    b2 <- constraints(ee)[1]
    expect_true(mean((b1-b2)^2)<1e-5)
    
    ## "Multiple group, constraints (non-linear in x)
    m <- lvm(y[m:v] ~ 1)
    addvar(m) <- ~x
    parameter(m) <- ~a+b
    constrain(m,m~a+b+x) <- function(z) z[1]+z[2]*z[3]
    ee <- estimate(list(m,m),list(d1[1:5,],d1[6:10,]))
    b1 <- coef(lm(y~x,d1[1:10,]))
    b2 <- coef(ee)[c("1@a","1@b")]
    expect_true(mean(b1-b2)^2<1e-4)

})






