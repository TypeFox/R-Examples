context("testing prior functions")

test_that("testing prior functions", {
  data(chelonia.simmap)
  tree <- chelonia.simmap$tree
  dat <- chelonia.simmap$dat
  emap <- chelonia.simmap$emap
  cache <- .prepare.ou.univariate(tree, dat)
  pars <- list(alpha=0.1, sig2=1, k=16, theta=c(3,4,5,6), sb=which(emap$sh==1), loc=emap$r1[which(emap$sh==1)])
  QGpars <- list(h2=0.1,P=1,w2=0.9,Ne=1,k=16,theta=c(3,4,5,6), sb=which(emap$sh==1), loc=emap$r1[which(emap$sh==1)])
  prior <- make.prior(tree,dists=list(dalpha="dunif",dsig2="dunif"),param=list(dalpha=list(min=0,max=1),dsig2=list(min=0,max=1),dsb=list(bmax=1,prob=1)),type="pars",plot.prior=TRUE)
  QGprior <- make.prior(tree,dists=list(dh2="dunif",dP="dunif",dw2="dunif",dNe="dunif"),param=list(dh2=list(min=0,max=1),dP=list(min=0,max=1),dw2=list(min=0,max=1),dNe=list(min=0,max=1)),model="QG",type="pars")
  expect_that(class(try(make.prior(tree,type="pars"),silent=TRUE))[1],equals("priorFn"))
  expect_that(class(try(make.prior(tree,type="emap"),silent=TRUE))[1],equals("priorFn"))
  expect_that(class(try(make.prior(tree,type="simmap"),silent=TRUE))[1],equals("priorFn"))
  expect_that(class(try(make.prior(tree,model="QG",type="pars"),silent=TRUE))[1],equals("priorFn"))
  expect_that(class(try(make.prior(tree,model="OUrepar",type="emap"),silent=TRUE))[1],equals("priorFn"))
  expect_that(class(try(make.prior(tree,model="QG",type="simmap"),silent=TRUE))[1],equals("priorFn"))
  expect_that(prior(pars),equals(-145.1539,tolerance=0.0001))
  expect_that(QGprior(QGpars),equals(prior(pars)))
  f1tree <- tree
  f1tree$edge.length[pars$sb] <- tree$edge.length[pars$sb]*100
  priorf1 <- make.prior(tree,dists=list(dalpha="dunif",dsig2="dunif"),param=list(dalpha=list(min=0,max=1),dsig2=list(min=0,max=1),dsb=list(bmax=Inf,prob=tree$edge.length)),type="pars")
  priorf2 <- make.prior(f1tree,dists=list(dalpha="dunif",dsig2="dunif"),param=list(dalpha=list(min=0,max=1),dsig2=list(min=0,max=1),dsb=list(bmax=Inf,prob=f1tree$edge.length)),type="pars")
  expect_that(priorf1(pars)< priorf2(pars), is_true())  
  expect_that(QGprior(pars),throws_error("Missing parameters:  h2 P w2 Ne"))
})
