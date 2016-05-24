test_that("estimateGraph works for different n.tot and d for all methods",{
  set.seed(1)
  expect_equivalent(c(0, 1.02621411114, 0), 
      estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), 
                             n.tot=100, method="LiuOwen")$tii[,1])
  fun <- function(x) x[,1]*x[,2]

  set.seed(1)
  expect_equivalent(0.006885859902, 
                    estimateGraph(f.mat=fun, d=2, n.tot=10000, method="LiuOwen")$tii[,1])
  
  set.seed(1)
  expect_equivalent(3.291240347 , 
                    estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), 
                                  n.tot=1000, method="FixFast")$tii[2,1])
  set.seed(1)
  expect_equivalent(0.007030367989 , 
                    estimateGraph(f.mat=fun, d=2, n.tot=10000, method="FixFast")$tii[,1])
  
  set.seed(1)
  expect_equivalent(2.621781995673, 
                    estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), 
                                  n.tot=100, method="PickFreeze")$tii[2,1])
  
  set.seed(1)
  expect_equivalent(-60.94382090, 
                    estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), 
                                  n.tot=100, method="RBD")$tii[2,1])
  set.seed(1)
  expect_equivalent(0.008837557797, 
                    estimateGraph(f.mat=fun, d=2, n.tot=10000, method="RBD")$tii[,1])
})

test_that("estimateGraph works for very small values",{
  fun <- function(x) ishigami.fun(x)/1000
  set.seed(1)
  expect_equivalent(0.00000102621411114, 
                    estimateGraph(f.mat=fun, d=3, q.arg=list(min=-pi,max=pi), 
                                  n.tot=100, method="LiuOwen")$tii[2,1])
})

test_that("totalIndex works",{
  set.seed(1)
  t1 <- totalIndex(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), n.mc=1000)
  expect_equivalent(7.129535241, t1[1])
  set.seed(1)
  t2 <- totalIndex(f.mat=sobol.fun, d=8, q.arg=list(min=0,max=1), n.mc=1000)
  expect_equivalent(0.3707765176, t2[1])
})

test_that("confint works for estimateGraph",{
  set.seed(1)
  g <- estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), 
                n.tot=1000, method="LiuOwen")
  expect_equal(c(3,4), dim(g$tii))
  expect_equivalent(0.8796950352725, g$tii[2,2])
})

test_that("overall variance and standard sobol indices work",{
  set.seed(1)
  g <- estimateGraph(f.mat=ishigami.fun, d=3, q.arg=list(min=-pi,max=pi), 
                     n.tot=100, method="LiuOwen")
  expect_equivalent(13.68750659404, g$V)        
  expect_equivalent(6.488175898052, g$i1[2,1])  
  
  set.seed(1)
  g.norm <- estimateGraph(f.mat=ishigami.fun, d=3, q="qnorm", q.arg=list(mean=0, sd=2), 
                               n.tot=100, method="LiuOwen")
  expect_equivalent(107.8001667819, g.norm$V)          
  expect_equivalent(6.817385649107, g.norm$i1[2,1])   
})

test_that("thresholdIdentification works", {
  set.seed(1)
  ############ simple 3-dimensional example with one interaction
  
  ### data (usually existing)
  x <- matrix(seq(0,1,,20), 20, 3)
  x <- apply(x,2,sample)
  y <- 2*(x[,1]-0.5) * (x[,2]-0.5) + 0.1*sin(10*x[,3])
  
  ### FANVOA graph (usually estimated from a meta model over the data)
  g <- list(d=3, 
            tii = matrix(c(0.0140, 0.0008, 0.0002)),
            V = 0.0222, 
            tii.scaled = matrix(c(0.6976, 0.0432, 0.0113))
  )
  class(g) <- "graphlist"
  

  ### Compare candidate thresholds on prediction performance
  comparison <- thresholdIdentification(g, x, y, n.cand = 1)
  expect_equivalent(c(0.328177717247, 0.218761863491, 0.990795265140), comparison$RMSE)
})
