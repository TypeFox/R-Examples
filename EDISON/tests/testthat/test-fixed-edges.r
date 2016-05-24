context('fixed-moves')

dataset = simulateNetwork(l=30, cps=c(10,20))

num.iter = 500

sim.data = dataset$sim_data
options = defaultOptions()
options$cp.fixed = TRUE
options$cp.init = list()
options$cp.init[[1]] = c(2,7,17,dim(sim.data)[2])
options$cp.init[[2]] = c(2,13,23,dim(sim.data)[2])

for(i in 3:dim(sim.data)[1]) {
  options$cp.init[[i]] = c(2,dim(sim.data)[2])
}

fixed.edges = matrix(-1, dim(sim.data)[1], dim(sim.data)[1])
fixed.edges[1,4] = 1
fixed.edges[7,5] = 0

result = EDISON.run(sim.data, num.iter=num.iter, options=options, 
                    fixed.edges=fixed.edges)
network = calculateEdgeProbabilities(result, cps=c(2, dataset$epsilon))

# Test fixed edges remain fixed (no changepoints) 
test_that('Fixed works without changepoints',
          {expect_that(network$probs.all[[1]][1,4], equals(1))
          expect_that(network$probs.all[[1]][7,5], equals(0))
          expect_that(network$probs.seg[[1]][1,4], equals(1))
          expect_that(network$probs.seg[[1]][7,5], equals(0))
          expect_that(network$probs.seg[[2]][1,4], equals(1))
          expect_that(network$probs.seg[[2]][7,5], equals(0))
          expect_that(network$probs.seg[[3]][1,4], equals(1))
          expect_that(network$probs.seg[[3]][7,5], equals(0))})

options$cp.fixed = FALSE
result = EDISON.run(sim.data, num.iter=num.iter, fixed.edges=fixed.edges, 
                    options=options)
network = calculateEdgeProbabilities(result, cps=c(2, dataset$epsilon))


test_that('Fixed works with changepoints',
          {expect_that(network$probs.all[[1]][1,4], equals(1))
          expect_that(network$probs.all[[1]][7,5], equals(0))
          expect_that(network$probs.seg[[1]][1,4], equals(1))
          expect_that(network$probs.seg[[1]][7,5], equals(0))
          expect_that(network$probs.seg[[2]][1,4], equals(1))
          expect_that(network$probs.seg[[2]][7,5], equals(0))
          expect_that(network$probs.seg[[3]][1,4], equals(1))
          expect_that(network$probs.seg[[3]][7,5], equals(0))})