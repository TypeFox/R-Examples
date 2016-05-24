require(distrTEst)
options("newDevice"=TRUE)

sim <- new("Simulation",
           seed = setRNG(),
           distribution = Norm(mean = 0, sd = 1),
           filename="sim_01",
           runs = 1000,
           samplesize = 30)

contsim <- new("Contsimulation",
               seed = setRNG(),
               distribution.id = Norm(mean = 0, sd = 1),
               distribution.c = Norm(mean = 0, sd = 9),
               rate = 0.1,
               filename="contsim_01",
               runs = 1000,
               samplesize = 30)

simulate(sim)
simulate(contsim)

print(sim)
summary(contsim)
plot(contsim)


psim <- function(theta,y,m0){
  mean(pmin(pmax(-m0, y - theta), m0))
}
mestimator <- function(x, m = 0.7) {
  uniroot(f = psim,
          lower = -20,
          upper = 20,
          tol = 1e-10,
          y = x,
          m0 = m,
          maxiter = 20)$root
}

result.id.mean <- evaluate(sim, mean)
result.id.mest <- evaluate(sim, mestimator)
result.id.median <- evaluate(sim, median)


result.cont.mean <- evaluate(contsim, mean)
result.cont.mest <- evaluate(contsim, mestimator)
result.cont.median <- evaluate(contsim, median)

elist <- EvaluationList(result.cont.mean,
                        result.cont.mest,
                        result.cont.median) 

print(elist)
summary(elist)
plot(elist,cex=0.7)
