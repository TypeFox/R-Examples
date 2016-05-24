require(distrSim)
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


