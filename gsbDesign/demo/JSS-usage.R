## Illustrates the usage of gsbDesign.
## See Section 4 of the vignette JSS-gsbDesign.pdf
## for more information.

## spacify trial design
design1 <- gsbDesign(nr.stages = 4, patients = c(10, 20), 
  sigma = c(7, 7), criteria.success = c(0, 0.8, 7, 0.5), 
  criteria.futility = c(2, 0.8), prior.difference = c(3, 5, 2))
names(design1)

design2 <- gsbDesign(nr.stages = 4, patients = c(10, 20), sigma = c(7, 7),
  criteria.success = c(0, 0.8, 7, 0.5), criteria.futility = c(2, 0.8),
  prior.control = c(3, 5), prior.treatment = c(6, 2))

## specify calculation / simulation setup
simulation1 <- gsbSimulation(truth = c(-10, 20, 60), 
  type.update = "treatment effect", method = "numerical integration")

simulation1a <- gsbSimulation(truth = c(-10, 20, 60), 
  type.update = "treatment effect", method = "simulation", 
  nr.sim = 50000, warnings.sensitivity = 100, seed = "generate")

simulation2 <- gsbSimulation(truth = list(seq(-5, 5, 3), seq(0, 5, 3)),
  type.update = "per arm", method = "simulation", grid.type = "table",
  nr.sim = 10000, warnings.sensitivity = 500, seed = "generate")

## calculate operating characteristics
oc1 <- gsb(design1, simulation1)
names(oc1)

summ.oc1 <- summary(oc1, atDelta = c(0, 2, 7))
summary(oc1, atDelta = c(0, 2, 7))

p1 <- plot(oc1, what = "cumulative all")
print(p1)
formals(gsbDesign:::plot.gsbMainOut)$what

p2 <- plot(oc1, what = "sample size")
p3 <- plot(oc1, what = "boundary")

print(p2, position = c(0,0, 0.5, 1), more = TRUE)
print(p3, position = c(0.5, 0,1, 1), more = FALSE)

tab(oc1, what = "cumulative success", atDelta = c(0, 2, 7), digits = 4, 
  export = FALSE)

## calculate operating characteristics
oc2 <- gsb(design2, simulation2)

tab(oc2, what = "sample size", digits = 0)

simulation2a <- gsbSimulation(truth = c(-5, 5, 0, 5, 50), 
  type.update = "per arm", method = "simulation", grid.type = "plot",
  nr.sim = 100000, warnings.sensitivity = 50, seed = "generate")
oc2a <- gsb(design2, simulation2a)

p12 <- plot(oc2a, what = "cumulative all", delta.grid = FALSE)
print(p12)
