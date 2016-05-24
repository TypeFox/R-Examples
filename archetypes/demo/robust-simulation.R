### Weighted and robust archetypal analysis: simulation study
###
### Analysis used in 'Weighted and Robust Archetypal Analysis' by
### Manuel J. A. Eugster and Friedrich Leisch.

library('archetypes')
library('ggplot2')


sim.file <- function(x) {
  system.file('opt', 'robust-simulation',
              sprintf('%s.Rdata', x),
              package = 'archetypes')
}


### Simulation 1: ####################################################

load(sim.file("sim1"))
str(sim1)


### Distances for one dimension and one number of outliers:

p1 <- subset(sim1, dim == 10 & perf %in% c("dist1", "dist2") & nout == 100)
p1 <- cast(p1, dim + n + nout + radius + sample + alg ~ perf)

ggplot(p1, aes(dist2, dist1)) +
  geom_point() +
  facet_grid(alg ~ radius)



### Distances for over the dimensions:

p2 <- subset(sim1, perf %in% c("dist1", "dist2") & nout == 100 & radius == 15)
p2 <- cast(p2, dim + n + nout + radius + sample + alg ~ perf)

ggplot(p2, aes(dist2, dist1)) +
  geom_point() +
  facet_grid(alg ~ dim)



### Number of iterations:

p3 <- subset(sim1, perf == "iters")
p3 <- ddply(p3, c("dim", "n", "nout", "radius", "perf", "alg"),
            function(x)
            c(iters = median(x$value)))

ggplot(p3, aes(dim, iters, group = alg, linetype = alg)) +
  geom_line() +
  facet_grid(nout ~ radius)



### Simulation 2: ####################################################

load(sim.file("sim2"))
str(sim2)


### Median distances and weighted RSS for the robust algorithm:

p4 <- subset(sim2, alg == "robust" & nout == 100 & radius == 15 &
             perf %in% c("dist1", "dist2", "wrss"))
p4 <- ddply(p4, c("dim", "n", "nout", "radius", "k", "alg", "perf"),
            function(p)
            c(value = median(p$value)))
p4$panel <- ifelse(p4$perf != "wrss", "Distance", "Weighted RSS")

ggplot(p4, aes(k, value, group = perf)) +
  geom_line() + geom_point() +
  facet_grid(panel ~ ., scales = "free_y")


