## comparison with the glasso package
library(simone)
if (!("glasso" %in% c(.packages(all.available=TRUE)))) {
  stop("\nYou need the 'glasso' package to run this demo")
}
require(glasso)

## Data set and graph generation
p    <- 80
g    <- rNetwork(p, pi=p)
data <- rTranscriptData(n=2*p,g)
attach(data)

cat("---------------------------------------------------------------",
    ""                                                               ,
    "SIMoNe embeds the GLasso algorithm of Friedman et al: when n>p ",
    "the very same algorithm is used to solve  the underlying  Lasso",
    "problems (the  'shooting' or  'pathwise coordinate' procedure).",
    ""                                                               , 
    "When n <= p,  the underlying  Lasso problems are solved with an",
    "active-set algorithm that produces sparser solutions."          ,
    "---------------------------------------------------------------",
    sep="\n")

cat("\n\n Check consistency between the glasso package and the simone package")

rho <-  0.2 # Fix the penalty level to 0.2
theta.glasso <- glasso(var(scale(X)), rho=rho)$wi # GLasso package

control      <- setOptions(penalties=rho, edges.steady="graphical.lasso", )
theta.simone <- simone(X, control=control)$networks[[1]] # SIMoNe

cat("\n norm2-error between the two solutions: ",
    sqrt(sum((theta.simone - theta.glasso)^2)))

cat("\n number of sign inconsistencies: ",
    sum(sign(theta.simone) != sign(theta.glasso)))

net.simone <- structure(list(A        = sign(theta.simone),
                             Theta    = theta.simone,
                             directed = FALSE,
                             clusters = factor(rep("N",p)),
                             name     = "simone package"), class="simone.network")

net.glasso <- structure(list(A        = sign(theta.glasso),
                             Theta    = theta.glasso,
                             directed = FALSE,
                             clusters = factor(rep("N",p)),
                             name     = "glasso package"), class="simone.network")
plot(net.glasso,net.simone)

cat("\n\n")
readline("Press enter for the full path (thanks to simone :) )")
control <- setOptions(edges.steady="graphical.lasso")
out     <- simone(X, control=control) # SIMoNe
plot(out, ask=FALSE)

detach(data)
