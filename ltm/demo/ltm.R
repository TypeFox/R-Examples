
## This demo presents some feature of the ltm() function


## Fit the Rasch model using ltm()
ltm(LSAT ~ z1, constr = cbind(1:length(LSAT), 2, 1))
rasch(LSAT, constr = cbind(length(LSAT) + 1, 1))


## ltm() fits latent trait models up to two latent variables
fit <- ltm(WIRS ~ z1 + z2)

## use the following to produce the plot of standardized loadings,
plot(fit, type = "loadings")

## and use the following to produce the Item Characteristic Surfaces
plot(fit, ticktype = "detailed", theta = 30, phi = 30, expand = 0.5, d = 2, cex = 0.7)


## ltm() can also be used to include nonlinear latent terms
fit <- ltm(Mobility ~ z1 * z2)

plot(fit, ticktype = "detailed", theta = 30, phi = 30, expand = 0.5, d = 2, cex = 0.7)
