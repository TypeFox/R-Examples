
## ------------------------------------------------------------------------
# analysis for fish data

library(nicheROVER)
data(fish) # 4 fish, 3 isotopes
aggregate(fish[2:4], fish[1], mean) # isotope means calculated for each species


## ----fig.width=7, fig.height=2.5-----------------------------------------
# fish data
data(fish)

# generate parameter draws from the "default" posteriors of each fish
nsamples <- 1e3
system.time({
 fish.par <- tapply(1:nrow(fish), fish$species,
                    function(ii) niw.post(nsamples = nsamples, X = fish[ii,2:4]))
})

# various parameter plots
clrs <- c("black", "red", "blue", "orange") # colors for each species

# mu1 (del15N), mu2 (del13C), and Sigma12
par(mar = c(4, 4, .5, .1)+.1, mfrow = c(1,3))
niche.par.plot(fish.par, col = clrs, plot.index = 1)
niche.par.plot(fish.par, col = clrs, plot.index = 2)
niche.par.plot(fish.par, col = clrs, plot.index = 1:2)
legend("topleft", legend = names(fish.par), fill = clrs)

# all mu (del15N, del13C, del34S)
niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = FALSE)
legend("topleft", legend = names(fish.par), fill = clrs)


## ----fig.width=7, fig.height=10------------------------------------------
# all mu and Sigma
par(mar = c(4.2, 4.2, 2, 1)+.1)
niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = TRUE)
legend("topright", legend = names(fish.par), fill = clrs)


## ----fig.width=7, fig.height=6-------------------------------------------
# 2-d projections of 10 niche regions
clrs <- c("black", "red", "blue", "orange") # colors for each species
nsamples <- 10
fish.par <- tapply(1:nrow(fish), fish$species,
                  function(ii) niw.post(nsamples = nsamples, X = fish[ii,2:4]))

# format data for plotting function
fish.data <- tapply(1:nrow(fish), fish$species, function(ii) X = fish[ii,2:4])

niche.plot(niche.par = fish.par, niche.data = fish.data, pfrac = .05,
          iso.names = expression(delta^{15}*N, delta^{13}*C, delta^{34}*S),
            col = clrs, xlab = expression("Isotope Ratio (per mil)"))


## ------------------------------------------------------------------------
# niche overlap plots for 95% niche region sizes
nsamples <- 1000
fish.par <- tapply(1:nrow(fish), fish$species,
                  function(ii) niw.post(nsamples = nsamples, X = fish[ii,2:4]))

# Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher accuracy.
# the variable over.stat can be supplied directly to the overlap.plot function

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1e3, alpha = c(.95, 0.99))

#The mean overlap metrics calculated across iteratations for both niche 
#region sizes (alpha = .95 and alpha = .99) can be calculated and displayed in an array.
over.mean <- apply(over.stat, c(1:2,4), mean)*100
round(over.mean, 2)
over.cred <- apply(over.stat*100, c(1:2, 4), quantile, prob = c(.025, .975), na.rm = TRUE)
round(over.cred[,,,1]) # display alpha = .95 niche region



## ----fig.width=7, fig.height=6-------------------------------------------
# Overlap plot.Before you run this, make sure that you have chosen your 
#alpha level.
clrs <- c("black", "red", "blue", "orange") # colors for each species
over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1e3, alpha = .95)
overlap.plot(over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
            xlab = "Overlap Probability (%) -- Niche Region Size: 95%")


