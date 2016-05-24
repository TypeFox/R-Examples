#################################################

# these are the source files for nicheROVER.
# this should be the easiest place to do all the editing; everything is in the same place.
# martin lysy, ashley statsko, heidi swanson.  july 2014.

#################################################

#'@title (Niche) (R)egion and Niche (Over)lap Metrics for Multidimensional Ecological Niches.
#'@description This package uses a probabilistic method to calculate niche regions and
#'pairwise niche overlap using multidimensional niche indicator data (e.g., stable
#'isotopes, environmental variables, etc.). The niche region is defined as the joint
#'probability density function of the multidimensional niche indicators at a user-defined
#'probability alpha (e.g., 95\%).  Uncertainty is accounted for in a Bayesian framework,
#'and the method can be extended to three or more indicator dimensions.  It provides
#'directional estimates of niche overlap, accounts for species-specific distributions in
#'multivariate niche space, and produces unique and consistent bivariate projections of
#'the multivariate niche region. See Swanson et al. (2014) for a detailed description and
#'worked example below using fish stable isotope data.
#'@references Heidi K. Swanson, Martin Lysy, Ashley D. Stasko, Michael Power, Jim D. Johnson, and James D. Reist (2014).  ``What Would Hutchinson Think?  A Probabilistic Quantification of Multidimensional Ecological Niches and Niche Overlap''.  \emph{Ecology: Statistical Reports} (accepted).
#'@examples
#'# analysis for fish data
#'
#'data(fish) # 4 fish, 3 isotopes
#'aggregate(fish[2:4], fish[1], mean) # isotope means per fish
#'
#'# random draws from posterior distribution with default prior
#'nsamples <- 500
#'fish.par <- tapply(1:nrow(fish), fish$species,
#'                   function(ii) niw.post(nsamples = nsamples, X = fish[ii,2:4]))
#'
#'# display p(mu | X) and p(Sigma | X) for each fish
#'clrs <- c("black", "red", "blue", "orange") # colors for each species
#'par(mar = c(4.2, 4.2, 2, 1)+.1)
#'niche.par.plot(fish.par, col = clrs)
#'legend(x = "topright", legend = names(fish.par), fill = clrs)
#'
#'# 2-d projections of 10 niche regions
#'nsamples <- 10
#'fish.par <- tapply(1:nrow(fish), fish$species,
#'                   function(ii) niw.post(nsamples = nsamples, X = fish[ii,2:4]))
#'
#'# format data for plotting function
#'fish.data <- tapply(1:nrow(fish), fish$species, function(ii) X = fish[ii,2:4])
#'
#'niche.plot(niche.par = fish.par, niche.data = fish.data, pfrac = .05,
#'           iso.names = expression(delta^{15}*N, delta^{13}*C, delta^{34}*S),
#'           col = clrs, xlab = expression("Isotope Ratio (\u2030)"))
#'
#'# niche overlap plots for 95% niche region sizes
#'
#'# overlap calculation.  use nsamples = nprob = 1e4 for higher accuracy.
#'nsamples <- 500
#'over.stat <- overlap(fish.par, nreps = nsamples, nprob = nsamples, alpha = .95)
#'
#'# overlap plot
#'overlap.plot(over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
#'             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")
#'@docType package
#'@importFrom mvtnorm rmvnorm
#'@name nicheROVER
NULL
