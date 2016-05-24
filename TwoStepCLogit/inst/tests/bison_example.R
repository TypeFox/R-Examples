# Some descriptive statistics about the data set:
ddim(formula = Y ~ strata(Strata) + cluster(Cluster), data = bison)
 
# Model 1: covariates meadow, biomass and biomass^2
# Random effects in front of biomass and biomass^2
# Main diagonal covariance structure for D
Fit1 <- Ts.estim(formula = Y ~ meadow + biomass + I(biomass^2) + 
        strata(Strata) + cluster(Cluster), data = bison, 
        random = ~ biomass + I(biomass^2), all.m.1=FALSE, D="UN(1)")

Fit1

# Model 2: only covariates biomass and biomass^2
# Random effects in front of biomass and biomass^2
# Main diagonal covariance structure for D
Fit2 <- Ts.estim(formula = Y ~ biomass + I(biomass^2) + strata(Strata) + 
        cluster(Cluster), data = bison, all.m.1=FALSE, D="UN(1)")
Fit2
# Results reported in Table 2 of Craiu et al. (2011).

