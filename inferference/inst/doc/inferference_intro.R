## ----, echo = TRUE-------------------------------------------------------
# load the package
library(inferference)

## ----, echo=TRUE, eval=FALSE---------------------------------------------
#  alphaz = seq(0, 1, by=.05) # needed to compute truth
#  theta.base <- c(.5, -0.788, -2.953, -0.098, -0.145, 0.351)
#  theta.sim <- c(0.2727, -0.0387, .2719, 1.0859)
#  sample_sim <- sim_interference(n = 3000, N = 250, nsims = 1,
#                                 base.parameters = theta.base,
#                                 parameters = theta.sim,
#                                 alphas = alphaz)
#  interference_sample <- sample_sim$sims[[1]]
#  #save(sample_data, file = "data/sample_data.rda")

## ------------------------------------------------------------------------
head(vaccinesim)

## ----, echo=TRUE---------------------------------------------------------
sample1 <- interference(
  allocations = c(.3, .45,  .6), # a minimum of 2 is required
  formula = y | A | B ~ X1 + X2 + (1|group) | group, 
  data = vaccinesim, # name of the data frame
  randomization = 2/3,
  method = 'simple' # speeds up grad()
  )

## ----, echo = TRUE-------------------------------------------------------
print(sample1)

## ----, echo = TRUE-------------------------------------------------------
# the matrix of weights
head(sample1$weights)

## ----, echo = TRUE-------------------------------------------------------
head(sample1$weightd[ , , 1]) # For the 1st allocation 

## ----, echo = TRUE-------------------------------------------------------
head(sample1$scores) 

## ----, echo = TRUE-------------------------------------------------------
sample1$estimates 

## ----, echo=TRUE---------------------------------------------------------
direct_effect(sample1, .3) #DE(.3)/ unfortunately de() is a function in utils.
ie(sample1, .3, .6) #IE(.3, .6) 
indirect_effect(sample1, .3, .6) # same
te(sample1, .3, .6) #TE(.3, .6)
oe(sample1, .3, .6) #OE(.3, .6)

# You can use the same functions to get all the effects. This will aid in 
# plotting (see below.)

# All direct effect of constrast Y(0) - Y(1)
direct_effect(sample1) 

# to get Y(1) - Y(0) use trt.lvl1 argument
direct_effect(sample1, trt.lvl1 = 1) 

# Other functions work similarly

# all indirect effects where treatment = 0 compare to an allocation of .3
ie(sample1, .3) 

## ----, echo=TRUE---------------------------------------------------------
sample2 <- interference(
  allocations = c(.3, .6), # a minimum of 2 is required
  formula = y | A | B ~ X1 + X2 | group,
  data = vaccinesim, # name of the data frame
  model_method = 'glm',
  model_options = list(family = binomial),
  randomization = 2/3,
  method = 'simple' # speeds up grad()
  )
print(sample2)

## ----, echo=TRUE---------------------------------------------------------
sample3 <- interference(
  allocations = c(.3, .6), # a minimum of 2 is required
  formula = y | A | B ~ X1 + X2 | group,
  data = vaccinesim, # name of the data frame
  model_method = 'oracle',
  model_options = list(fixed.effects = c(0.2727, -0.0387, .2719),
                       random.effects = NULL),
  causal_estimation_options = list(variance_estimation = 'naive'),
  randomization = 2/3,
  method = 'simple' # speeds up grad()
  )
print(sample3)

## ----, echo=TRUE---------------------------------------------------------
sample4 <- interference(
  allocations = c(.3, .6), # a minimum of 2 is required
  formula = y | A ~ X1 + X2 | group,
  data = vaccinesim, # name of the data frame
  model_method = 'glm',
  causal_estimation_options = list(variance_estimation = 'robust'),
  method = 'simple' # speeds up grad()
  )
print(sample4)

## ----, echo = TRUE-------------------------------------------------------
myFUN <- function(b, x, pos){
  return(.5 * dnorm(b))
}

integrate(myFUN, l = -Inf, u = Inf) # returns 0.5

sample5 <- interference(
  propensity_integrand = 'myFUN',
  allocations = c(.3, .6), # a minimum of 2 is required
  formula = y | A ~ X1 + X2 | group,
  data = vaccinesim, # name of the data frame
  model_method = 'glm',
  causal_estimation_options = list(variance_estimation = 'naive'),
  # won't work for 'robust' unless addition arguments are defined
  method = 'simple' # speeds up grad()
  )
print(sample5)

## ----, echo = TRUE-------------------------------------------------------
sample6 <- interference(
  allocations = seq(.3, .6, by = .01), 
  formula = y | A ~ X1 + X2 | group,
  data = vaccinesim, # name of the data frame
  model_method = 'glm',
  causal_estimation_options = list(variance_estimation = 'robust'),
  method = 'simple' # speeds up grad()
  )

# Look at the direct effects
deff <- direct_effect(sample6)

# Plot the point estimates
plot(deff$alpha1, deff$estimate, type = 'l')

