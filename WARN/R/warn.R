# =====R code=========================
# Weaning age reconstruction using d15N.
# ==============================
# 2013-01-14: Tsutaya T: All files are combined and variable names are integrated.
# 2013-01-22: Tsutaya T: Added warpper function for optim.
# 2013-02-09: Tsutaya T: Prior was changed.
# 2013-04-11: Tsutaya T: Added conditionnig on tolerances.
# 2014-11-03: Tsutaya T: Deleted call to 'MASS' and added "MASS::" to width.SJ.
# ==============================
# OBJECTIVE ==========
# This program performs Apporoximate Bayesian Computation with SMC
#  adopting Partial Rejection Control (PRC)
#  in order to calculates posterior distributuions
#  of the weaning ages and weaning parameters.
# ==============================
# FUNCTION DEFINITIONS ==========
# warn --------------------
## Process assigned data.
SubtractAgeResidue <- function(age){
# Calculate age before one unit time and its residue for a given age vector.
#
# Args:
#  age: A vector of ages.
#
# Returns:
#  Dataframe of original ages, their residues and their subtracts.
#
# Exception handling:
#  When age is 0 years, its subtract is converted to 0.
	residue <- age %% 1.0
	residue <- ifelse(residue == 0, 1.0, residue)
	subtract <- age - residue
  subtract <- ifelse(subtract == -1.0, 0, subtract)
  results <- matrix(c(age, residue, subtract), length(age), 3)
  colnames(results) <- c("age", "residue", "subtract")
	# Returns the results.
	return(results)
}

## Set the prior distributions.
SamplePar <- function(mu.t1, sigma.t1, mu.t2, sigma.t2,
  mu.enrich, sigma.enrich, mu.wnfood, sigma.wnfood,
  mu.epsilon, sigma.epsilon){
# Sample parameters from the prior normal distribution
#  defined with thier hyper parameters.
#
# args:
#  Hyper parameters (mu: mean, sigma: 1SD):
#   t1: Age at the beginning of weaning.
#   t2: Age at the completion of weaning.
#   enrich: Enrichment factor between mother and infant.
#   wnfood: d15N value of weaning food.
#   epsilon: Hyper parameter (1SD) for the individual epsilon.
#
# returns:
#  Dataframe of all parameters.
  while(TRUE){
    par.t1 <- rnorm(1, mu.t1, sigma.t1)
    par.t2 <- rnorm(1, mu.t2, sigma.t2)
    ifelse(par.t1 < par.t2 & par.t1 >= 0, break, FALSE)
  }
  par.enrich <- rnorm(1, mu.enrich, sigma.enrich)
  par.wnfood <- rnorm(1, mu.wnfood, sigma.wnfood)
  par.epsilon <-  rnorm(1, mean = mu.epsilon, sd = sigma.epsilon)

  return(c(t1 = par.t1, t2 = par.t2, enrich = par.enrich,
    wnfood = par.wnfood, epsilon = par.epsilon))
}

GenerateEpsilon <- function(sigma, num.individual){
# Sample individual epsilons under the given sigma (distribution)
#  which was defined by the hyper parameter of epsilon.
#
# args:
#  sigma: Hyper parameter of SD for distribution of epsilon.
#   Individual epsilons are given as N~(0, sigma.epsilon^2).
#  num.individual: The number of individuals in observed data.
#
# returns:
#  A vector of individual epsilons.
  epsilons <- rnorm(n = num.individual, mean = 0, sd = abs(sigma))

  return(epsilons)
}

## Calculate collagen turnover rates.
SetColTurnover <- function(x){
# Set turnover rate of bone collagen at given age.
#  This function is only applicable for nonadults under 20 years of age.
#  Turnover rate is represented as that of 1 year from t - 1.0 to t year.
#
# args:
#  x: An age at which we want to calculate the rate
  return(1.778e+00 - 4.121e-01 * x + 5.029e-02 * x^2 - 
    2.756e-03 * x^3 + 5.325e-05 * x^4)
}

ArrangeColTurnover <- function(age, subtract){
# Calculate and arrange "collagen" turnover rate (modeling + remodeling).
#  If turnover rate is larger than 1, it must be 1.
#  If an objective individual is 0 years old, turnover rate of it must be 1.
#  Turnover rates during period of [subtract -> age] will be calculated.
#
# args:
#  age: A vector of reference or observed ages
#   with which we want to calculate the turnover rates.
#  subtract: A vector of their subtracts (unit time is <= 1.0 year).
#
# returns:
#  A vector of combined bone "collagen" turnover rate 
#   (modeling: "growth" + remodeling: "turnover") 
#   for the given age dataframe used in further computation.
#  Calculated reference turnover[1] corresponds to ages between 0 to 0 years,
#   reference turnover[2] corresponds to ages between 0 to 1.0 years, and so on.
  age.ceil <- ceiling(age)
  int.all <-  mapply(function(a, b){ # from a to b
    return(integrate(SetColTurnover, a, b)$value)
  },
  a = subtract, b = age.ceil)
  int.age <-  mapply(function(a, b){ # from a to b
    return(integrate(SetColTurnover, a, b)$value)
  },
  a = subtract, b = age)

  turnover <- SetColTurnover(age.ceil) * int.age / int.all

  turnover <- ifelse(turnover >= 1, 1, turnover)
  turnover <- ifelse(age == 0, 1, turnover)

  return(turnover)
}

## Weaning model.
CalcNonmilkProp <- function(age, t1, t2, form = "parabolic"){
# Calculate proportion of wnfood food for a given age and weaning ages.
#
# args:
#  age: An age at which we want to calculate the proportion.
#  t1 and t2: Weaning ages.
#  form: The form of change in d15N values during weaning.
#
# returns:
#  The proportion of wnfood food.
	if(form == "linear"){
		proportion <- (age - t1) / (t2 - t1)  # linear form
	}else if(form == "parabolic"){
		proportion <- ((age - t1) / (t2 - t1))^2  # parabolic form
	}else if(form == "reverse"){
		proportion <- 1 - ((t2 - age) / (t2 - t1))^2  # reverse parabolic form
	}else if(form == "sigmoid"){
		front <- 1/2 + (1/2) * ((2 * age - t1 - t2) / (t2 - t1))^(1/3)
		front <- ifelse(is.nan(front), 0, front)
		post <- 1/2 - (1/2) * ((2 * age - t1 - t2) / (t1 - t2))^(1/3)
		post <- ifelse(is.nan(post), 0, post)
		post <- ifelse(post== (1/2), 0, post)	
		proportion <- (front + post) # sigmoid form
	}
	return(proportion)
}

IntNonmilkProp <- function(age, residue, subtract, t1, t2, form = "parabolic"){
# Integrate proportions of wnfood food for given ages and weaning ages.
#
# args:
#  age: A vector of reference or observed ages
#   with which we want to calculate the proportions.
#  residue: A vector of their residues (unit time is <= 1.0 years).
#  subtract: A vector of their subtracts (unit time is <= 1.0 years).
#  t1 and t2: Weaning ages.
#  form: The form of change in d15N values during weaning.
#
# returns:
#  The integrated proportions of wnfood food.
#   The return[1] corresponds to the period between 0 and 0 years,
#   the return[2] corresponds to the period between 0 and 1.0 years and so on.
#
# functions:
#  CalcNonmilkProp.
  Int <- function(age, residue, subtract, t1, t2, form){
  	if(age == 0){
	  	return(0 / residue)
	  }else if(age < t1){
		  return(0 / residue)
	  }else if(age >= t1 && age < t2 && subtract < t1){
		  inted.value <- integrate(CalcNonmilkProp, t1, age,
        t1 = t1, t2 = t2, form = form)
  		return(inted.value$value / residue)
	  }else if(subtract < t1 && age >= t2){
		  inted.value <- integrate(CalcNonmilkProp, t1, t2,
        t1 = t1, t2 = t2, form = form)
  		return(((age - t2) + inted.value$value) / residue)
	  }else if(subtract >= t1 && age < t2){
		  inted.value <- integrate(CalcNonmilkProp, subtract,
        age, t1 = t1, t2 = t2, form = form)
  		return(inted.value$value / residue)
	  }else if(subtract >= t1 && subtract < t2 && age >= t2){
		  inted.value <- integrate(CalcNonmilkProp, subtract, t2,
        t1 = t1, t2 = t2, form = form)
  		return(((age - t2) + inted.value$value) / residue)
	  }else if(subtract >= t2){
		  return(residue / residue)
	  }
  }
  mapply(Int, age = age, residue = residue, subtract = subtract,
    t1 = t1, t2 = t2, form = form)
}

CalcRefBone <- function(age, residue, subtract, inted.wnfood.prop, enrich,
  female.mean, n.wnfood, turnover){
# Calculate reference bone d15N values.
#
# args:
#  age: A vector of reference ages with which we want to calculate the proportions.
#   This must start with 0 and be sequential numbers by unit time (e.g. 1.0 years).
#  residue: A vector of their residues (unit time is <= 1.0 years).
#  subtract: A vector of their subtracts (unit time is <= 1.0 years).
#  inted.wnfood.prop: A vector of integrated wnfood proportions.
#   for the reference ages.
#  enrich: Enrichment factor.
#  female.mean: Mean d15N value of adult females.
#  n.wnfood: Hypothesized d15N value of weaning foods.
#  turnover: A vector of integrated turnover ratis for the reference ages.
#
# returns:
#  Vector of d15N values for the reference bone.
  n.diets <- (1 - inted.wnfood.prop) * (enrich + female.mean) +
    inted.wnfood.prop * n.wnfood
  # The n.diets at 0 years old equal to female.mean + enrich.

	# Sequentially calculate the modeled bone d15N values for reference ages.
	n.bones <- numeric(length(age))
	n.bones[1] <- female.mean
	for(i in 2:length(age)){
		n.bones[i] <- n.bones[i - 1] * (1 - turnover[i]) + n.diets[i] * turnover[i]
  }

	return(n.bones)
}

SimulateBone <- function(age, residue, subtract, inted.wnfood.prop,
  enrich, female.mean, n.wnfood, n.bones.ref, epsilons, turnover){
# Simulate bone d15N values for the observed ages.
#
# args:
#  age: A vector of observed ages with which we want to calculate the proportions.
#  residue: A vector of their residues (unit time is <= 1.0 years).
#  subtract: A vector of their subtracts (unit time is <= 1.0 years).
#  inted.wnfood.prop: A vector of integrated wnfood proportions.
#   for the observed ages.
#  inted.turnover: A vector of integrated turnover ratis for the observed ages.
#  enrich: Enrichment factor.
#  female.mean: Mean d15N value of adult females.
#  n.wnfood: Hypothesized d15N value of weaning foods.
#  n.bones.ref: Reference bone d15N values for its unit time.
#  epsilon: A vector of individual error terms.
#  turnover: A vector of integrated turnover ratis for the observed ages.
#
# returns:
#  Vector of d15N values for the simulated bone.
  n.diets <- (1 - inted.wnfood.prop) * (enrich + female.mean) +
    inted.wnfood.prop * n.wnfood

	# Calculate the modeled bone d15N values for observed ages.
  unit.time <- 1.0
  n.bones.prev <- n.bones.ref[(subtract / unit.time) + 1]
	n.bones <- n.bones.prev * (1 - turnover) + n.diets * turnover
  # If age is 0 year, its bone must be equal to the d15N value of mother.
  n.bones[which(age == 0)] <- female.mean

	return(n.bones + epsilons)
}

## Optimization.
ConditionAge <- function(t1, t2){
# Condition weaning ages as 0 <= t1 < t2.
#
# args:
#  t1: The canditate age at the beginning of weaning.
#  t2: The canditate age at the completion of weaning.
#
# returns:
#  A vector of weaning ages.
	# Conditions: t1 >= 0 and t2 >= 0.
	if(t1 <= 0){	t1 <- 0}
	if(t2 <= 0){	t2 <- 0}

	# Condition: t1 < t2.
	if(t1 > t2){ A <- t1; t1 <- t2; t2 <- A}

	# Condition: t1 != t2.
	if(t1 == t2){	t2 <- t2 + 0.000001}
	
	return(c(t1 = t1, t2 = t2))
}

CalcDistanceForOptim <- function(age.ref, residue.ref, subtract.ref,
  age.obs, residue.obs, subtract.obs, female.mean, epsilon = 0, 
  turnover.ref, turnover.obs, d15N, form = "parabolic"){
# Calculate the distance between observed and simulated data for optimization.
#  Distance is represented as a mean of least square distances.
#
# args:
#  age: A vector of observed ages with which we want to calculate the distance.
#  residue: A vector of their residues (unit time is <= 1.0 years).
#  subtract: A vector of their subtracts (unit time is <= 1.0 years).
#  female.mean: Mean d15N value of adult females.
#  epsilon: SD for the distribution of individual epsilons.
#  turnover.ref: A vector of tuenover rate for reference ages.
#  turnover.obs: A vector of tuenover rate for observed ages.
#  d15N: A vector of observed d15N values for the nonadults.
#  form: The form of change in d15N values during weaning.
#
# returns:
#  The distance.
#
# functions:
#  ConditionAge
#  CalcDistance
  return(function(par, ...){ # , age, residue, subtract,
#    female.mean, epsilon = 0, turnover.ref, turnover.obs, d15N){
    t1 = par[1]
    t2 = par[2]
    enrich = par[3]
    wnfood = par[4]

    weaning.age <- ConditionAge(t1 = t2, t2 = t2)
    t1 <- weaning.age[1]
    t2 <- weaning.age[2]

    return(CalcDistance(age.ref = age.ref, residue.ref = residue.ref,
      subtract.ref = subtract.ref, age.obs = age.obs,
      residue.obs = residue.obs, subtract.obs = subtract.obs,
      t1 = t1, t2 = t2, enrich = enrich, wnfood = wnfood,
      female.mean = female.mean, epsilon = 0, turnover.ref = turnover.ref,
      turnover.obs = turnover.obs, d15N = d15N, form = form))
  })
}

OptimizePar <- function(par.initial, age.ref, residue.ref, subtract.ref,
  age.obs, residue.obs, subtract.obs, female.mean, epsilon = 0,
  turnover.ref, turnover.obs, d15N, form = "parabolic"){
# Perform optimization for a given dataset.
#
# args:
#  par.initial: Initial values for the parameters (t1, t2, enrich, and wnfood).
#  age: A vector of observed ages with which we want to calculate the distance.
#  residue: A vector of their residues (unit time is <= 1.0 years).
#  subtract: A vector of their subtracts (unit time is <= 1.0 years).
#  female.mean: Mean d15N value of adult females.
#  epsilon: SD for the distribution of individual epsilons.
#  turnover.ref: A vector of tuenover rate for reference ages.
#  turnover.obs: A vector of tuenover rate for observed ages.
#  d15N: A vector of observed d15N values for the nonadults.
#  form: The form of change in d15N values during weaning.
#
# returns:
#  List of optimized parameters and other results.
#
# functions:
#  CalcDistanceForOptime
#  ConditionAge
	# optimize (minimize) the parameters
	par.optimized <- optim(par = par.initial, 
		fn = CalcDistanceForOptim(
      age.ref = age.ref,
      residue.ref = residue.ref,
      subtract.ref = subtract.ref,
      age.obs = age.obs,
      residue.obs = residue.obs,
      subtract.obs = subtract.obs,
      female.mean = female.mean,
      epsilon = 0,
      turnover.ref = turnover.ref,
      turnover.obs = turnover.obs,
      d15N = d15N,
      form = form
    ),
		control = list(maxit = 10000, ndeps = 1e-2, reltol = 1e-7)
	)

	# conditioning the weaning ages
	weaning.age <- ConditionAge(
    t1 = par.optimized$par[1], t2 = par.optimized$par[2])
	par.optimized$par[1] <- weaning.age[1]
	par.optimized$par[2] <- weaning.age[2]

	return(par.optimized)
}

WrapperOptim <- function(age, d15N, female.mean,
  par.initial = c(0.5, 3, 1.9, female.mean), 
  form = "parabolic", ...){
# Wrapper function for optimization of the weaning parameters.
#
# args:
#  age: A vector of observed ages for the subadults.
#  d15N: A vector of observed d15N values for the subadults.
#  female.mean: Female mean and SD d15N values.
#  par: Initial values for parameters to be optimized over.
#  form: The form of change in d15N values during weaning.
#  ...: Additional arguments passed to optim().
#
# returns:
#  Returned values of the function optim().
#
# functions:
#  CalcDistanceForOptim
#  SubtractAgeResidue
#  ArrangeTurnoverRate
  # Arrange ages.
  age.ref = seq(0, 10, 1.0)

  age.reference <- SubtractAgeResidue(age.ref)
  age.observed <- SubtractAgeResidue(age)
 
  age.ref <- age.reference[ , 1]
  residue.ref <- age.reference[ , 2]
  subtract.ref <- age.reference[ , 3]

  age.obs <- age.observed[ , 1]
  residue.obs <- age.observed[ , 2]
  subtract.obs <- age.observed[ , 3]

  # Arrange turnover rate.
  turnover.ref <- ArrangeColTurnover(age.ref, subtract.ref)
  turnover.obs <- ArrangeColTurnover(age.obs, subtract.obs)

	# Optimize (minimize) the parameters.
	par.optimized <- optim(par = par.initial, 
		fn = CalcDistanceForOptim(
      age.ref = age.ref,
      residue.ref = residue.ref,
      subtract.ref = subtract.ref,
      age.obs = age.obs,
      residue.obs = residue.obs,
      subtract.obs = subtract.obs,
      female.mean = female.mean,
      epsilon = 0,
      turnover.ref = turnover.ref,
      turnover.obs = turnover.obs,
      d15N = d15N,
      form = form
    ), ...
    # control = list(maxit = 10000, ndeps = 1e-2, reltol = 1e-7)
	)

	return(par.optimized)
}

## Perform ABC-PRC.
CalcDistance <- function(age.ref, residue.ref, subtract.ref,
  age.obs, residue.obs, subtract.obs, t1, t2, enrich, wnfood,
  female.mean, epsilon, turnover.ref, turnover.obs, d15N,
  dist.method = "square", zero.omit = FALSE, form = "parabolic"){
# Calculate the distance between observed and simulated data.
#  Distance is represented as a mean of least square distances.
#
# args:
#  age.ref: A vector of reference ages
#   with which we want to calculate the distance.
#  residue.ref: A vector of their residues (unit time is <= 1.0 years).
#  subtract.ref: A vector of their subtracts (unit time is <= 1.0 years).
#  age.obs: A vector of observed ages
#   with which we want to calculate the distance.
#  residue.obs: A vector of their residues (unit time is <= 1.0 years).
#  subtract.obs: A vector of their subtracts (unit time is <= 1.0 years).
#  t1 and t2: Weaning ages.
#  enrich: Enrichment factor.
#  wnfood: d15N value of wnfood food.
#  female.mean: Mean d15N value of adult females.
#  epsilon: SD for the distribution of individual epsilons.
#  turnover.ref: A vector of tuenover rate for reference ages.
#  turnover.obs: A vector of tuenover rate for observed ages.
#  d15N: A vector of observed d15N values for the nonadults.
#  dist.method: Which method ("square" or "as.is")
#   should be applied to calculate the distance.
#   The former calculates mean square distance and latter Euclidean.
#  zero.omit: Omit or not individuals in 0 years of age.
#  form: The form of change in d15N values during weaning.
#
# returns:
#  When distance = "square": Mean of square distances.
#  When distance = "as.is": 1SD of Euclidean distances.
#
# functions:
#  IntNonMilkProp
#  CalcRefBone
#  GenerateEpsilon
#  SimulateBone
  inted.wnfood.ref <- IntNonmilkProp(age = age.ref, residue = residue.ref,
    subtract = subtract.ref, t1 = t1, t2 = t2, form = form)
  inted.wnfood.obs <- IntNonmilkProp(age = age.obs, residue = residue.obs,
    subtract = subtract.obs, t1 = t1, t2 = t2, form = form)

  n.bone.ref <- CalcRefBone(age.ref, residue.ref, subtract.ref,
    inted.wnfood.ref, enrich, female.mean, wnfood, turnover.ref)

  epsilons <- GenerateEpsilon(epsilon, length(age.obs))

  n.bone.obs <- SimulateBone(age.obs,
    residue.obs,
    subtract.obs,
    inted.wnfood.obs,
    enrich,
    female.mean,
    wnfood,
    n.bone.ref,
    epsilons,
    turnover.obs)

  if(zero.omit == TRUE){
    is.zero <- age.obs == 0
  }else{
    is.zero <- FALSE
  }

  if(dist.method == "square"){
    distance <- mean(((d15N - n.bone.obs)^2)[!is.zero])
  }else if(dist.method == "as.is"){
    distance <- sd((d15N - n.bone.obs)[!is.zero])
  }else if(dist.method == "mean"){
    distance <- mean((d15N - n.bone.obs)[!is.zero])
  }

  return(distance)
}

CalcDistanceMDE <- function(age.ref = seq(0, 10, 1.0), par.mde,
  age.obs, d15N, female.mean, form = "parabolic"){
# Calcurate variation (scatterness) under the MLE framework
#  for a given observed dataset and MDE parameters.
#
# args:
#  age.ref: A vector of reference ages (unit time).
#  par.mde: MDEs of each parameters, c(t1, t2, enrich, wnfood).
#  pars: Hyper parameters in the order of c(mu.t1, sigma.t1,
#  age.obs: A vector of observed ages
#   with which we want to calculate the posterior.
#  d15N: A vector of observed d15N values for the nonadults.
#  female.mean: Mean d15N value of adult females.
#  form: The form of change in d15N values during weaning.
#
# returns:
#  A mean of least square distances.
#   
# functions:
#  SubtractAgeResidue
#  ArrangeColTurnover
#  OptimizePar
#  CalcDistance
  age.reference <- SubtractAgeResidue(age.ref)
  age.observed <- SubtractAgeResidue(age.obs)
 
  age.ref <- age.reference[ , 1]
  residue.ref <- age.reference[ , 2]
  subtract.ref <- age.reference[ , 3]

  age.obs <- age.observed[ , 1]
  residue.obs <- age.observed[ , 2]
  subtract.obs <- age.observed[ , 3]

  turnover.ref <- ArrangeColTurnover(age.ref, subtract.ref)
  turnover.obs <- ArrangeColTurnover(age.obs, subtract.obs)

  square.opt <- CalcDistance(age.ref = age.ref,
    residue.ref = residue.ref,
    subtract.ref = subtract.ref,
    age.obs = age.obs,
    residue.obs = residue.obs,
    subtract.obs = subtract.obs,
    t1 = par.mde[1],
    t2 = par.mde[2],
    enrich = par.mde[3],
    wnfood = par.mde[4],
    female.mean = female.mean,
    epsilon = 0,
    turnover.ref = turnover.ref,
    turnover.obs = turnover.obs,
    d15N = d15N,
    dist.method = "square",
    zero.omit = FALSE,
    form = form)

  return(square.opt)
}

SetTransition <- function(mu){
# Transition kernel used in Sequential Monte Carlo.
#
# args:
#  mu: A previous parameter from which we want to calculate perturbed parameter.
#
# returns:
#  Next canditate parameter value.
  rnorm(1, mu, 0.10)
}

CalcTransition <- function(from, to){
# Clculate probability densitity of the transition in SMC.
#
# args:
#  from: Previous parameter value.
#  to: Perturbed next parameter value.
#
# returns:
#  Density of the transition.
  dnorm(to - from, 0, 0.10)
}

PerformABCPRC <- function(threshold, tolerances, num.particle, mu.t1, sigma.t1,
  mu.t2, sigma.t2, mu.enrich, sigma.enrich, mu.wnfood, sigma.wnfood,
  mu.epsilon, sigma.epsilon, age.ref, residue.ref, subtract.ref,
  age.obs, residue.obs, subtract.obs, female.mean,
  turnover.ref, turnover.obs, d15N, form = "parabolic"){
# Perform ABC-SMC with Partial Rejection Control (PRC).
#
# args:
#  threshold: Mean least square distances between observed and optimized model.
#  tolerances: A vector of sequentially-decreasing tolerances.
#  num.particle: The number of particles.
#  Hyper parameters (mu: mean, sigma: 1SD):
#   t1: Age at the beginning of weaning.
#   t2: Age at the completion of weaning.
#   enrich: Enrichment factor between mother and infant.
#   wnfood: d15N value of weaning food.
#   epsilon: Hyper parameter (1SD) for the individual epsilon.
#  age.ref: A vector of reference ages
#   with which we want to calculate the distance.
#  residue.ref: A vector of their residues (unit time is <= 1.0 years).
#  subtract.ref: A vector of their subtracts (unit time is <= 1.0 years).
#  age.obs: A vector of observed ages
#   with which we want to calculate the distance.
#  residue.obs: A vector of their residues (unit time is <= 1.0 years).
#  subtract.obs: A vector of their subtracts (unit time is <= 1.0 years).
#  female.mean: Mean d15N value of adult females.
#  turnover.ref: A vector of tuenover rate for reference ages.
#  turnover.obs: A vector of tuenover rate for observed ages.
#  d15N: A vector of observed d15N values for the nonadults.
#  form: The form of change in d15N values during weaning.
#
# returns:
#  A list of target theta particle and its weight.
#
# functions:
  # Indicator.
  cat("Remains: ")
  # Set a number of parallelization and hyper parameters.
  parallel <- 1000
  mu.all <- c(mu.t1, mu.t2, mu.enrich, mu.wnfood, mu.epsilon)
  sigma.all <- c(sigma.t1, sigma.t2, sigma.enrich, sigma.wnfood, sigma.epsilon)
  # Start SMC sampling with PRC.
  t <- 1
  while(t <= length(tolerances)){
    # Indicator.
    cat(length(tolerances) - t + 1, "")
    # Set initial values.
    theta.chain <- array(0, dim = c(num.particle + parallel, 5))
    weight.chain <- array(0, dim = c(num.particle, 5))
    i <- 1
    # Computation for a first population.
    if(t ==  1){
      while(i <= num.particle){
        while(TRUE){
          # Set initial value.
          theta.aa <- array(0, dim = c(parallel, 5))
          # Sample from prior.
          for(j in 1:parallel){
            theta.aa[j, ] <- SamplePar(
              mu.t1, sigma.t1, mu.t2, sigma.t2, mu.enrich, sigma.enrich,
              mu.wnfood, sigma.wnfood, mu.epsilon, sigma.epsilon)
          }
          # Calculate S(theta^**) and check the rho.
          distance.aa <- mapply(function(t1, t2, enrich, wnfood, epsilon){
            CalcDistance(
            age.ref = age.ref,
            residue.ref = residue.ref, subtract.ref = subtract.ref,
            age.obs = age.obs, residue.obs = residue.obs,
            subtract.obs = subtract.obs, t1 = t1, t2 = t2,
            enrich = enrich, wnfood = wnfood,
            female.mean =  female.mean, epsilon = epsilon,
            turnover.ref = turnover.ref, turnover.obs = turnover.obs,
            d15N = d15N, form = form)},
            t1 = theta.aa[ , 1], t2 = theta.aa[ , 2],
            enrich = theta.aa[ , 3], wnfood = theta.aa[ , 4],
            epsilon = theta.aa[ , 5]
          )
          theta.tf <- ifelse(distance.aa - threshold <= tolerances[t], TRUE, FALSE)
          # Retry if there is no accepted theta.
          ifelse(sum(theta.tf) <= 0, TRUE, break)
        }
        # Assign accepted theta.
        theta.chain[i:(i + sum(theta.tf) - 1), ] <- theta.aa[theta.tf, ]
        # Increment.
        i <- i + sum(theta.tf)
      }
      # Set weight.
      weight.chain[ , ] <- 1
    }else{
    # Computation for the remaining population(s).
      while(i <= num.particle){
        while(TRUE){
          # Set initial value.
          theta.a <- array(0, dim = c(parallel, 5))
          # Sample canditate theta from previous particle with its weights.
          for(j in 1:5){
            theta.a[ , j] <- try(sample(x = theta.prev[ , j],
              size = parallel, replace = TRUE, prob = weight.prev[ , j]))
          }
          # Propagate theta with transition kernel.
          theta.aa <- apply(theta.a, c(1, 2), SetTransition)
          # Condition for t1 and t2.
          while(TRUE){
            t.tf <- ifelse(theta.aa[ , 1] < theta.aa[ , 2] & theta.aa[ , 1] > 0,
              FALSE, TRUE)
            # Check.
            ifelse(sum(t.tf) == 0, break, TRUE)
            for(j in 1:2){
              theta.a[t.tf, j] <- try(sample(x = theta.prev[ , j],
                size = sum(t.tf), replace = TRUE, prob = weight.prev[ , j]))
            }
            # Propagate theta with transition kernel.
            theta.aa[t.tf, 1:2] <- mapply(SetTransition, mu = theta.a[t.tf, 1:2])
          }
          # Calculate S(theta^**) and check the rho.
          distance.aa <- mapply(function(t1, t2, enrich, wnfood, epsilon){
            CalcDistance(
            age.ref = age.ref,
            residue.ref = residue.ref, subtract.ref = subtract.ref,
            age.obs = age.obs, residue.obs = residue.obs,
            subtract.obs = subtract.obs, t1 = t1, t2 = t2,
            enrich = enrich, wnfood = wnfood,
            female.mean =  female.mean, epsilon = epsilon,
            turnover.ref = turnover.ref, turnover.obs = turnover.obs,
            d15N = d15N, form = form)},
            t1 = theta.aa[ , 1], t2 = theta.aa[ , 2],
            enrich = theta.aa[ , 3], wnfood = theta.aa[ , 4],
            epsilon = theta.aa[ , 5]
          )
          theta.tf <- ifelse(distance.aa - threshold <= tolerances[t],
            TRUE, FALSE)
          # Retry if there is no accepted theta.
          ifelse(sum(theta.tf) <= 0, TRUE, break)
        }
        # Assign accepted theta.
        theta.chain[i:(i + sum(theta.tf) - 1), ] <- theta.aa[theta.tf, ]
        # Increment.
        i <- i + sum(theta.tf)
      }
      # If t = T, weighting does not be performed.
      if(t < length(tolerances)){
        # After required particles are filled: Calculate Density of theta.
        theta.probability <- t(apply(theta.chain[1:num.particle, ], 1,
            function(x){dnorm(x, mean = mu.all, sd = sigma.all)}))
        # Calculate the weight of theta.
        weight.denom <- array(0, dim = c(num.particle, 5))
        for(j in 1:5){
          for(k in 1:num.particle){
            weight.denom[k, j] <- sum(weight.prev[ , j] * 
              CalcTransition(theta.prev[ , j], theta.chain[k, j]))
          }
        }
        weight.chain <- theta.probability / weight.denom
      }
    }
    # Trim theta chain and assign theta.
    theta.prev <- theta.chain[1:num.particle, ]
    # If t = T, weighting does not be performed.
    ifelse(t == length(tolerances), break, TRUE)
    # Normarile the weights.
    weight.sum <- apply(weight.chain, 2, sum)
    weight.tf <- weight.sum == 0
    weight.chain[ , weight.tf] <- 1
    weight.sum[weight.tf] <- num.particle
    weight.normalized <- t(apply(weight.chain, 1, function(y){ y / weight.sum}))
    # Perform PRC in case of skewed distribution in theta (ESS < N/2).
    weight.prev <- array(0, dim = c(num.particle, 5))
    ess <- (1 / apply(weight.normalized^2, 2, sum))
    for(l in 1:5){
      if(ess[l] < (num.particle / 2)){
        theta.prev[ , l] <- try(sample(x = theta.prev[ , l],
          size = num.particle, prob = weight.normalized[ , l]))
        weight.prev[ , l] <- (1 / num.particle)
      }else{
        weight.prev[ , l] <- weight.normalized[ , l]
      }  
    }
    # Increment.
    t <- t + 1
  }
  # Indicator.
  cat("\n")
  return(list(theta = theta.prev, weight = weight.prev))
}

CalcPostDensity <- function(posterior.x, posterior.y = NA,
  is.2d = FALSE, decimal = 1){
# Apply 1d- or 2d-KDE to given posterior particle(s).
#
# args:
#  posterior.x: Posterior particle 1.
#  posterior.y: Posterior particle 2.
#  is.2d: If TRUE, 2d-KDE will be applyed.
#  decimal: Decimal point to be rounded. 1 is recommennded.
#
# dependence:
#  library(MASS): If 2d-KDE is applyed.
#
# returns:
#  A list returned from density() or kde2d().
  times <- 10^(decimal)

  if(is.2d){
    # Preliminary density.
    d1 <- kde2d(x = posterior.x, y = posterior.y,
      h = c(MASS::width.SJ(posterior.x, method = "dpi"),
        MASS::width.SJ(posterior.y, method = "dpi")))

    # Adjust density to given decimal points.
    min.dx <- floor(min(d1$x) * times) / times
    max.dx <- ceiling(max(d1$x) * times) / times
    n.dx <- (max.dx - min.dx) * times + 1

    min.dy <- floor(min(d1$y) * times) / times
    max.dy <- ceiling(max(d1$y) * times) / times
    n.dy <- (max.dy - min.dy) * times + 1

    # Adjested density.
    d2 <- kde2d(x = posterior.x, y = posterior.y,
      h = c(MASS::width.SJ(posterior.x, method = "dpi"),
        MASS::width.SJ(posterior.y, method = "dpi")),
      n = c(n.dx, n.dy), lims = c(min.dx, max.dx, min.dy, max.dy))
  }else{
    # Preliminary density.
    d1 <- density(posterior.x,
      bw = MASS::width.SJ(posterior.x, method = "dpi"))

    # Adjust density to given decimal points.
    min.dx <- floor(min(d1$x) * times) / times
    max.dx <- ceiling(max(d1$x) * times) / times
    n.dx <- (max.dx - min.dx) * times + 1

    # Adjested density.
    d2 <- density(posterior.x,
      bw = MASS::width.SJ(posterior.x, method = "dpi"),
      n = n.dx, from = min.dx, to = max.dx)
  }

  return(d2)
}

WrapperDensity <- function(result.d, is.2d = FALSE){
# Wrapper function returns MDE and its probability.
#
# args:
#  result.d: Returned result from density() or kde2d().
#  is.2d: If TRUE, function is applied to 2d-KDE result.
#
# returns:
#  An array contains [c(x, y), c(mde, marginal.probability)], if is.2d = TRUE.
#  A vector contains c(mde, marginal.probability), if is.2d = FALSE.
#
# functions:
#  CalcPostDensity
  if(is.2d){
    mde.ind <- which(result.d$z == max(result.d$z), arr.ind = TRUE)
    mde <- c(result.d$x[mde.ind[1]], result.d$y[mde.ind[2]])
    mde.prob.x <- sum(result.d$z[mde.ind[1], ]) / sum(result.d$z)
    mde.prob.y <- sum(result.d$z[ , mde.ind[2]]) / sum(result.d$z)
    estimator <- cbind(mde, c(mde.prob.x, mde.prob.y))
  }else{
    mde.ind <- which(result.d$y == max(result.d$y))
    mde <- result.d$x[mde.ind]
    mde.prob <- result.d$y[mde.ind] / sum(result.d$y)
    estimator <- c(mde, mde.prob)
  }

  return(estimator)
}

WrapperABCPRC <- function(age.ref = seq(0, 10, 1.0), prior,
  age.obs, d15N, female.mean, num.particle = 10000, form = "parabolic",
  tolerances = c(2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0)){
# Set threshold (optimized value) and perform ABC-PRC
#  for a given observed dataset and prior distributions.
#
# args:
#  age.ref: A vector of reference ages (unit time).
#  prior: Hyper parameters in the order of c(mu.t1, sigma.t1,
#   mu.t2, sigma.t2, mu.enrich, sigma.enrich, mu.wnfood, sigma.wnfood,
#   mu.epsilon, sigma.epsilon).
#   mu: mean, sigma: 1SD.
#   t1: Age at the beginning of weaning.
#   t2: Age at the completion of weaning.
#   enrich: Enrichment factor between mother and infant.
#   wnfood: d15N value of weaning food.
#   epsilon: Hyper parameter (1SD) for the individual epsilon.
#  age.obs: A vector of observed ages
#   with which we want to calculate the posterior.
#  d15N: A vector of observed d15N values for the nonadults.
#  female.mean: Mean d15N value of adult females.
#  num.particle: The number of extracted parameters included in one particle.
#  form: The form of change in d15N values during weaning.
#
# returns:
#  A list of target theta particle and its weight.
#   
# functions:
#  SubtractAgeResidue
#  ArrangeTurnoverRate
#  OptimizePar
#  PerformABCPRC
  age.reference <- SubtractAgeResidue(age.ref)
  age.observed <- SubtractAgeResidue(age.obs)
 
  age.ref <- age.reference[ , 1]
  residue.ref <- age.reference[ , 2]
  subtract.ref <- age.reference[ , 3]

  age.obs <- age.observed[ , 1]
  residue.obs <- age.observed[ , 2]
  subtract.obs <- age.observed[ , 3]

  turnover.ref <- ArrangeColTurnover(age.ref, subtract.ref)
  turnover.obs <- ArrangeColTurnover(age.obs, subtract.obs)

  par.initial <- c(prior[c(1, 3, 5, 7)])
  par.opt <- OptimizePar(par.initial = par.initial,
    age.ref = age.ref,
    residue.ref = residue.ref,
    subtract.ref = subtract.ref,
    age.obs = age.obs,
    residue.obs = residue.obs,
    subtract.obs = subtract.obs,
    female.mean = female.mean,
    epsilon = 0,
    turnover.ref = turnover.ref,
    turnover.obs = turnover.obs,
    d15N = d15N,
    form = form)

  threshold <- par.opt$value
#  tolerances <- tolerances 

  prc <- PerformABCPRC(threshold, tolerances, num.particle,
    prior[1], prior[2], prior[3], prior[4],
    prior[5], prior[6], prior[7], prior[8], 0, 1,
    age.ref = age.ref, residue.ref = residue.ref,
    subtract.ref = subtract.ref,
    age.obs = age.obs, residue.obs = residue.obs,
    subtract.obs = subtract.obs, female.mean = female.mean,
    turnover.ref = turnover.ref,
    turnover.obs, d15N = d15N, form = form)

  return(prc)
}

WARN <- function(age, d15N, female.mean, female.sd = NA,
  prior = c(0.5, 3, 3, 3, 1.9, 0.9, female.mean, 3), 
  num.particle = 10000, form = "parabolic",
  tolerances = c(2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0)){
# Wrapper function for weaning ages reconstruction using ABC-SMC-PRC.
#  Weaning Ages Reconstruction from Nitrogen isotope ratios.
#
# args:
#  age: A vector of observed ages for the subadults.
#  d15N: A vector of observed d15N values for the subadults.
#  female.mean, female.sd: Female mean and SD d15N values.
#  prior: Hyper parameters in the order of c(mu.t1, sigma.t1,
#   mu.t2, sigma.t2, mu.enrich, sigma.enrich, mu.wnfood, sigma.wnfood,
#   mu.epsilon, sigma.epsilon).
#   mu: mean, sigma: 1SD.
#   t1: Age at the beginning of weaning.
#   t2: Age at the completion of weaning.
#   enrich: Enrichment factor between mother and infant.
#   wnfood: d15N value of weaning food.
#  num.particle: The number of extracted parameters included in one particle.
#  form: The form of change in d15N values during weaning.
#
# returns:
#  A list of result which will be attributed class "warn".
#
# functions:
#  WrapperABCPRC
#  DrawAll
#  DrawContour
#  WrapperDensity
#  CalcDistanceMDE
  # ABC-PRC.
  prior <- c(prior, c(0, 1)) # Prior parameters for epsilon.
  posterior <- WrapperABCPRC(age.ref = seq(0, 10, 1.0), prior = prior,
    age.obs = age, d15N = d15N, female.mean = female.mean,
    num.particle = num.particle, form = form,
    tolerances = tolerances)$theta

  # KDE.
  kde.age <- CalcPostDensity(posterior.x = posterior[ , 1],
    posterior.y = posterior[ , 2],
    is.2d = TRUE, decimal = 1)
  kde.enrich <- CalcPostDensity(posterior.x = posterior[ , 3],
    is.2d = FALSE, decimal = 1)
  kde.wnfood <- CalcPostDensity(posterior.x = posterior[ , 4],
    is.2d = FALSE, decimal = 1)

  kdes <- list(age = kde.age, enrich = kde.enrich, wnfood = kde.wnfood)

  # Calculate KDE-MDEs and its probabilities.
  estimator <- rbind(WrapperDensity(result.d = kde.age, is.2d = TRUE),
    WrapperDensity(result.d = kde.enrich, is.2d = FALSE),
    WrapperDensity(result.d = kde.wnfood, is.2d = FALSE))
  rownames(estimator) <- c("t1", "t2", "enrich", "wnfood")
  colnames(estimator) <- c("mde", "probability")

  # Calculate probability of the combination of weaning ages.
  prob.2d.age <- CalcProb2D(kde = kde.age,
    range.x = rep(estimator[1, 1], 2),
    range.y = rep(estimator[2, 1], 2))$prob

  # Calculate mean squared distance between observed and modeled data.
  dist.mde <- CalcDistanceMDE(par.mde = estimator[ , 1],
    age.obs = age,
    d15N = d15N,
    female.mean = female.mean,
    form = form)

  return(list(mde = estimator,
    prob.2d.age = prob.2d.age,
    dist.mde = dist.mde,
    kde.age = kde.age,
    kde.enrich = kde.enrich,
    kde.wnfood = kde.wnfood,
    posterior = posterior[ , -5],
    age = age,
    d15N = d15N,
    female.mean = female.mean,
    female.sd = female.sd,
    prior = prior,
    particle = num.particle,
    form = form))
}

# warnProb --------------------
CalcProb1D <- function(kde, range.x){
# Calculate KDE-probability for given weaning parameter ranges.
#
# args:
#  kde: Product of density(). Usually a list contains x and y.
#  range.x: A range of x with which we want to calculate probability.
#   Fractional point lower than e-002 is omitted.
#
# returns:
#  A list contains $kde, $prob, and $range c(from.x, to.x).
#
# note:
#  This function works for "enrichment" and "wnfood".
  # Edit the ranges.
  range.x1 <- round(range.x, 1)
  range.x2 <- range.x1
  kde.x <- round(kde$x, 1)

  range.x2[1] <- ifelse(range.x2[1] < min(kde.x), min(kde.x), range.x2[1])
  range.x2[2] <- ifelse(range.x2[2] > max(kde.x), max(kde.x), range.x2[2])

  # Calculate the indexes of the ranges.
  ind.x <- c(which(kde.x == range.x2[1]), which(kde.x == range.x2[2]))

  # Calculate the probability of targetted ranges.
  target <- sum(kde$y[((ind.x[1]):(ind.x[2]))])
  prob <- target / sum(kde$y)

  result <- list(kde = kde,
    probability = prob,
    range = c(from.x = range.x1[1], to.x = range.x1[2]))

  return(result)
}

CalcProb2D <- function(kde, range.x, range.y){
# Calculate KDE-probability for given weaning age ranges.
#
# args:
#  kde: Product of kde2d(). Usually a list contains x, y, z.
#  range: Weaning age range with which we want to calculate probability.
#   Fractional point lower than e-002 is omitted.
#   x: t1, and y: t2.
#
# returns:
#  A list contains $kde, $prob, and $range c(from.x, to.x, from.y, to.y).
#
# note:
#  This function works for "age".
  # Edit the ranges.
  range.x1 <- round(range.x, 1)
  range.y1 <- round(range.y, 1)
  range.x2 <- range.x1
  range.y2 <- range.y1
  kde.x <- round(kde$x, 1)
  kde.y <- round(kde$y, 1)

  range.x2[1] <- ifelse(range.x2[1] < min(kde.x), min(kde.x), range.x2[1])
  range.x2[2] <- ifelse(range.x2[2] > max(kde.x), max(kde.x), range.x2[2])
  range.y2[1] <- ifelse(range.y2[1] < min(kde.y), min(kde.y), range.y2[1])
  range.y2[2] <- ifelse(range.y2[2] > max(kde.y), max(kde.y), range.y2[2])

  # Calculate the indexes of the ranges.
  ind.x <- c(which(kde.x == range.x2[1]), which(kde.x == range.x2[2]))
  ind.y <- c(which(kde.y == range.y2[1]), which(kde.y == range.y2[2]))

  # Calculate the probability of targetted ranges.
  target <- sum(kde$z[((ind.x[1]):(ind.x[2])), (ind.y[1]):(ind.y[2])])
  prob <- target / sum(kde$z)

  result <- list(kde = kde,
    probability = prob,
    range = c(from.x = range.x1[1], to.x = range.x1[2],
      from.y = range.y1[1], to.y = range.y1[2]))

  return(result)
}

# graphics --------------------
DrawProb1D <- function(kde, range.x1 = NA, mde = NA,
  is.legend = TRUE, is.prior = FALSE, hyper.par = NA, ...){
# Draw KDE-probability with rectangle showing a given weaning parameter range.
#
# args:
#  kde: Product of density(). Usually a list contains x and y.
#  range.x: A range of x with which we want to calculate probability.
#   Fractional point lower than e-002 is omitted.
#  mde: MDE, which is indicated on the figure.
#  is.legend: If FALSE, a default legend is not drawn.
#  is.prior: If TRUE, prior distribution is also drawn.
#  ...: Additional arguments are passed to plot().
#
# returns:
#  A figure indicating KDE-probability.
  x <- kde$x
  probability <- kde$y / sum(kde$y)
  Probability <- probability

  # Probability.
  plot(x, Probability, type = "l", ...)

  # Prior.
  if(is.prior){
    range.prior <- seq(range.x1[1] - 10, range.x1[2] + 10, 0.02)
    prob.prior <- dnorm(range.prior, mean = hyper.par[1], sd = hyper.par[2])
    points(range.prior, prob.prior / 10, type = "l", lty = "dotted")
  }

  # MDE.
  points(mde, max(probability), pch = 16, col = "red")

  # Range.
  rect(range.x1[1] - 0.05, 0 - max(probability) / 50,
    range.x1[2] + 0.05, max(probability) * 51 / 50,
    border = "red")

  # Legend.
  if(is.legend){
    pch.l <- c(NA, NA, 19, NA)
    lty.l <- c("solid", "dotted", NA, "solid")
    col.l <- c("black", "black", "red", "red")
    legend.l <- c("Posterior", "Prior", "MDE", "Range")

    if(is.prior){
      ind.l <- 1:4
    }else{
      ind.l <- c(1, 3, 4)
    }

    legend(min(x), max(probability),
      legend = legend.l[ind.l],
      pch = pch.l[ind.l],
      lty = lty.l[ind.l],
      col = col.l[ind.l],
      cex = 0.7, xjust = 0, yjust = 1)
  }

}

DrawProb2D <- function(kde, range.x1 = NA, range.y1 = NA, mde = c(NA, NA),
  is.legend = TRUE, is.contour = TRUE, is.image = FALSE, ...){
# Draw KDE-probability with rectangle showing given weaning age ranges.
#
# args:
#  kde: Product of kde2d(). Usually a list contains x, y, z.
#  range: Lower and upper renges with which we want to indicate rectangle.
#   Fractional point lower than e-002 is omitted.
#  mde: MDE, which is indicated on the figure, c(t1, t2).
#  is.image: If TRUE, results are also descrived by image().
#  is.legend: If TRUE, a legend is drawn.
#  is.legend: If FALSE, a default legend is not drawn.
#  ...: Additional arguments are passed to plot().
#
# returns:
#  A figure indicating KDE-probability.
  t1 <- c(0, max(kde$x))
  t2 <- c(0, max(kde$y))
  plot(t1, t2, type = "n", ...)

  # image().
  if(is.image){
    image(kde, col =  gray((100:0)/100), add = TRUE)
  }

  # contour().
  if(is.contour){
    mde.freq <- max(kde$z)
    mde.levels <- (1:5 * (mde.freq / 5)) - mde.freq / 10
    contour(kde, add = TRUE, levels = mde.levels,
      labels = c(0.1, 0.3, 0.5, 0.7, 0.9))
  }

  # Line dividing t1 and t2.
  points(c(-10, 20), c(-10, 20), type = "l")

  # MDE.
  points(mde[1], mde[2], pch = 16, col = "red")

  # Range.
  rect(range.x1[1] - 0.05, range.y1[1] - 0.05,
    range.x1[2] + 0.05, range.y1[2] + 0.05,
    border = "red")

  # Legend.
  if(is.legend){
    legend(diff(t1), 0, legend = c("MDE", "Range"),
      pch = c(19, NA), lty = c(NA, "solid"), col = "red",
      cex = 0.7, xjust = 1, yjust = 0)
  }
}

DrawMDE <- function(par.mde, d15N, age,
  female.mean = NA, female.sd = 0,
  form = "parabolic",
  hline.female = TRUE,
  hline.adult = FALSE,
  adult.mean = NA, adult.sd = 0, 
  is.legend = TRUE,
  is.female = TRUE,
  cex = 1, ...){
# Draw the KDE-MDE results.
#
# args:
#  par.mde: A vector of the MDEs: c(t1, t2, enrich, wnfood).
#  d15N: A vector of observed d15N values for the nonadults.
#  age: A vector of observed ages for the subadults.
#  female.mean: Female mean d15N values.
#  female.sd: SD of female d15N values.
#  form: The form of change in d15N values during weaning.
#  hline: If true 1SD ranges are indicated by horizontal line.
#  adult.mean: Adult mean d15N values.
#  adult.sd: SD of adult d15N values.
#  is.legend: If TRUE, a legend is drawn.
#  is.female: If FALSE, female mean is not drawn.
#  ...: Arguments passed to plot().
#
# returns:
#  A figure of the KDE-MDE result.
#
# function:
#  ArrangeColTurnover
#  CalcNonmilkProp
#  IntNonmilkProp
#  CalcRefBone
  xmax <- ceiling(max(age) * 2) / 2 + 1
  age.reference <- SubtractAgeResidue(seq(0, xmax - 1, 0.5))
  age.ref <- age.reference[ , 1]
  residue <- age.reference[ , 2]
  subtract <- age.reference[ , 3]
  turnover.ref <- ArrangeColTurnover(age.ref, subtract)

  t1 <- par.mde[1]
  t2 <- par.mde[2]
  enrich <- par.mde[3]
  n.milk <- female.mean + enrich
  n.wnfood <- par.mde[4]

  # Change in the diet d15N value.
	before.weaning <- c(0, t1)
	during.weaning <- seq(t1, t2, length = 100)
	after.weaning <- c(t2, xmax - 1)
	model.age <- c(before.weaning, during.weaning, after.weaning)

  wnfood.mod <- mapply(CalcNonmilkProp,
    age = during.weaning, t1 = t1, t2 = t2, form = form)
  model.diets <- (1 - wnfood.mod) * (n.milk - n.wnfood) + n.wnfood
  model.diets <- c(n.milk, n.milk, model.diets, n.wnfood, n.wnfood)

  # Change in the bone d15N value.
  inted.wnfood.ref <- IntNonmilkProp(age.ref, residue, subtract, t1, t2, form)

  ref.bones <- CalcRefBone(age.ref, residue, subtract,
    inted.wnfood.ref, enrich, female.mean, n.wnfood, turnover.ref)

  # Draw the figure.
  Age <- c(0, xmax)
  delta_15N <- c(min(d15N) - 1, max(d15N) + 1)
  plot(Age, delta_15N, type = "n", ...)

  # Subadults.
	points(age, d15N, pch = 23, cex = cex)

	# Measured d15N mean and SD of females.
  female.mean <- ifelse(is.female, female.mean, NA)
  female.sd <- ifelse(is.female, female.sd, NA)

  female.age <- xmax - 0.4
	arrows(female.age, (female.mean + female.sd),
    female.age, (female.mean - female.sd),
    angle = 90, length = 0.025, code = 3)
	points(female.age, female.mean, pch = 19, col = "white", cex = cex)
	points(female.age, female.mean, pch = 21, cex = cex)
  
	# Measured d15N mean and SD of adults.
  adult.age <- xmax
	arrows(adult.age, (adult.mean + adult.sd),
    adult.age, (adult.mean - adult.sd),
    angle = 90, length = 0.025, code = 3)
	points(adult.age, adult.mean, pch = 4, cex = cex)

  # Horizontal line of adults.
  if(hline.adult){
    lines(c(-10, 100), rep(adult.mean + adult.sd, 2),
      type = "l", lty = "dotted", lwd = 0.5)
    lines(c(-10, 100), rep(adult.mean - adult.sd, 2),
      type = "l", lty = "dotted", lwd = 0.5)
  }
  if(hline.female){
    lines(c(-10, 100), rep(female.mean + female.sd, 2),
      type = "l", lty = "dotted", lwd = 0.5)
    lines(c(-10, 100), rep(female.mean - female.sd, 2),
      type = "l", lty = "dotted", lwd = 0.5)
  }

	# Optimized d15N of diet.
	points(model.age, model.diets, type = "l", cex = cex)

	# Optimized d15N of bone.
	points(age.ref, ref.bones, pch = 18, cex = cex)

  # Legend.
  if(is.legend){
    pch.l <- c(23, 18, NA, 4, 21)
    lty.l <- c(NA, NA, "solid", NA, NA)
    legend.l <- c("Measured d15N", "Modeled d15N", "Modeled diet",
      "Total adult", "Adult female")
    ind.l <- 1:5
    if(is.na(adult.mean)){
      ind.l <- ind.l[-(which(ind.l == 4))]
    }
    if(is.na(female.mean)){
      ind.l <- ind.l[-(which(ind.l == 5))]
    }

    legend(xmax, delta_15N[2], legend = legend.l[ind.l],
    pch = pch.l[ind.l], lty = lty.l[ind.l], cex = 0.7,
    xjust = 1, yjust = 1)
  }
}

# warnCI --------------------
CalcCI2D <- function(kde, mde.x, mde.y,
  threshold = 0.95){
# Calcurate KDE-probability  and its range for given CI from kde data.
#  This function works for "age" of the weaning parameters.
#
# args:
#  kde: Product of kde2d(). Usually a list contains x, y, z.
#  mde.x/y: Weaning age of MDE estimates.
#   Fractional point lower than e-002 is omitted.
#   x: t1, and y: t2.
#  threshold: A scalar indicating threshold of CI. From 0 to 1.
#
# returns:
#  A list contains $kde, $probability,
#   and $range c(from.x, to.x, from.y, to.y).
#
# depend:
#  CalcProb2D
#
  # Warnings.
  if(threshold > 1){
    stop(message="c(0, 1) for the range of 'threshold'")
  }
  # Variables
  nega <- c(-0.1, 0)
  posi <- c(0, 0.1)

  # Initial values.
  range.new.x <- rep(mde.x, 2)
  range.new.y <- rep(mde.y, 2)
  prob.mde <- CalcProb2D(
    kde = kde,
    range.x = range.new.x,
    range.y = range.new.y)
  probability.new <- prob.mde$probability

  # If MDE is Ok, return MDE.
  if(probability.new > threshold){
    result <- list(
#      kde = kde,
      probability = probability.new,
      range = c(
        prob.mde$range[1],
        prob.mde$range[2],
        prob.mde$range[3],
        prob.mde$range[4]))
    return(result)
  }

  # Searching
  while(probability.new < threshold){
    names(range.new.x) <- NULL
    names(range.new.y) <- NULL

    # Searching for (expanding to) the 4 directions.
    prob.dash <- list(
      prob.nx = CalcProb2D(
        kde = kde,
        range.x = range.new.x + nega,
        range.y = range.new.y),
      prob.px = CalcProb2D(
        kde = kde,
        range.x = range.new.x + posi,
        range.y = range.new.y),
      prob.xn = CalcProb2D(
        kde = kde,
        range.x = range.new.x,
        range.y = range.new.y + nega),
      prob.xp = CalcProb2D(
        kde = kde,
        range.x = range.new.x,
        range.y = range.new.y + posi))

    # Select new range with max robability.
    probability.dash <- c(
      prob.dash[[1]]$probability,
      prob.dash[[2]]$probability,
      prob.dash[[3]]$probability,
      prob.dash[[4]]$probability)
    ind.max <- match(max(probability.dash), probability.dash)
    prob.new <- prob.dash[[ind.max]]

    # Variables for the next searching.
    probability.new <- prob.new$probability
    range.new.x <- prob.new$range[1:2]
    range.new.y <- prob.new$range[3:4]
  }

  # Results.
  result <- list(
#    kde = kde,
    probability = probability.new,
    range = c(range.new.x, range.new.y))

  return(result)
}


CalcCI1D <- function(kde, mde.x,
  threshold = 0.95){
# Calcurate KDE-probability  and its range for given CI from kde data.
#  This function works for "enrich" and "wnfood" of the weaning parameters.
#
# args:
#  kde: Product of kde2d(). Usually a list contains x, y, z.
#  mde.x/y: Weaning age of MDE estimate.
#   Fractional point lower than e-002 is omitted.
#  threshold: A scalar indicating threshold of CI. From 0 to 1.
#
# returns:
#  A list contains $kde, $probability,
#   and $range c(from.x, to.x).
#
# depend:
#  CalcProb1D
#
  # Warnings.
  if(threshold > 1){
    stop(message="c(0, 1) for the range of 'threshold'")
  }

  # Variables
  nega <- c(-0.1, 0)
  posi <- c(0, 0.1)

  # Initial values.
  range.new.x <- rep(mde.x, 2)
  prob.mde <- CalcProb1D(
    kde = kde,
    range.x = range.new.x)
  probability.new <- prob.mde$probability

  # If MDE is Ok, return MDE.
  if(probability.new > threshold){
    result <- list(
#      kde = kde,
      probability = probability.new,
      range = c(
        prob.mde$range[1],
        prob.mde$range[2]))
    return(result)
  }

  # Searching
  while(probability.new < threshold){
    names(range.new.x) <- NULL

    # Searching for (expanding to) the 4 directions.
    prob.dash <- list(
      prob.n = CalcProb1D(
        kde = kde,
        range.x = range.new.x + nega),
      prob.p = CalcProb1D(
        kde = kde,
        range.x = range.new.x + posi))

    # Select new range with max robability.
    probability.dash <- c(
      prob.dash[[1]]$probability,
      prob.dash[[2]]$probability)
    ind.max <- match(max(probability.dash), probability.dash)
    prob.new <- prob.dash[[ind.max]]

    # Variables for the next searching.
    probability.new <- prob.new$probability
    range.new.x <- prob.new$range
  }

  # Results.
  result <- list(
#    kde = kde,
    probability = probability.new,
    range = range.new.x)

  return(result)
}

# ==============================
# NOTE ----------
# ABC-SMC-PRC algorithm was quoted from
#  Sisson, S. Fan, Y. Tanaka, Mark M. 2007, 2009.
#  "Sequential Monte Carlo without likelihoods."
#  PNAS 104:1760-1765.
# ==============================
# END
