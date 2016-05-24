###
# Varying model formulations (other than formulas). Using and subsetting the
# 'varying' variables.
###
library(nlme)
library(MuMIn)

# from example(corSpher)
fm1BW.lme <- lme(weight ~ Time * Diet, BodyWeight, random = ~ Time)


# generate model selection table:
fm1BW.dd <- dredge(fm1BW.lme,
	# fix all terms in all models:
	fixed = TRUE,
	varying = list(
		# vary correlation structure:
		correlation = alist(exp = corExp(form = ~ Time),
							spher = corSpher(form = ~ Time),
							NULL ),
		# vary heteroscedasticity structure:
		weights = alist(vPower = varPower(),
						none = NULL )
	),
	# additional constraint (regardless of whether it makes sense or not): 
	#  include either heteroscedasticity or correlation structure (but not both).
	# Note use of 'is.null' for unnamed item, and "none" when named.
	subset = xor(is.null(V(correlation)), V(weights) == "none"),
	# global model was fitted with method = "REML" (the default), but for model 
	# selection we use AICc of a ML model. This additional argument is passed to
	# AICc.
	REML = FALSE)


print(fm1BW.dd)