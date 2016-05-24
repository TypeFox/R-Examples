#     This file is part of FAiR, a program to conduct Factor Analysis in R
#     Copyright 2008 Benjamin King Goodrich
#
#     FAiR is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     FAiR is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Affero General Public License for more details.
#
#     You should have received a copy of the GNU Affero General Public License
#     along with FAiR.  If not, see <http://www.gnu.org/licenses/>.

## These are new S4 generics that probably only make sense within FAiR.
## To add a model to FAiR, methods for some or all of these generics are needed.

## make_manifest is a constructor for an object that inherits from class manifest
setGeneric("make_manifest", def = 
function(x, data, covmat, ...) standardGeneric("make_manifest"),
useAsDefault = FALSE, package = "FAiR"

# showMethods("make_manifest")
# Function: make_manifest (package FAiR)
# x="data.frame", data="missing", covmat="missing"
# x="formula", data="data.frame", covmat="missing"
# x="matrix", data="missing", covmat="missing"
# x="missing", data="data.frame", covmat="missing"
# x="missing", data="matrix", covmat="missing"
# x="missing", data="missing", covmat="CovMcd"
# x="missing", data="missing", covmat="hetcor"
# x="missing", data="missing", covmat="list"
# x="missing", data="missing", covmat="matrix"

)

## make_restrictions is a constructor for an object that inherits from class restrictions
# setGeneric("make_restrictions", def = 
# function(manifest, ...) standardGeneric("make_restrictions"),
# useAsDefault = FALSE, package = "FAiR"
# 
# # showMethods("make_restrictions")
# # Function: make_restrictions (package FAiR)
# # manifest="manifest.basic"
# 
# )

setGeneric("make_restrictions", def = 
function(manifest, Omega, beta, Phi, Delta, Xi, ...) standardGeneric("make_restrictions"),
useAsDefault = FALSE, package = "FAiR"

)

## restrictions2model fills the free parameters within a restrictions object
setGeneric("restrictions2model", def = 
function(par, restrictions, manifest, lower, ...) standardGeneric("restrictions2model"),
useAsDefault = FALSE, signature = c("restrictions", "manifest"), package = "FAiR"

# showMethods("restrictions2model")
# Function: restrictions2model (package FAiR)
# restrictions="restrictions.1storder.EFA", manifest="manifest.basic"
# restrictions="restrictions.1storder", manifest="manifest.basic"
# restrictions="restrictions.2ndorder", manifest="manifest.basic"
# restrictions="restrictions.general", manifest="manifest.basic"
# restrictions="restrictions.independent", manifest="manifest.basic"
# restrictions="restrictions.orthonormal", manifest="manifest.basic"

)

## fitS4 produces a vector of criteria to be lexically minimized
setGeneric("fitS4", def = 
function(par, restrictions, manifest, lower, ...) standardGeneric("fitS4"), 
useAsDefault = FALSE, signature = c("restrictions", "manifest"), package = "FAiR"

# showMethods("fitS4")
# Function: fitS4 (package FAiR)
# restrictions="restrictions.factanal", manifest="manifest.basic"
# restrictions="restrictions.independent", manifest="manifest.basic"
# restrictions="restrictions", manifest="manifest.basic"

)

## bfgs_fitS4 produces a scalar criterion, usually the last criterion from fitS4
setGeneric("bfgs_fitS4",  def = 
function(par, restrictions, manifest, helper, lower, ...) standardGeneric("bfgs_fitS4"), 
useAsDefault = FALSE, signature = c("restrictions", "manifest"), package = "FAiR"

# showMethods("bfgs_fitS4")
# Function: bfgs_fitS4 (package FAiR)
# restrictions="restrictions", manifest="manifest.basic"

)

## gr_fitS4 produces the gradient of bfgs_fitS4
setGeneric("gr_fitS4", def = 
function(par, restrictions, manifest, helper, lower, ...) standardGeneric("gr_fitS4"),
useAsDefault = FALSE, signature = c("restrictions", "manifest"), package = "FAiR"

# showMethods("gr_fitS4")
# Function: gr_fitS4 (package FAiR)
# restrictions="restrictions.1storder", manifest="manifest.basic"
# restrictions="restrictions.2ndorder", manifest="manifest.basic"
# restrictions="restrictions.factanal", manifest="manifest.basic"
# restrictions="restrictions.general", manifest="manifest.basic"
# restrictions="restrictions", manifest="manifest.basic"
# restrictions="restrictions.orthonormal", manifest="manifest.basic"

)

## bfgs_helpS4 provides helper material (in a list) to bfgs_fitS4 and gr_fitS4
setGeneric("bfgs_helpS4", def= function(initial, restrictions, manifest, done, lower, ...)
                standardGeneric("bfgs_helpS4"), useAsDefault = FALSE, 
		signature = c("restrictions", "manifest"), package = "FAiR",

# showMethods("bfgs_helpS4")
# Function: bfgs_helpS4 (package FAiR)
# restrictions="restrictions", manifest="manifest.basic"

)

## create_FAobject is a constructor for objects that are or inherit from class "FA"
## this function is called at the end of Factanal() to bundle all the estimates, etc.
setGeneric("create_FAobject", def = function(restrictions, manifest, opt, call, scores,
	  			lower, analytic, ...) standardGeneric("create_FAobject"),
useAsDefault = FALSE, signature = c("restrictions", "manifest"), package = "FAiR"

# showMethods("create_FAobject")
# Function: create_FAobject (package FAiR)
# restrictions="restrictions.1storder", manifest="manifest.basic"
# restrictions="restrictions.2ndorder", manifest="manifest.basic"
# restrictions="restrictions.factanal", manifest="manifest.basic"
# restrictions="restrictions.general", manifest="manifest.basic"
# restrictions="restrictions.independent", manifest="manifest.basic"
# restrictions="restrictions", manifest="manifest.basic"
# restrictions="restrictions.orthonormal", manifest="manifest.basic"

)

## create_start creates a matrix of starting values and can be called from the prompt.
## However, it is typically called from within Factanal() before the call to genoud()
setGeneric("create_start", def=function(number, start, restrictions, manifest, reject,...)
standardGeneric("create_start"), useAsDefault = TRUE, 
signature = c("restrictions", "manifest"), package = "FAiR"

# showMethods("create_start")
# Function: create_start (package FAiR)
# restrictions="restrictions.1storder", manifest="manifest.basic"
# restrictions="restrictions.2ndorder", manifest="manifest.basic"
# restrictions="restrictions.factanal", manifest="manifest.basic"
# restrictions="restrictions.general", manifest="manifest.basic"
# restrictions="restrictions", manifest="manifest.basic"
# restrictions="restrictions.orthonormal", manifest="manifest.basic"

)

## restrictions2draws creates a matrix of simulated parameters a la King (1991) and would
## typically be called indirectly by FA2draws() at the command prompt
setGeneric("restrictions2draws", def = 
function(restrictions, manifest, vcov, nsim, covariances, standardized, ...)
	standardGeneric("restrictions2draws"), useAsDefault = TRUE, 
	signature = c("restrictions", "manifest"), package = "FAiR"

# showMethods("restrictions2draws")
# Function: restrictions2draws (package FAiR)
# restrictions="restrictions.1storder.EFA", manifest="manifest.basic"
# restrictions="restrictions.1storder", manifest="manifest.basic"
# restrictions="restrictions.2ndorder", manifest="manifest.basic"
# restrictions="restrictions.general", manifest="manifest.basic"
# restrictions="restrictions.orthonormal", manifest="manifest.basic"

)

## convert from restrictions to format for reticular action model
setGeneric("restrictions2RAM", useAsDefault = TRUE, signature = "restrictions", def = 
function(restrictions, ...) standardGeneric("restrictions2RAM"), package = "FAiR"

# showMethods("restrictions2RAM")
# Function: restrictions2RAM (package FAiR)
# restrictions="restrictions.1storder"
# restrictions="restrictions.2ndorder"
# restrictions="restrictions.general"
# restrictions="restrictions.orthonormal"

)

## cormat is a "get" function for a correlation matrix and could be called by the user
setGeneric("cormat", def = function(object, ...) standardGeneric("cormat"),
		useAsDefault = TRUE, signature = "object", package = "FAiR"

# showMethods("cormat")
# Function: cormat (package FAiR)
# object="FA"
# object="FA.2ndorder"
# object="manifest.basic"
# object="parameter.cormat"
# object="restrictions"
# object="restrictions.2ndorder"
# object="restrictions.factanal"
# object="restrictions.orthonormal"

)

## uniquenesses is a "get" function for the uniquenesses and could be called by the user
setGeneric("uniquenesses", def = function(object, ...) standardGeneric("uniquenesses"),
		useAsDefault = TRUE, signature = "object", package = "FAiR"

# showMethods("uniquenesses")
# Function: uniquenesses (package FAiR)
# object="FA"
# object="FA.general"
# object="restrictions"
# object="restrictions.factanal"

)

## make_parameter fills in an object of parameter class with par (theta)
setGeneric("make_parameter", signature = "object", def = function(par, object, ...) 
standardGeneric("make_parameter"), useAsDefault = TRUE, package = "FAiR"

# showMethods("make_parameter")
# Function: make_parameter (package FAiR)
# object="parameter.coef"
# object="parameter.coef.nl"
# object="parameter.coef.SEFA"
# object="parameter.coef.SEFA.nl"
# object="parameter.cormat"
# object="parameter.scale"

)

## loadings is a "get" function for the "loadings" and could be called by the user
setGeneric("loadings", def = function(x, ...) standardGeneric("loadings")

# Function: loadings (package FAiR)
# x="ANY"
# x="FA"
# x="FA.general"
# x="restrictions"
# x="restrictions.factanal"
# x="restrictions.general"

)

## in FAiR, this function does a Thurstone plot, only defined for x = "FA"
setGeneric("pairs", def = function(x, ...) standardGeneric("pairs"))

## in FAiR, this function yields a covariance matrix as produced by the model
setGeneric("fitted", def = function(object, ...) standardGeneric("fitted")

# Function: fitted (package stats)
# object="ANY"
# object="FA"
# object="restrictions"
# object="restrictions.factanal"
# object="restrictions.independent"
# object="restrictions.orthonormal"

)

## in FAiR, this function yields residual correlations or covariances (defined for "FA")
setGeneric("residuals", def = function(object, ...) standardGeneric("residuals"))

## in FAiR, this function yields standardized residual covariances (defined for "FA")
setGeneric("rstandard", def = function(model, ...) standardGeneric("rstandard"))

## in FAiR, this function yields (approximate) weights used (defined for "FA")
setGeneric("weights", def = function(object, ...) standardGeneric("weights"))

## in FAiR, this function yields some function of the weighted residuals(defined for "FA")
setGeneric("influence", def = function(model, ...) standardGeneric("influence"))

## gets degrees of freedom (defined for "FA" and "restrictions")
setGeneric("df.residual", def = function(object, ...) standardGeneric("df.residual"))

## in FAiR, gets the value of the discrepancy function (from a "FA" object)
setGeneric("deviance", def = function(object, ...) standardGeneric("deviance"))

## in FAiR, gets the sample covariance or correlation matrix
setGeneric("model.matrix", def = function(object, ...) standardGeneric("model.matrix")

# Function: model.matrix (package stats)
# object="ANY"
# object="FA"
# object="manifest.basic"

)

## in FAiR, draws values for the covariance matrix
setGeneric("simulate", def = function(object, nsim = 1, seed = NULL, ...)
standardGeneric("simulate"))
