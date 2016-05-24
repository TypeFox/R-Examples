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

## This file defines formal S4 classes used by FAiR

## NOTE: This file is meant to be read with 90 columns and 8 space tabs


## this critical function is for executing the mapping rule (for exact zeros) in a SEFA
## this function is not a class and is not even a generic, but is needed for a prototype
mapping_rule <- 
function(coefs, cormat, zeros = rep(ncol(coefs), ncol(coefs)),
	row_complexity = NA_integer_, quasi_Yates = FALSE, weak_Thurstone = FALSE, 
	Butler = FALSE, viral = FALSE, mve = FALSE, communality = FALSE) {

	# Make reference structure matrix
	D_inv <- sqrt(diag(solve(cormat)))
	Upsilon <- t(t(coefs) / D_inv)

	if(!is.na(row_complexity[1])) {            # some zeros by rows
		if(length(row_complexity) == 1) {  # all rows have the same complexity
			Upsilon <- t(apply(abs(Upsilon), 1, FUN = function(x) {
				ranks <- order(x, decreasing = TRUE)
				x[ranks > row_complexity] <- 0 # zero out some cells
				return(x)
			}))
		}
		else {  # rows have different complexities
			joined <- cbind(row_complexity, abs(Upsilon))
			Upsilon <- t(apply(joined, 1, FUN = function(x) {
				row_complexity <- x[1]
				x <- x[-1]
				ranks <- order(x, decreasing = TRUE)
				x[ranks > row_complexity] <- 0 # zero out some cells
				return(x)
			}))
		}
		coefs <- t(t(Upsilon) * D_inv)
	} # have to keep going to check zeros by column condition

	else if(weak_Thurstone) { # r zeros per factor, all of complexity r - 1
		Upsilon.abs <- abs(Upsilon)
		orderer <- order(c(Upsilon.abs))
		changed <- rep(0, ncol(coefs))
		for(i in 1:length(orderer)) { ## explain logic here
			dims <- which(coefs == coefs[orderer[i]], arr.ind = TRUE)
			if(changed[dims[1,2]] == ncol(coefs)) next
			if(any(!coefs[dims[1,1], ])) next
			coefs[orderer[i]] <- 0
			changed[dims[1,2]] <- changed[dims[1,2]] + 1
			if(all(changed == ncol(coefs))) break
		}

		if(all(zeros <= ncol(coefs))) return(coefs) # can add more zeros otherwise
		Upsilon <- t(t(coefs) / D_inv)
	}

	else if(Butler) { # r tests of complexity 1
		biggest <- apply(abs(coefs), 2, which.max)
		if(length(unique(biggest)) < ncol(coefs)) return(coefs) # will fail check
		for(i in 1:ncol(coefs)) coefs[biggest[i], -i] <- 0
		if(all(zeros < ncol(coefs))) return(coefs) # can add more zeros otherwise
		Upsilon <- t(t(coefs) / D_inv)
	}

	else if(quasi_Yates) { # cohyperplanarity with minimal collinearity
		FC <- coefs * (coefs %*% cormat)
		for(i in 1:ncol(coefs)) {
			diffs <- FC - FC[,i]
			coefs[apply(diffs, 2, which.max)[-i], i] <- 0
		}
		if(all(zeros < ncol(coefs))) return(coefs) # can add more zeros otherwise
		Upsilon <- t(t(coefs) / D_inv)
	}

	else if(viral) { # r zeros per column bunched in tests of complexity r - 2
		Upsilon_abs <- abs(Upsilon)
		for(i in 1:(ncol(coefs) - 1)) {
			ranks <- order(Upsilon_abs[,i])
			marks <- rep(FALSE, ncol(coefs))
			marks[i] <- TRUE
			for(j in 1:nrow(coefs)) { ## explain logic here; WTF is the logic?
				ranks2 <- order(Upsilon_abs[ranks[j],])
				if(i == ranks2[1] && !marks[ranks2[2]]) {
					Upsilon_abs[ranks[j], ranks2[2]] <- 0
					marks[ranks2[2]] <- TRUE
					Upsilon_abs[ranks[j], i] <- 0
				}
				else if(!marks[ranks2[1]]) {
					Upsilon_abs[ranks[j], ranks2[1]] <- 0
					marks[ranks2[1]] <- TRUE
					Upsilon_abs[ranks[j], i] <- 0
				}
				if(all(marks)) break
			}
		}
		coefs[!Upsilon_abs] <- 0
		if(all(zeros < ncol(coefs))) return(coefs) # can add more zeros otherwise
		Upsilon <- t(t(coefs) / D_inv)
	}
# 	else if(mve) { ## FIXME mve
# 		stop("this should not have happened, please contact maintainer")
# 		worst <- FAiR_mve(coefs) # need to write a new function for this
# 		mark  <- apply(abs(coefs[worst,]), 2, which.min)
# 		coefs[worst,mark] <- 0
# 		Upsilon <- sweep(coefs, 2, D_inv, "/", check.margin = FALSE)
# 		for(i in 1:ncol(coefs)) {
# 			worst <- FAiR_mve(coefs[,-i], zeros[i])
# 			coefs[worst,i] <- 0
# 		}
# 		return(coefs)
# 	}
	else if(communality) { # One zero per factor on a high communality outcome
		FC <- coefs * (coefs %*% cormat)
		ratios <- apply(FC, 1, FUN = FAiR_AM_GM_ratio)
		smallest <- apply(abs(Upsilon), 1, which.min)
		orderer <- order(ratios, decreasing = TRUE)
		unchanged <- rep(TRUE, ncol(coefs))
		for(i in 1:length(orderer)) {
			so <- smallest[orderer[i]]
			if(unchanged[so]) {
				coefs[orderer[i], so] <- 0
				unchanged[so] <- FALSE
			}
			if(!any(unchanged)) break
		}
		Upsilon <- t(t(coefs) / D_inv)
	}

	coefs.abs <- abs(coefs)
	for(i in 1:ncol(coefs)) { # make sure there are enough zeros per column
		ordered <- order(coefs.abs[,i])
		coefs[ordered[1:zeros[i]],i] <- 0
	}
	return(coefs)
}

## An object of parameter.* class inherits from the parameter class and is the 
## thing to be estimated by Factanal()
 
setClass("parameter", representation(  x = "ANY", free = "ANY", num_free = "integer", 
					select = "logical", invalid = "numeric",
					Domains = "array" ), package = "FAiR"

# Slots:
#                                                    
# Name:         x     free num_free   select  invalid
# Class:      ANY      ANY  integer  logical  numeric
# 
# Known Subclasses: 
# Class "parameter.cormat", directly
# Class "parameter.coef", directly
# Class "parameter.scale", directly
# Class "parameter.coef.nl", by class "parameter.coef", distance 2
# Class "parameter.coef.SEFA", by class "parameter.coef", distance 2
# Class "parameter.coef.SEFA.nl", by class "parameter.coef", distance 3

)   ## end setClass("parameter")

## validity checker for parameter
valid_parameter <-
function(object) {
	out <- TRUE
	if(length(object@free) != length(object@x)) {
		out <- "the 'free' slot must be same length as the 'x' slot"
	}
	else if(object@num_free != sum(object@free)) {
		out <- paste("the 'free' slot must sum to the 'num_free' slot")
	}
	else if(!(sum(object@select) %in% c(0, object@num_free))) {
		out <- paste("the 'select' slot must sum to zero or sum to the 'num_free'
				 slot")
	}
	return(out)
}

setValidity("parameter", valid_parameter)

## parameter.cormat is a correlation matrix among primary factors
setClass("parameter.cormat", representation("parameter", x = "matrix", free = "matrix"), 
	package = "FAiR"

# Slots:
#                                                    
# Name:         x     free num_free   select  invalid
# Class:   matrix   matrix  integer  logical  numeric
# 
# Extends: "parameter"

)   ## end setClass("parameter.cormat")

## validity checker for parameter.cormat
valid_parameter.cormat <-
function(object) {
	out <- TRUE
	if(nrow(object@x) != ncol(object@x)) {
		out <- "the 'x' slot must be a square matrix"
	}
	else if(any(diag(object@x) != 1)) {
		out <- "the 'x' slot must have 1.0 along its diagonal"
	}
	return(out)
}

setValidity("parameter.cormat", valid_parameter.cormat)

## parameter.coef is a primary pattern matrix (without nonlinear exact restrictions)
setClass("parameter.coef", representation("parameter", x = "matrix", free = "matrix",
						equalities = "list"), package = "FAiR"

# Slots:
#                                                                         
# Name:           x       free equalities   num_free     select    invalid
# Class:     matrix     matrix       list    integer    logical    numeric
# 
# Extends: "parameter"
# 
# Known Subclasses: 
# Class "parameter.coef.nl", directly
# Class "parameter.coef.SEFA", directly
# Class "parameter.coef.SEFA.nl", by class "parameter.coef.SEFA", distance 2

)   ## end setClass("parameter.coef")

## validity checker for parameter.coef
valid_parameter.coef <-
function(object) {
	out <- TRUE
	if(l <- length(object@equalities)) for(i in 1:l) {
		if(!is(object@equalities[[i]], "equality_restriction")) {
			out <- paste("each element of the 'equalities' slot must be",
					"an object of class 'equality_restriction'")
		}
	}
# 	else if(!is(object, "parameter.coef.SEFA")) {
# 		stop("FIXME: check Howe conditions")
# 	}
	return(out)
}

setValidity("parameter.coef", valid_parameter.coef)

## primary pattern matrix with nonlinear exact restrictions
setClass("parameter.coef.nl", representation("parameter.coef",
						nonlinearities = "function"), 
	prototype(nonlinearities = function(x, ...) {
		stop("you must replace the body of the function in the 'nonlinearities'",
			" slot of the object of class 'parameter.coef.nl'")
	}), package = "FAiR"

# Slots:
#                                                                   
# Name:  nonlinearities              x           free     equalities
# Class:       function         matrix         matrix           list
#                                                    
# Name:        num_free         select        invalid
# Class:        integer        logical        numeric
# 
# Extends: 
# Class "parameter.coef", directly
# Class "parameter", by class "parameter.coef", distance 2

)   ## end setClass("parameter.coef.nl")

## validity checker
valid_parameter.coef.nl <-
function(object) {
	out <- TRUE
	if(!("x" %in% names(formals(object@nonlinearities)))) {
		out <- paste("the function in the 'nonlinearities' slot must have an",
				"argument called 'x'")
	}
	else if(!any(is.infinite(object@x))) {
		out <- paste("when imposing nonlinear exact restrictions, at least one",
				"cell of the coefficient matrix must be Inf and this",
				"cell must be changed by the function body")
	}
	return(out)
}

setValidity("parameter.coef.nl", valid_parameter.coef.nl)

## primary pattern matrix with a mapping rule (SEFA) but no nonlinear exact restrictions
setClass("parameter.coef.SEFA", representation("parameter.coef", 
						rankcheck = "character",
						mapping_rule = "function",
						squashed = "matrix"),
	prototype(rankcheck = "reiersol", mapping_rule = mapping_rule), 
	package = "FAiR"

# Slots:
#                                                                        
# Name:  mapping_rule            x         free   equalities     num_free
# Class:     function       matrix       matrix         list      integer
#                                 
# Name:        select      invalid
# Class:      logical      numeric
# 
# Extends: 
# Class "parameter.coef", directly
# Class "parameter", by class "parameter.coef", distance 2
# 
# Known Subclasses: "parameter.coef.SEFA.nl"

)   ## end setClass("parameter.coef.SEFA")

## validity checker for parameter.coef.SEFA
valid_parameter.coef.SEFA <-
function(object) {
	out <- TRUE
	args <- formals(object@mapping_rule)
	if(!all(c("coefs", "cormat", "zeros") %in% names(args))) {
		out <- paste("the function in the 'mapping rule' slot must have the",
			"following arguments: 'coefs', 'cormat', 'zeros'")
	}
	if(ncol(object@x) <= 1) {
		out <- "SEFA is only sensible with multiple factors"
	}
	else if(object@rankcheck == "reiersol") {
		if(is.numeric(args$zeros) && any(args$zeros < ncol(object@x))) {
			out <- paste("The Reiersol (1950) theorem requires at least",
				ncol(object@x), "zeros per factor in the coefficient",
				"matrix. Your mapping rule has too few.")
		}
	}
	else if(object@rankcheck == "howe") {
		if(is.numeric(args$zeros) && any(args$zeros < ncol(object@x) - 1)) {	
			out <- paste("The Howe (1955) theorem requires at least",
				ncol(object@x) - 1, "zeros per factor in the coefficient",
				"matrix. Your mapping rule has too few.")
		}
	}
	else out <- paste("'rankcheck' must be either 'reiersol' or 'howe'")

	if(!isTRUE(out)) return(out)
	else if(!identical(names(args), names(formals(mapping_rule)))) return(out)

	if(isTRUE(args$Butler) && ncol(object@x) == 2) {
		out <- paste("The Butler mapping rule is inapplicable when the number of",
				"factors is 2")
	}
	else if(isTRUE(args$viral) && ncol(object@x) == 2) {
		out <- paste("The viral mapping rule is inapplicable when the number of",
				"factors is 2")
	}
	else if(isTRUE(args$weak_Thurstone) && ncol(object@x)^2 > nrow(object@x)) {
		out <- paste("There must be at least", ncol(object@x)^2, "outcomes to",
				"use the Thurstone mapping rule with", ncol(object@x),
				"factors")
	}
	return(out)
}

setValidity("parameter.coef.SEFA", valid_parameter.coef.SEFA)

## primary pattern matrix with a mapping rule (SEFA) and nonlinear exact restrictions
setClass("parameter.coef.SEFA.nl", representation("parameter.coef.SEFA", 
						nonlinearities = "function"), 
	prototype(nonlinearities = function(x, ...) {
		stop("you must replace the body of the function in the 'nonlinearities'",
			" slot of the object of class 'parameter.coef.SEFA.nl'")
	}, rankcheck = "reiersol", mapping_rule = mapping_rule), package = "FAiR"


# Slots:
#                                                                   
# Name:  nonlinearities   mapping_rule              x           free
# Class:       function       function         matrix         matrix
#                                                                   
# Name:      equalities       num_free         select        invalid
# Class:           list        integer        logical        numeric
# 
# Extends: 
# Class "parameter.coef.SEFA", directly
# Class "parameter.coef", by class "parameter.coef.SEFA", distance 2
# Class "parameter", by class "parameter.coef.SEFA", distance 3

)   ## end setClass("parameter.coef.SEFA.nl")

## validity checker for parameter.coef.SEFA.nl
valid_parameter.coef.SEFA.nl <-
function(object) {
	out <- TRUE
	if(!("x" %in% names(formals(object@nonlinearities)))) {
		out <- paste("the function in the 'nonlinearities' slot must have an",
				"argument called 'x'")
	}
	else if(!any(is.infinite(object@x))) {
		out <- paste("when imposing nonlinear exact restrictions, at least one",
				"cell of the coefficient matrix must be Inf and this",
				"cell must be changed by the function body")
	}
	return(out)
}

setValidity("parameter.coef.SEFA.nl", valid_parameter.coef.SEFA.nl)

## diagonal scale matrix, with estimated manifest standard deviations along the diagonal
setClass("parameter.scale", representation("parameter", x = "numeric", free = "logical"), 
	package = "FAiR"

# Slots:
#                                                    
# Name:         x     free num_free   select  invalid
# Class:  numeric  logical  integer  logical  numeric
# 
# Extends: "parameter"

)   ## end setClass("parameter.scale")

## An object of class restrictions holds a lot of necessary information when using
## Factanal(). Several restrictions.* classes inherit from restrictions and cause
## Factanal() to estimate different factor analysis models.
setClass("restrictions", representation(factors   = "integer",   # number of factors
					nvars     = "integer",   # see genoud
					dof       = "integer",   # degrees of freedom
					Domains   = "matrix",    # bounds on parameters
					model     = "character", # what model?
				   discrepancy    = "character", # discrepancy function
					     free = "logical"),  # which free
	package = "FAiR"

# Slots:
#                                                                               
# Name:      factors       nvars         dof     Domains       model discrepancy
# Class:     integer     integer     integer      matrix   character   character
#                   
# Name:         free
# Class:     logical
# 
# Known Subclasses: 
# Class "restrictions.independent", directly
# Class "restrictions.factanal", directly
# Class "restrictions.orthonormal", directly
# Class "restrictions.1storder", directly
# Class "restrictions.1storder.EFA", by class "restrictions.1storder", distance 2
# Class "restrictions.general", by class "restrictions.1storder", distance 2
# Class "restrictions.2ndorder", by class "restrictions.general", distance 3

)   ## end setClass("restrictions")

## validity checker for restrictions
valid_restrictions <-
function(object) {
	out <- TRUE
	if(length(object@factors) != 2) {
		out <- paste("The length of the 'factors' slot must be exactly two.",
			"You probably intended zero second-order factors")
	}
	else if(any(object@factors < 0)) {
		out <- "both elements of the 'factors' slot must be non-negative"
	}
	else if(object@nvars != nrow(object@Domains)) {
		out <- paste("the 'nvars' slot must equal the number of rows in the",
				"'Domains' slot")
	}
	else if(object@dof < 0) {
		out <- paste("the 'dof' slot is negative, which is either a mistake or",
				"indicates that the model is not estimatable")
	}
	else if(object@nvars != sum(object@free)) {
		out <- paste("the 'free' slot must sum to the 'nvars' slot")
	}
	return(out)
}

setValidity("restrictions", valid_restrictions)

## null model with no factors
setClass("restrictions.independent", representation("restrictions",
							Omega = "parameter.scale", 
							criteria = "list"), 
	prototype(factors = c(0L, 0L), model = "CFA"), package = "FAiR"


# Slots:
#                                                                       
# Name:            scale        criteria         factors           nvars
# Class: parameter.scale            list         integer         integer
#                                                                       
# Name:              dof         Domains           model     discrepancy
# Class:         integer          matrix       character       character
#                       
# Name:             free
# Class:         logical
# 
# Extends: "restrictions"

)   ## end setClass("restrictions.independent")

## validity checker
valid_restrictions.independent <-
function(object) {
	out <- TRUE
	if(length(object@criteria) != 1) {
		out <- paste("the 'criteria' slot must contain exactly one",
				"(discrepancy) function")
	}
	else if(!is.function(object@criteria[[1]])) {
		out <- paste("the 'criteria' slot must contain exactly one",
				"(discrepancy) function")
	}
	else if(!all(object@factors == 0)) {
		out <- "the independence model must have zero factors"
	}
	else if(length(object@Omega@select) != object@nvars) {
		out <- paste("the length of 'select' slot in the 'Omega' slot must be",
				"equal to the number of parameters to estimate")
	}
	return(out)
}

setValidity("restrictions.independent", valid_restrictions.independent)

## EFA model with the Lawley & Maxwell restrictions
setClass("restrictions.factanal", representation("restrictions"), package = "FAiR"

# Slots:
#                                                                               
# Name:         fast     factors       nvars         dof     Domains       model
# Class:     logical     integer     integer     integer      matrix   character
#                               
# Name:  discrepancy        free
# Class:   character     logical
# 
# Extends: "restrictions"

)   ## end setClass("restrictions.factanal")

## EFA model with the LISREL restrictions (upper triangle of beta is zero)
setClass("restrictions.orthonormal", representation("restrictions",
							beta  = "parameter.coef",
							Omega = "parameter.scale",
						     criteria = "list"),
	package = "FAiR"

# Slots:
#                                                                       
# Name:             beta           scale        criteria         factors
# Class:  parameter.coef parameter.scale            list         integer
#                                                                       
# Name:            nvars             dof         Domains           model
# Class:         integer         integer          matrix       character
#                                       
# Name:      discrepancy            free
# Class:       character         logical
# 
# Extends: "restrictions"

)   ## end setClass("restrictions.orthonormal")

## validity checker for restrictions.orthonormal
valid_restrictions.orthonormal <-
function(object) {
	out <- TRUE
	Lambda <- coef(object)
	factors <- ncol(Lambda)
	if(length(object@criteria) != 1) {
		out <- paste("the 'criteria' slot must contain exactly one",
				"(discrepancy) function")
	}
	else if(!is.function(object@criteria[[1]])) {
		out <- paste("the 'criteria' slot must contain exactly one",
				"(discrepancy) function")
	}
	else if(any(Lambda[upper.tri(Lambda)] != 0)) {
		out <- paste("the upper triangle of the 'x' slot in the 'beta' slot",
				"must contain all zeros")
	}
	else if(object@beta@num_free != length(Lambda) - 0.5 * factors * (factors - 1)) {
		out <- paste("all the cells of the 'x' slot in the 'beta' slot must be",
				"free except those in the upper triangle")
	}
	else if(length(object@Omega@x) != nrow(Lambda)) {
		out <- paste("the length of 'Omega@x' is inconsistent with the number",
				"of rows in 'Omega@beta'")
	}
	else if(length(object@Omega@select) != object@nvars) {
		out <- paste("the length of 'select' slot in the 'Omega' slot must be",
				"equal to the number of parameters to estimate")
	}
	return(out)
}

setValidity("restrictions.orthonormal", valid_restrictions.orthonormal)

## model with correlated factors bur no second-order factors (implicitly SEFA or CFA)
setClass("restrictions.1storder", representation("restrictions",
							Phi   = "parameter.cormat",
							beta  = "parameter.coef",
							Omega = "parameter.scale",
						     criteria = "list"),
	package = "FAiR"

# Slots:
#                                                                           
# Name:               Phi             beta            scale         criteria
# Class: parameter.cormat   parameter.coef  parameter.scale             list
#                                                                           
# Name:         rankcheck          factors            nvars              dof
# Class:        character          integer          integer          integer
#                                                                           
# Name:           Domains            model      discrepancy             free
# Class:           matrix        character        character          logical
# 
# Extends: "restrictions"
# 
# Known Subclasses: 
# Class "restrictions.1storder.EFA", directly
# Class "restrictions.general", directly
# Class "restrictions.2ndorder", by class "restrictions.general", distance 2

)   ## end setClass("restrictions.1storder")

## validity checker
valid_restrictions.1storder <-
function(object) {
	out <- TRUE
	if(!all(sapply(object@criteria, is.function))) {
		out <- paste("the 'criteria' slot must be a list of functions")
	}
	else if(ncol(coef(object)) != nrow(cormat(object))) {
		out <- paste("the number of factors in the coefficient matrix is",
				"inconsistent with the number of factors in the",
				"intercorrelation matrix")
	}
	else if(length(object@Omega@x) != nrow(coef(object))) {
		out <- paste("the length of 'Omega@x' is inconsistent with the number",
				"of rows in 'Omega@beta'")
	}
	else if(length(object@Omega@select) != object@nvars) {
		out <- paste("the length of 'select' slot in the 'Omega' slot must be",
				"equal to the number of parameters to estimate")
	}
	else if(length(object@beta@select) != object@nvars) {
		out <- paste("the length of 'select' slot in the 'beta' slot must be",
				"equal to the number of parameters to estimate")
	}
	else if((l <- length(object@Phi@select)) && l != object@nvars) {
		out <- paste("the length of 'select' slot in the 'Phi' slot must be",
				"equal to the number of parameters to estimate or zero")
	}
	return(out)
}

setValidity("restrictions.1storder", valid_restrictions.1storder)

## EFA model after that factors have been transformed
setClass("restrictions.1storder.EFA", representation("restrictions.1storder",
							Lambda = "matrix",
							orthogonal = "logical",
							Tmat = "matrix",
							Tcriteria = "list"),
	package = "FAiR"

# Slots:
#                                                                           
# Name:            Lambda       orthogonal             Tmat        Tcriteria
# Class:           matrix          logical           matrix             list
#                                                                           
# Name:               Phi             beta            scale         criteria
# Class: parameter.cormat   parameter.coef  parameter.scale             list
#                                                                           
# Name:         rankcheck          factors            nvars              dof
# Class:        character          integer          integer          integer
#                                                                           
# Name:           Domains            model      discrepancy             free
# Class:           matrix        character        character          logical
# 
# Extends: 
# Class "restrictions.1storder", directly
# Class "restrictions", by class "restrictions.1storder", distance 2

)   ## end setClass("restrictions.1storder.EFA")

## validity checker
valid_restrictions.1storder.EFA <-
function(object) {
	out <- TRUE
	factors <- object@factors[1]
	if(ncol(object@Lambda) != factors) {
		out <- paste("the number of columns of the matrix in the 'Lambda' slot",
				"must be", factors)
	}
	else if(nrow(object@Lambda) != nrow(coef(object))) {
		out <- paste("the number of rows of the matrix in the 'Lambda' slot",
				"must be", nrow(coef(object)))
	}
	else if(!identical(dim(object@Tmat), c(factors, factors))) {
		out <- paste("the matrix in the 'Tmat' slot must be square and be of",
				"order", factors)
	}
	else if(!isTRUE(all.equal(rep(1, factors), colSums(object@Tmat^2),
			check.attributes = FALSE))) {
		out <- paste("the matrix in the 'Tmat' slot must have unit-length",
				"columns")
	}
	else if(!all(sapply(object@Tcriteria, is.function)) &&
		!(length(object@Tcriteria) == 1 & is.character(object@Tcriteria[[1]]))) {
		out <- paste("the 'Tcriteria' slot must be a list of functions",
				"or a list with exactly one character string")
	}
	return(out)
}

setValidity("restrictions.1storder.EFA", valid_restrictions.1storder.EFA)

## SEFA or CFA model with one second-order factor
setClass("restrictions.general", representation("restrictions.1storder",
							Delta = "parameter.coef"),
	package = "FAiR"

# Slots:
#                                                                           
# Name:             Delta              Phi             beta            scale
# Class:   parameter.coef parameter.cormat   parameter.coef  parameter.scale
#                                                                           
# Name:          criteria        rankcheck          factors            nvars
# Class:             list        character          integer          integer
#                                                                           
# Name:               dof          Domains            model      discrepancy
# Class:          integer           matrix        character        character
#                        
# Name:              free
# Class:          logical
# 
# Extends: 
# Class "restrictions.1storder", directly
# Class "restrictions", by class "restrictions.1storder", distance 2
# 
# Known Subclasses: "restrictions.2ndorder"

)   ## end setClass("restrictions.general")

## validity checker
valid_restrictions.general <-
function(object) {
	out <- TRUE
	Delta <- loadings(object, level = 2)
	onefactor <- object@factors[2] == 1
	if(onefactor & ncol(Delta) != 1) {
		out <- paste("The second-order coefficient matrix must have exactly one",
				"column in a restrictions.general object\n. If you want",
				"more second-order factors, use restrictions.2ndorder")
	}
	else if(nrow(Delta) != ncol(coef(object))) {
		out <- paste(  "The number of rows in 'object@Delta@x' must equal the",
				"number of columns in 'object@beta@x'")
	}
	else if(length(object@Delta@select) != object@nvars) {
		out <- paste("the length of 'select' slot in the 'Delta' slot must be",
				"equal to the number of parameters to estimate")
	}
	return(out)
}

setValidity("restrictions.general", valid_restrictions.general)

## SEFA or CFA model with more than one second-order factor
setClass("restrictions.2ndorder", representation("restrictions.general",
							Xi = "parameter.cormat"),
	package = "FAiR"

# Slots:
#                                                                           
# Name:                Xi            Delta              Phi             beta
# Class: parameter.cormat   parameter.coef parameter.cormat   parameter.coef
#                                                                           
# Name:             scale         criteria        rankcheck          factors
# Class:  parameter.scale             list        character          integer
#                                                                           
# Name:             nvars              dof          Domains            model
# Class:          integer          integer           matrix        character
#                                         
# Name:       discrepancy             free
# Class:        character          logical
# 
# Extends: 
# Class "restrictions.general", directly
# Class "restrictions.1storder", by class "restrictions.general", distance 2
# Class "restrictions", by class "restrictions.general", distance 3

)   ## end setClass("restrictions.2ndorder")

## validity checker
valid_restrictions.2ndorder <-
function(object) {
	out <- TRUE
	if(ncol(loadings(object, level = 2)) != nrow(cormat(object, level = 2))) {
		out <- paste("The number of columns of the second-order coefficient",
				"matrix must equal the number of rows in the",
				"second-order intercorrelation matrix")
	}
	else if(length(object@Xi@select) != object@nvars) {
		out <- paste("the length of 'select' slot in the 'Xi' slot must be",
				"equal to the number of parameters to estimate")
	}
	return(out)
}

setValidity("restrictions.2ndorder", valid_restrictions.2ndorder)

## class for equality restrictions among coefficients
setClass("equality_restriction", representation(free  = "integer", 
						fixed = "integer",
						dims  = "integer",
					     rownames = "character",
						level = "integer"),
	package = "FAiR"

# Slots:
#                                                         
# Name:       free     fixed      dims  rownames     level
# Class:   integer   integer   integer character   integer

)

## validity checker
valid_equality_restriction <-
function(object) {
	out <- TRUE
	if(length(object@free) != 1) {
		out <- "the 'free' slot must be a single positive integer"
	}
	else if(length(object@fixed) == 0) {
		out <- paste("the 'fixed' slot must be an integer vector (possibly of",
				"length one)")
	}
	else if(length(object@dims) != 2) {
		out <- paste("the 'dims' slot must be of length two, so as to indicate",
				" the rows and columns of the primary pattern matrix")
	}
	else if(length(object@rownames) && (length(object@rownames) != object@dims[1])) {
		out <- "'rownames' must be the same length as the first dimension"
	}
	else if(length(object@level) != 1) {
		out <- "the 'level' slot must be 1L or 2L"
	}
	else if(object@free <= 0) {
		out <- "the 'free' slot must be a single positive integer"
	}
	else if(object@free > prod(object@dims)) {
		out <- paste("the 'free' slot must be a single positive integer that",
				"is less than or equal to", prod(object@dims))
	}
	else if(any(object@fixed <= 0)) {
		out <- "the 'fixed' slot must contain positive integers only"
	}
	else if(any(object@fixed > prod(object@dims))) {
		out <- paste("the 'fixed' slot must only contain positive integers that",
				"are less than or equal to", prod(object@dims))
	}
	return(out)
}

setValidity("equality_restriction", valid_equality_restriction)

setOldClass("hetcor") ## from library(polycor) by John Fox

## An object of class manifest holds a lot of necessary information when using
## Factanal(). Several manifest.* classes inherit from manifest.
setClass("manifest", representation(n.obs = "integer", how = "character"), 
	package = "FAiR"

# Slots:
#                           
# Name:      n.obs       how
# Class:   integer character
# 
# Known Subclasses: 
# Class "manifest.basic", directly
# Class "manifest.basic.userW", by class "manifest.basic", distance 2
# Class "manifest.data", by class "manifest.basic", distance 2
# Class "manifest.data.ordinal", by class "manifest.data", distance 3
# Class "manifest.data.ranks", by class "manifest.data", distance 3
# Class "manifest.data.mcd", by class "manifest.data", distance 3

)   ## end setClass("manifest")

## only the covariance matrix is available and no W matrix supplied by the user
setClass("manifest.basic", representation("manifest",   # object of class manifest
					cov = "matrix", # estimated covariance matrix
					cor = "matrix", # estimated correlation matrix
					sds = "numeric", # estimated standard deviations
					center = "numeric",# estimated mean
					diag = "logical"), # diagonal included?
	package = "FAiR"

# Slots:
#                                                                             
# Name:        cov       cor       sds    center      diag     n.obs       how
# Class:    matrix    matrix   numeric   numeric   logical   integer character
# 
# Extends: "manifest"
# 
# Known Subclasses: 
# Class "manifest.basic.userW", directly
# Class "manifest.data", directly
# Class "manifest.data.ordinal", by class "manifest.data", distance 2
# Class "manifest.data.ranks", by class "manifest.data", distance 2
# Class "manifest.data.mcd", by class "manifest.data", distance 2

)   ## end setClass("manifest.basic")

valid_manifest.basic <- 
function(object) {
	out <- TRUE
	if(nrow(object@cov) != ncol(object@cov)) {
		out <- "covariance matrix not square"
	}
	else if(nrow(object@cor) != ncol(object@cor)) {
		out <- "correlation matrix not square"
	}
	else if(length(object@sds) != ncol(object@cov)) {
		out <- "length of 'sds' is not equal to the number of manifest variables"
	}
	else if((l <- length(object@center)) && (l != ncol(object@cov))) {
		out <- paste("length of 'center' is not equal to the number of manifest",
				"variables")
	}
	else if(length(object@diag) != 1) {
		out <- "'diag' should be either TRUE or FALSE"
	}
	else if(is.null(rownames(cormat(object)))) {
		out <- "'manifest@cormat' must have rownames"
	}
	ev <- eigen(object@cor, TRUE, TRUE)$values
	if(ev[length(ev)] < sqrt(.Machine$double.eps)) {
		warning("sample correlation matrix is not positive definite\n",
			"this makes maximum likelihood estimation impossible and",
			" is problematic for other estimators as well")
	}
	return(out)
}

setValidity("manifest.basic", valid_manifest.basic)

## only the covariance matrix and the W matrix are available
setClass("manifest.basic.userW", representation("manifest.basic", acov = "dMatrix"),
	package = "FAiR"

# Slots:
#                                                                        
# Name:          acov          cov          cor          sds       center
# Class: ddenseMatrix       matrix       matrix      numeric      numeric
#                                              
# Name:          diag        n.obs          how
# Class:      logical      integer    character
# 
# Extends: 
# Class "manifest.basic", directly
# Class "manifest", by class "manifest.basic", distance 2

)

## validity checker for manifest.userW
valid_manifest.userW <- 
function(object) {
	out <- TRUE
	n <- nrow(cormat(object))
	if(object@diag) n_star <- 0.5 * n * (n + 1)
	else            n_star <- 0.5 * n * (n - 1)

	if(nrow(object@acov) != ncol(object@acov)) {
		out <- "the 'acov' slot must contain a square matrix"
	}
	else if(nrow(object@acov) != n_star) {
		out <- paste("the order of the matrix in the 'acov' slot must be",
				n_star)
	}
	return(out)
}

setValidity("manifest.basic.userW", valid_manifest.userW)

## raw data are available, implicitly all numeric
setClass("manifest.data", representation("manifest.basic",# extends manifest.basic
					X = "matrix",     # data matrix
					wt = "numeric",   # numeric vector of weights
					acov = "dMatrix"), # covmat of covmat
	package = "FAiR"

# Slots:
#                                                                        
# Name:             X           wt         acov          cov          cor
# Class:       matrix      numeric ddenseMatrix       matrix       matrix
#                                                                        
# Name:           sds       center         diag        n.obs          how
# Class:      numeric      numeric      logical      integer    character
# 
# Extends: 
# Class "manifest.basic", directly
# Class "manifest", by class "manifest.basic", distance 2
# 
# Known Subclasses: "manifest.data.ordinal", "manifest.data.ranks", "manifest.data.mcd"

)   ## end setClass("manifest.data")

## validity checker for manifest.data
valid_manifest.data <- 
function(object) {
	out <- TRUE
	n <- nrow(cormat(object))
	if(object@diag) n_star <- 0.5 * n * (n + 1)
	else            n_star <- 0.5 * n * (n - 1)

	if(nrow(object@acov) != ncol(object@acov)) {
		out <- "the 'acov' slot must contain a square matrix"
	}
	else if(nrow(object@acov) > 0 & nrow(object@acov) != n_star) {
		out <- paste("the order of the matrix in the 'acov' slot must be",
				n_star, "or (worse) zero")
	}
	else if(length(object@wt) != nrow(object@X)) {
		out <- paste("length of 'wt' must be the same as the number of rows in",
				"'X'")
	}
	return(out)
}

setValidity("manifest.data", valid_manifest.data)

## raw data are available but some are ordinal
setClass("manifest.data.ordinal", representation("manifest.data"), package = "FAiR"

# Slots:
#                                                                        
# Name:             D            X           wt         acov          cov
# Class:       matrix       matrix      numeric ddenseMatrix       matrix
#                                                                        
# Name:           cor          sds       center         diag        n.obs
# Class:       matrix      numeric      numeric      logical      integer
#                    
# Name:           how
# Class:    character
# 
# Extends: 
# Class "manifest.data", directly
# Class "manifest.basic", by class "manifest.data", distance 2
# Class "manifest", by class "manifest.data", distance 3

)

## raw data are available but were converted to ranks
setClass("manifest.data.ranks", representation("manifest.data"), 
	package = "FAiR"

# Slots:
#                                                                        
# Name:             X           wt         acov          cov          cor
# Class:       matrix      numeric ddenseMatrix       matrix       matrix
#                                                                        
# Name:           sds       center         diag        n.obs          how
# Class:      numeric      numeric      logical      integer    character
# 
# Extends: 
# Class "manifest.data", directly
# Class "manifest.basic", by class "manifest.data", distance 2
# Class "manifest", by class "manifest.data", distance 3

)

## raw data are available and the MCD estimator was used to calculate their covariance
setClass("manifest.data.mcd", representation("manifest.data", # extends manifest.data
					     mcd = "CovMcd"), # from library(rrcov)
	package = "FAiR"

# Slots:
#                                                                        
# Name:           mcd            X           wt         acov          cov
# Class:       CovMcd       matrix      numeric ddenseMatrix       matrix
#                                                                        
# Name:           cor          sds       center         diag        n.obs
# Class:       matrix      numeric      numeric      logical      integer
#                    
# Name:           how
# Class:    character
# 
# Extends: 
# Class "manifest.data", directly
# Class "manifest.basic", by class "manifest.data", distance 2
# Class "manifest", by class "manifest.data", distance 3

)   ## end setClass("manifest.data.mcd")

## FA is the basic class for objects produced by Factanal() and Rotate(). It holds
## estimates and information about the completed estimation process. There are several
## inherited classes.
setClass("FA", representation(  loadings      = "array",     # variables x factors x 5
                                correlations  = "array",     # factors   x factors x 3
                                uniquenesses  = "numeric",   # vector of length(variables)
					scale = "numeric",   # vector of length(variables)
                                restrictions  = "restrictions", # see restrictions-class
                                Jacobian      = "matrix",    # C wrt parameters
					vcov  = "matrix",    # parameters x parameters
				scores        = "matrix",    # observations x factors
                                manifest      = "manifest",  # see manifest-class
                                optimization  = "list",      # list of lists
					call  = "language",  # call to function
					seeds = "matrix"     # PRNG seeds for genoud()
                             ),       package = "FAiR"

# Slots:
#                                                                        
# Name:      loadings correlations uniquenesses        scale restrictions
# Class:        array        array      numeric      numeric restrictions
#                                                                        
# Name:      Jacobian         vcov       scores     manifest optimization
# Class:       matrix       matrix       matrix     manifest         list
#                                 
# Name:          call        seeds
# Class:     language       matrix
# 
# Known Subclasses: 
# Class "FA.EFA", directly
# Class "FA.general", directly
# Class "FA.2ndorder", by class "FA.general", distance 2

)   ## end setClass("FA")

## validity checker for FA
valid_FA <- 
function(object) {
	out <- TRUE
	n <- nrow(cormat(object@manifest))
	factors <- object@restrictions@factors[1]
	if(is(object@restrictions, "restrictions.independent")) return(out)

	if(dim(object@loadings)[1] != n) {
		out <- paste("the number of rows in the 'loadings' slot must equal the",
				"number of outcome variables,", n)
	}
	else if(dim(object@loadings)[2] != factors) {
		out <- paste("the number of columns in the 'loadings' slot must equal",
				"the number of factors,", factors)
	}
	else if(dim(object@loadings)[3] != 5) {
		out <- paste("the array in the 'loadings' slot must stack 5 matrices:",
				"primary pattern, reference structure, primary structure",
				"reference pattern, and factor contribution")
	}
	else if(!all(rownames(object@loadings) == rownames(cormat(object@manifest)))) {
		out <- paste("the row names of the array in the 'loadings' slot must",
				"match the names of the outcome variables")
	}
	else if(!all(dimnames(object@loadings)[[3]] == c("PP", "RS", "PS", "RP", "FC"))) {
		out <- paste("the names of the third dimension of the array in the",
				"'loadings' slot must be: 'PP', 'RS', 'PS', 'RP', 'FC'")
	}
	else if(dim(object@correlations)[1] != factors) {
		out <- paste("the number of rows in the 'correlations' slot must equal",
				"the number of factors,", factors)
	}
	else if(dim(object@correlations)[2] != factors) {
		out <- paste("the number of columns in the 'correlations' slot must",
				"equal the number of factors,", factors)
	}
	else if(dim(object@correlations)[3] != 3) {
		out <- paste("the array in the 'correlations' slot must stack 3",
				"intercorrelation matrices: primary factors, reference",
				"factors, and (diagonal) primary-reference")
	}
	else if(!all(dimnames(object@correlations)[[3]] == c("PF", "RF", "PR"))) {
		out <- paste("the names of the third dimension of the array in the",
				"'correlations' slot must be: 'PF', 'RF', 'PR'")
	}
	else if(factors > 1 && !isTRUE(all.equal(diag(object@correlations[,,1]), 
				rep(1,factors), check.attributes = FALSE))) {
		out <- paste("the diagonal elements of the first matrix in the",
				"'correlations' slot must all be 1.0")
	}
	else if(factors > 1 && !isTRUE(all.equal(diag(object@correlations[,,2]), 
				rep(1, factors), check.attributes = FALSE))) {
		out <- paste("the diagonal elements of the second matrix in the",
				"'correlations' slot must all be 1.0")
	}
	else if(length(object@uniquenesses) != n) {
		out <- paste("the length of the 'uniquenesses' slot must equal the",
				"number of manifest variables,", n)
	}
	else if(any(object@uniquenesses < 0)) {
		out <- paste("uniquenesses must be non-negative")
	}
	else if(length(object@scale) != n) {
		out <- paste("the length of the 'scale' slot must equal the",
				"number of manifest variables,", n)
	}
	else if(any(object@scale < 0)) {
		out <- paste("scales must be non-negative")
	}

	if(!isTRUE(out)) return(out)
	if(factors == 1) return(out)

	PF <- cormat(object, matrix = "PF")
	RF <- cormat(object, matrix = "RF")
	D <- cormat(object, matrix = "PR")
	if(!isTRUE(all.equal(1/sqrt(diag(solve(RF))), diag(D)))) {
		out <- paste("the correlation matrix among reference factors is",
			"inconsistent with the correlation matrix between primary and",
			"reference factors")
		return(out)
	}
	PP <- loadings(object, matrix = "PP")
	RS <- loadings(object, matrix = "RS")
	if(!isTRUE(all.equal(PP %*% D, RS[,1:ncol(RS)]))) {
		out <- paste("the primary pattern matrix and reference structure matrix",
				"are inconsistent")
		return(out)
	}
	PS <- loadings(object, matrix = "PS")
	RP <- loadings(object, matrix = "RP")
	if(!isTRUE(all.equal(PS[,1:ncol(PS)], RP %*% D))) {
		out <- paste("the primary structure matrix and reference pattern",
				"matrix are inconsistent")
		return(out)
	}
	FC <- loadings(object, matrix = "FC")
	if(!isTRUE(all.equal(FC, RS * RP))) {
		out <- paste("the factor contribution matrix is inconsistent with the",
				"reference axis solution")
		return(out)
	}
	if(!isTRUE(all.equal(FC, PS * PP))) {
		out <- paste("the factor contribution matrix is inconsistent with the",
				"primary axis solution")
		return(out)
	}
	return(out)
}

setValidity("FA", valid_FA)

## tailored for EFA results
setClass("FA.EFA", representation("FA",
				rotated = "logical",   # TRUE or FALSE
				Lambda = "matrix",     # variables x factors
				trans_mats = "array"), # factors x factors x 3
	package = "FAiR"

# Slots:
#                                                                        
# Name:       rotated       Lambda   trans_mats     loadings correlations
# Class:      logical       matrix        array        array        array
#                                                                        
# Name:  uniquenesses        scale restrictions     Jacobian         vcov
# Class:      numeric      numeric restrictions       matrix       matrix
#                                                                        
# Name:        scores     manifest optimization         call        seeds
# Class:       matrix     manifest         list     language       matrix
# 
# Extends: "FA"

)

## validity checker for FA.EFA
valid_FA.EFA <- 
function(object) {
	out <- TRUE
	n <- nrow(cormat(object@manifest))
	factors <- object@restrictions@factors[1]

	if(dim(object@Lambda)[1] != n) {
		out <- paste("the number of rows in the 'Lambda' slot must equal the",
				"number of outcome variables,", n)
	}
	else if(dim(object@loadings)[2] != factors) {
		out <- paste("the number of columns in the 'Lambda' slot must equal",
				"the number of factors,", factors)
	}
	else if(dim(object@trans_mats)[1] != factors) {
		out <- paste("the number of rows in the 'trans_mats' slot must equal",
				"the number of factors,", factors)
	}
	else if(dim(object@trans_mats)[2] != factors) {
		out <- paste("the number of columns in the 'trans_mats' slot must equal",
				"the number of factors,", factors)
	}
	else if(dim(object@trans_mats)[3] != 3) {
		out <- paste("the array in the 'trans_mats' slot must stack 3",
				"transformation matrices: that for the primary factors",
				"that for the reference factors and the transformation",
				"matrix found during the optimization phase")
	}
	else if(!all(dimnames(object@trans_mats)[[3]] == c("primary", "reference", "T"))){
		out <- paste("the names of the third dimension of the array in the",
				"'trans_mats' slot must be: 'primary', 'reference', 'T'")
	}
	else if(factors == 1) return(out)
	else if(!isTRUE(all.equal(primary <- object@trans_mats[,,"primary"], 
			  t(solve(object@trans_mats[,,"T"])), check.attributes = FALSE))) {
		out <- paste("transformation matrix for the primary factors must equal",
				"the transpose of the inverse of the T matrix")
	}
	else if(!isTRUE(all.equal(object@trans_mats[,,"reference"],
			sweep(primary, 2, sqrt(colSums(primary^2)), "/"),
				check.attributes = FALSE))) {
		out <- paste("transformation matrix of the reference factors must be",
				"column-wise proportion to that for the primary factors")
	}
	return(out)
}

setValidity("FA.EFA", valid_FA.EFA)

## tailored for results with one second-order factor
setClass("FA.general", representation("FA",
					restrictions = "restrictions.general",
					loadings_2nd = "array",
					uniquenesses_2nd = "numeric"), 
	package = "FAiR"

# Slots:
#                                                                      
# Name:          restrictions         loadings_2nd     uniquenesses_2nd
# Class: restrictions.general                array              numeric
#                                                                      
# Name:              loadings         correlations         uniquenesses
# Class:                array                array              numeric
#                                                                      
# Name:                 scale             Jacobian                 vcov
# Class:              numeric               matrix               matrix
#                                                                      
# Name:                scores             manifest         optimization
# Class:               matrix             manifest                 list
#                                                 
# Name:                  call                seeds
# Class:             language               matrix
# 
# Extends: "FA"
# 
# Known Subclasses: "FA.2ndorder"

)   ## end setClass("FA.general")

## validity checker for FA.general
valid_FA.general <- 
function(object) {
	out <- TRUE
	factors <- object@restrictions@factors

	if(dim(object@loadings_2nd)[1] != factors) {
		out <- paste("the number of rows in the 'loadings_2nd' slot must equal",
				"the number of first-order factors,", factors[1])
	}
	else if(dim(object@loadings_2nd)[2] != factors[2]) {
		out <- paste("the number of columns in the 'loadings_2nd' slot must",
				"equal the number of second-order factors,", factors[2])
	}
	else if(dim(object@loadings_2nd)[3] != 5) {
		out <- paste("the array in the 'loadings_2nd' slot must stack 5",
				"matrices: primary pattern, reference structure, primary",
				"structure, reference pattern, and factor contribution")
	}
	else if(!all(dimnames(object@loadings_2nd)[[3]] == 
			c("PP", "RS", "PS", "RP", "FC"))) {
		out <- paste("the names of the third dimension of the array in the",
			"'loadings_2nd' slot must be: 'PP', 'RS', 'PS', 'RP', 'FC'")
	}
	else if(length(object@uniquenesses_2nd) != factors[1]) {
		out <- paste("the length of the 'uniquenesses_2nd' slot must equal",
				"the number of first-order factors", factors[1])
	}
	else if(any(object@uniquenesses_2nd < 0)) {
		out <- "uniquenesses at level 2 must be non-negative"
	}
	return(out)
}

setValidity("FA.general", valid_FA.general)

## tailored for results with multiple second-order factors
setClass("FA.2ndorder", representation("FA.general",
					restrictions = "restrictions.2ndorder",
					correlations_2nd = "array"), 
	package = "FAiR"

# Slots:
#                                                                         
# Name:           restrictions      correlations_2nd          loadings_2nd
# Class: restrictions.2ndorder                 array                 array
#                                                                         
# Name:       uniquenesses_2nd              loadings          correlations
# Class:               numeric                 array                 array
#                                                                         
# Name:           uniquenesses                 scale              Jacobian
# Class:               numeric               numeric                matrix
#                                                                         
# Name:                   vcov                scores              manifest
# Class:                matrix                matrix              manifest
#                                                                         
# Name:           optimization                  call                 seeds
# Class:                  list              language                matrix
# 
# Extends: 
# Class "FA.general", directly
# Class "FA", by class "FA.general", distance 2

)   ## end setClass("FA.2ndorder")

## validity checker for FA.2ndorder
valid_FA.2ndorder <-
function(object) {
	out <- TRUE
	factors <- object@restrictions@factors[2]
	if(dim(object@correlations_2nd)[1] != factors) {
		out <- paste("the number of rows in the 'correlations_2nd' slot must",
				"equal the number of second-order factors,", factors)
	}
	else if(dim(object@correlations_2nd)[2] != factors) {
		out <- paste("the number of columns in the 'correlations_2nd' slot must",
				"equal the number of second-order factors,", factors)
	}
	else if(dim(object@correlations_2nd)[3] != 3) {
		out <- paste("the array in the 'correlations' slot must stack 3",
				"intercorrelation matrices at level 2: primary factors,",
				"reference factors, and (diagonal) primary-reference")
	}
	else if(!all(dimnames(object@correlations_2nd)[[3]] == c("PF", "RF", "PR"))) {
		out <- paste("the names of the third dimension of the array in the",
				"'correlations_2nd' slot must be: 'PF', 'RF', 'PR'")
	}
	else if(!isTRUE(all.equal(diag(object@correlations_2nd[,,1]), rep(1, factors),
				check.attributes = FALSE))) {
		out <- paste("the diagonal elements of the first matrix in the",
				"'correlations_2nd' slot must all be 1.0")
	}
	else if(!isTRUE(all.equal(diag(object@correlations_2nd[,,2]), rep(1, factors),
				check.attributes = FALSE))) {
		out <- paste("the diagonal elements of the second matrix in the",
				"'correlations_2nd' slot must all be 1.0")
	}

	if(!isTRUE(out)) return(out)

	PP <- loadings(object, matrix = "PP", level = 2)
	PS <- loadings(object, matrix = "PS", level = 2)
	RP <- loadings(object, matrix = "RP", level = 2)
	RS <- loadings(object, matrix = "RS", level = 2)
	FC <- loadings(object, matrix = "FC", level = 2)
	if(!isTRUE(all.equal(FC, RP * RS))) {
		out <- paste("factor contributions at level 2 are inconsistent with",
				"the reference axis solution")
	}
	else if(!isTRUE(all.equal(FC, PP * PS))) {
		out <- paste("factor contributions at level 2 are inconsistent with",
				"the primary axis solution")
	}
	return(out)
}

setValidity("FA.2ndorder", valid_FA.2ndorder)

## A class for the output produced by summary(FAobject)
setClass("summary.FA", representation(    restrictions = "restrictions",
						 draws = "list",
						 order = "integer",
					      polarity = "integer",
					    conf.level = "numeric",
					    orthogonal = "logical",
					  standardized = "logical",
						  call = "call"), package = "FAiR"

# Slots:
#                                                                        
# Name:  restrictions        draws        order     polarity   conf.level
# Class: restrictions         list      integer      integer      numeric
#                                              
# Name:    orthogonal standardized         call
# Class:      logical      logical         call

)   ## end setClass("summary.FA")
