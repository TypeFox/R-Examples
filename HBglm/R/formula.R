##############################################################################
#  Formula Parsing Routine  
#  
#  A formula for hbglm is of type:
#
#  response ~   lower level pooled covariates   (group varying)
#             | upper level covariates  
#             | name of col with group IDs 
#                           
###############################################################################
# Argument: 
#   f - formula object or anything coercable to a Formula object
#
# Output: 
#  list object with following keys:
#    formula:   the input formula as an object of class Formula
#    var.names: name of all variables in formula (inlcuding response)
#    response:  name of response variable
#    intercept: logical vector of length 2 indicating whether Intercept   
#               is included in the first & second part of formula respectively
#    lower.cov: vector of names of lower level pooled covariates
#    upper.cov: vector of names of upper level covariates
#    pooled:    TRUE only if lower level covariates are pooled
#    grpID.col: name of column with groupID 
#
# Note: 
#    1. Presence of '-1' or '0' indicates intercept is turned of 
#    2. For NO lower level pooling, use '-1' or '0' in the 2nd part of the RHS
#    3. 2nd part can NOT be omitted from formula.
#    4. ONLY intercept in the 2nd part means that each pool is drawn from a
#       normal distribution with unknown mean and variance.
#    5. This function will NOT expand categories or interactions
###############################################################################
parseFormula <- function(f)
{
   if (is.null(f)) return(NULL)
   if (!is.Formula(f)) f <- Formula(f)
   call <- formula(f)
   parsed.fm <- list("var.names" = all.vars(f))
   parsed.fm$formula <- f 
   numLHS <- length(f)[1]
   numRHS <- length(f)[2]

   # Checking
   if (numLHS != 1 || numRHS != 3)
       stop("Invalid formula supplied: RHS must have 3 slots and LHS 1.")
   
   # Get response
   lhsTerms <- terms(f, lhs=1, rhs=0)
   response <- all.vars(attr(lhsTerms, "variables"))
   if (length(response) != 1)
       stop("Invalid formula: response (LHS) can have only 1 term.")
   parsed.fm$response <- response 

   # Lower level covariates 
   vars <- terms(f, lhs=0, rhs=1)
   x <- all.vars(attr(vars, "variables"))
   parsed.fm$lower.cov <- if (length(x) > 0) x
                          else NULL
   parsed.fm$intercept <- c(attr(vars, "intercept"), FALSE)
   
   # Upper level covariates 
   vars <- terms(f, lhs=0, rhs=2)
   x <- all.vars(attr(vars, "variables"))
   parsed.fm$upper.cov <- if (length(x) > 0) x
                          else NULL
   parsed.fm$intercept[2] <- attr(vars, "intercept")
   # Flag whether lower level is unpooled
   parsed.fm$pooled <- parsed.fm$intercept[2] || 
                       !is.null(parsed.fm$upper.cov)

   # Group ID column
   vars <- terms(f, lhs=0, rhs=3)
   x <- all.vars(attr(vars, "variables"))
   parsed.fm$grpID.col <- x[1]
   return(parsed.fm)
}
###############################################################################
## Formula examples
###############################################################################
## a and b are group-varying & completely unpooled
## print(parseFormula(y ~ a + b | 0 | county))
## print(parseFormula(y ~ a + b | -1 | county))

## Intercept in lower & upper levels
## print(parseFormula(y ~ a+b | c+d | group))

## Only intercep in the lower level
## print(parseFormula(y ~ 1 | c+d | group))

## a and b are group-varying & pooled; no upper-level explanatory variables
## Only intercept in upper level
## print(parseFormula(y ~ a+b | 1 | group))

## More complicated formula
## print(parseFormula(y ~ a * b | c + d | group))
###############################################################################
