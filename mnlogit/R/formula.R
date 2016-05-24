##############################################################################
#  Formula Parsing Routine  
#  
#  A formula for mnlogit is of type:
#
#  response ~   choice specific variables with generic coefficients
#             | individual specific variables with choice specific coefficients
#             | choice specific variables with individual variation and choice
#               specific coefficiients                    
#                           
###############################################################################
# Argument: 
#   f - formula object or anything coercable to a Formula object
#
# Output: 
#  Input formula object with attributes:
#    varNames: name of all variables in formula (inlcuding response)
#    response: name of response variable
#    Intercept: logical variable indicating whether an additional 
#              'intercept' parameter (which is choice sp) is estimated.
#    csvGenCoeff: vector of choice specific variables with generic coeff
#    indSpVar: vector of indivicual specific variables 
#    csvChCoeff: vector of choice specific variables with choice sp coeff
#
# Note: 
#    1. Presence of '-1' or '0' indicates intercept is turned of 
#    2. To NOT include any variable type, use a '1' as a placeholder
###############################################################################
parseFormula <- function(f)
{
   if (!is.Formula(f)) f <- Formula(f)
   call <- formula(f)
   attr(call, "varNames") <- all.vars(f)
   numLHS <- length(f)[1]
   numRHS <- length(f)[2]

   # Checking
   if (numLHS != 1 && numRHS >= 1 && numRHS <= 3)
       stop("Invalid formula supplied.")
   
   # Get response
   lhsTerms <- terms(f, lhs=1, rhs=0)
   response <- all.vars(attr(lhsTerms, "variables"))
   if (length(response) != 1)
       stop("Invalid formula: response (LHS) can have only 1 term.")

   interceptON <- TRUE
   # Choice specific with generic coeffs
   vars <- terms(f, lhs=0, rhs=1)
   x <- attr(vars, "term.labels")
   attr(call, "csvGenCoeff") <- if (length(x) > 0) x
                                else NULL
   interceptON <- (interceptON && attr(vars, "intercept"))
   
   # Individual specific vars
   if (numRHS > 1) {
       vars <- terms(f, lhs=0, rhs=2)
       x <- attr(vars, "term.labels")
       attr(call, "indSpVar") <- if (length(x) > 0) x
                                 else NULL
       interceptON <- (interceptON && attr(vars, "intercept"))
   }

   # Choice specific with choice specific coeffs
   if (numRHS > 2) {
       vars <- terms(f, lhs=0, rhs=3)
       x <- attr(vars, "term.labels")
       attr(call, "csvChCoeff") <- if (length(x) > 0) x
                                   else NULL
       interceptON <- (interceptON && attr(vars, "intercept"))
   }
    attr(call, "Intercept") <- interceptON 
    attr(call, "response")  <- response 
    return(call)
}
