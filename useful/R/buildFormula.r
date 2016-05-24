#' Formula Builder
#'
#' Formula Builder
#'
#' Builds a formula easily given the left and right hand sides.  Right now it only handles additive formulas and not interactions unless that is specified in the character.
#'
#' @param lhs Character vector for left side of formula
#' @param rhs Character vector for right side of formula
#' @author Jared P. Lander www.jaredlander.com
#' @aliases build.formula
#' @export build.formula
#' @importFrom stats as.formula
#' @seealso formula as.formula
#' @return A formula object
#' @examples
#'
#' build.formula("Y", "X")
#' build.formula(c("Y", "Z"), "X")
#' build.formula("Z", c("X", "Q"))
#' build.formula(c("Y", "Z"), c("X", "Q"))
#'
build.formula <- function(lhs, rhs)
{
    if(is.null(lhs) && is.null(rhs))
    {
        return(NULL)
    }
    
    # build a formula for aggregation
    theFormula <- as.formula(sprintf("%s ~ %s", paste(lhs, collapse=" + "), paste(rhs, collapse=" + ")))
    return(theFormula)
}
