


#' Interpret the results of the Oxford Immunotec TSPOT.TB assay for latent tuberculosis infection.
#'
#'
#' Given vectors of nil, TB antigen (panels A and B), and mitogen results
#' in spots, this function computes TSPOT qualitative interpretations.
#' The function uses the Oxford Immunotec North America criterion by default;
#' alternative criteria sets can be created as methods for the
#' tspots.criteria function
#'
#'
#' @include equal.lengths.r
#' @include tspot.cens.r
#' @include is.wholenumber.r
#' @include trim.output.r
#'
#' @param nil A vector of nil results (in spots)
#' @param panel.a A vector of Panel A TB antigen (ESAT-6) results (in spots)
#' @param panel.b A vector of Panel B TB antigen (CFP10) results (in spots)
#' @param mito A vector of mitogen results (in spots)
#' @param criteria The name of the desired result criteria (defaults to the Oxford Immunotec criteria for North America).
#' @param verbosity The level of verbosity ("onechar", "terse", "verbose") of the output.
#' @param ... Other arguments passed to the crtieria evaluation function chosen by the "criteria" argument.
#'
#'
#'
#' @return The function returns a vector of qualitative results.  The verbosity of results depends on the argument passed to "verbosity":
#' \item{onechar }{Returns a single character indicating the result (N for Negative, B for Borderline, P for Positive, I for Indeterminate).}
#' \item{terse }{Returns a single word indicating the result (Negative, Borderline, Positive, Indeterminate).}
#' \item{verbose }{Returns the same results as "terse", with the addition of a short comment indicating the reason for an "Indeterminate" result.}
#' 
#' Multiple criteria sets are available.  The function defaults to the
#' standard Oxford North American criteria (\code{criteria = "oxford.usa"}),
#' but other currently available options include: 
#' \item{criteria = "oxford.global"}{The Oxford global criteria, for which the
#' criterion for positivity is lowered from an 8-spot difference between the 
#' antigen and nil panels and which does not include the borderline qualitative
#' result;}
#' \item{criteria = "10spot"}{A criteria set in which the borderline result is
#' extended to include differences of 5 to 9 spots and only differences of 10 or
#' more spots indicate a positive result.}
#'
#' @details All spot counts greater than 20 are automatically censored to 20
#' for the purposes of calculating qualitative results, following Oxford's
#' interpretation instructions.
#'
#' @references Oxford Immunotec <http://www.oxfordimmunotec.com/>
#'
#' @note This function is provided purely as a convenience and is not a replacement for manual interpretation, manufacturer-provided software, or common sense.  Absolutely not for clinical use. 
#'
#' @seealso \code{\link{qft.interp}} for Quantiferon interpretation. 
#' 
#' @export
#' 
#' @examples
#' 
#' # Calculate results
#' test.tspots$result.check <- with(test.tspots, 
#'                                  tspot.interp(nil = nil, 
#'                                               panel.a = panel.a,
#'                                               panel.b = panel.b,
#'                                               mito = mito))
#' 
#' # Compare lab and calculated results
#' with(test.tspots, table(lab.result, result.check, exclude = NULL))
#' 
#' # Compare different levels of verbosity
#' test.tspots$verbose.check <- with(test.tspots, 
#'                                   tspot.interp(nil = nil, 
#'                                                panel.a = panel.a,
#'                                                panel.b = panel.b,
#'                                                mito = mito,
#'                                                verbosity = "verbose"))
#'
#' test.tspots$onechar.check <- with(test.tspots, 
#'                                   tspot.interp(nil = nil, 
#'                                                panel.a = panel.a,
#'                                                panel.b = panel.b,
#'                                                mito = mito,
#'                                                verbosity = "onechar"))
#'
#' unique(test.tspots[ , c("lab.result", "result.check", 
#'                       "verbose.check", "onechar.check")])



tspot.interp <- function(nil, panel.a, panel.b, mito,
                         criteria = "oxford.usa",
                         verbosity = "terse",
                         ...){



    # Check for equal vector lengths - throw error if not equal
    equal.lengths(nil, panel.a, panel.b, mito)


    # Check for numeric results - throw error if non-numeric
    if(any(!is.numeric(nil),
           !is.numeric(panel.a),
           !is.numeric(panel.b),
           !is.numeric(mito))){stop(
           "The vectors of TB, nil, and mitogen values must all be numeric.")}


    # Check that input values are positive - warn if negative
    if(any(nil < 0, na.rm = TRUE)){warning("One or more nil values are negative - that probably shouldn't happen!")}
    if(any(panel.a < 0, na.rm = TRUE)){warning("One or more panel.a values are negative - that probably shouldn't happen!")}
    if(any(panel.b < 0, na.rm = TRUE)){warning("One or more panel.b values are negative - that probably shouldn't happen!")}
    if(any(mito < 0, na.rm = TRUE)){warning("One or more mito values are negative - that probably shouldn't happen!")}

    # Check for non-integer results - warn if decimal
    if(any(!is.wholenumber(nil), na.rm = TRUE)){warning("One or more nil values aren't integers - that probably shouldn't happen!")}
    if(any(!is.wholenumber(panel.a), na.rm = TRUE)){warning("One or more panel.a values aren't integers - that probably shouldn't happen!")}
    if(any(!is.wholenumber(panel.b), na.rm = TRUE)){warning("One or more panel.b values aren't integers - that probably shouldn't happen!")}
    if(any(!is.wholenumber(mito), na.rm = TRUE)){warning("One or more mito values aren't integers - that probably shouldn't happen!")}


    # Censor to 20 spots
    nil.cens <- tspot.cens(nil)
    panel.a.cens <- tspot.cens(panel.a)
    panel.b.cens <- tspot.cens(panel.b)
    mito.cens <- tspot.cens(mito)


    # Set up the interpretation object
    interp.this <- data.frame(nil = nil.cens,
                              panel.a = panel.a.cens,
                              panel.b = panel.b.cens,
                              mito = mito.cens
    )

    # Set the class so that the generic function knows which method to call
    class(interp.this) <- c("data.frame", criteria)

    # Call the generic function to apply the appropriate criteria
    res <- tspot.criteria(interp.this)

    # Pare down the output as requested
    res.out <- trim.output(res, verbosity)


    return(res.out)

}




################################################################################
# Define the generic function
# The criteria argument from tspot.interp sets the class of interp.this,
# which in turn determines which criteria set is dispatched by tspot.criteria
tspot.criteria <- function(interp.this){
    UseMethod("tspot.criteria", interp.this)
}



