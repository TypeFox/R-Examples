# A version of glmerControl, from lme4,
# with some options changed from their defaults
lme4_control <- function()
{
  lme4::glmerControl(check.nobs.vs.rankZ = "ignore",
                     check.nobs.vs.nlev = "ignore",
                     check.nlev.gtreq.5 = "ignore",
                     check.nlev.gtr.1 = "ignore",
                     check.nobs.vs.nRE = "ignore",
                     check.rankX = c("message+drop.cols","silent.drop.cols", "warn+drop.cols",
                                     "stop.deficient", "ignore"),
                     check.scaleX  = "warning",
                     check.formula.LHS = "stop",
                     check.response.not.const = "ignore")
}

# code modified from use of control in optim
find_control_with_defaults <- function(control, method)
{
  if( length(method) == 0 )
    stop("You must specify which method to use for likelihood approximation", call. = FALSE)

  conLaplace <- list()
  conAGQ <- list(nAGQ = 15)
  conSR <- list(nSL = 3)
  conIS <- list(nIS = 1000)
  con_tot <- c(conAGQ, conSR, conIS)

  con <- switch(method,
                "Laplace" = conLaplace,
                "AGQ" = conAGQ,
                "SR" = conSR,
                "IS" = conIS,
                stop(paste("The method", method, "is not recognised"), call. = FALSE))

  are_known <- names(control) %in% names(con_tot)
  names_known <- names(control)[are_known]
  names_unknown <- names(control)[!are_known]

  are_needed <- names(control[are_known]) %in% names(con)
  names_not_needed <- names_known[!are_needed]

  if ( length(names_unknown) > 0 )
    warning("unknown names in control: ", paste(names_unknown, collapse = ", "),
            call. = FALSE)

  if ( length(names_not_needed) > 0 )
    warning("For method = ", method, "parts of control were ignored: ",
            paste(names_not_needed, collapse = ", "), call. = FALSE)

  con[names(control)] <- control
  if(method == "SR") {
    con$nAGQ <- 2^(con$nSL + 1) - 1
  }

  con
}
