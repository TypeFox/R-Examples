setup.softinfo <- function (N = 10, order=2, warn = FALSE, ...) 
{
#
# setup.softinfo()
#
if (missing (N)) stop ("N is a required argument to setup.softinfo()")

LARGE_POSITIVE_NUMBER <- 1000
LARGE_NEGATIVE_NUMBER <- -1000

pars <- list (...)

expecting <- c("upperboundsupport", "upperboundsk", "lowerboundsk",
"upperboundak0", "lowerboundak0", "upperboundakp", "lowerboundakp",
"lsc", "usc", "continuous", "continuousDiff", 
"pointwiseFisherUpper", "pointwiseFisherLower", "M", 
"monotone", "unimodal", "unimodaluppertail", "unimodallowertail", 
"lowerdensityvalue",  "lowerdensityvalueEndpt", "lowerdensityvalueSpecific",
"upperdensityvalue", "upperdensityvalueEndpt",
"KLDivergenceUpper","KLDivergenceLower", "KLDensity", "KLDensityParams",
"upperbound1moment", "upperbound2moment") 
#
# Issue a warning about missing ones, if asked
#
if (warn) {
    misses <- names (pars)[!is.element (names (pars), expecting)]
    if (length (misses) > 0)
        warning (paste ("Input list included unknown entries", paste (misses, collapse=",")))
}
#
# Set any default values
#
if (!any (names (pars) == "upperboundsk"))
    pars$upperboundsk <- rep (LARGE_POSITIVE_NUMBER, N+1)
if (!any (names (pars) == "lowerboundsk"))
    pars$lowerboundsk <- rep (LARGE_NEGATIVE_NUMBER, N+1)
if (!any (names (pars) == "upperboundak0"))
    pars$upperboundak0 <- rep (LARGE_POSITIVE_NUMBER, N)
if (!any (names (pars) == "lowerboundak0"))
    pars$lowerboundak0 <- rep (LARGE_NEGATIVE_NUMBER, N)
if (!any (names (pars) == "upperboundakp"))
    pars$upperboundakp <- rep (LARGE_POSITIVE_NUMBER, N)
if (!any (names (pars) == "lowerboundakp"))
    pars$lowerboundakp <- rep (LARGE_NEGATIVE_NUMBER, N)

if (!any (names (pars) == "M"))
    pars$M <- 5
#
# If density values are specified (and are other than the defaults),
# we will need to enforce "integrate to one" explicitly.
#

integrateToOne <- FALSE

if (any (names (pars) == "lowerdensityvalue")) {
    if (any (pars$lowerdensityvalue != rep (0, N)))
        integrateToOne <- TRUE
} else pars$lowerdensityvalue <- rep (0, N)
if (any (names (pars) == "lowerdensityvalueEndpt")) {
    if (any (pars$lowerdensityvalueEndPt != rep (0, N + 1)))
        integrateToOne <- TRUE
} else pars$lowerdensityvalueEndpt <- rep (0, N + 1)
if (any (names (pars) == "upperdensityvalue")) {
    if (any (pars$upperdensityvalue != rep (LARGE_POSITIVE_NUMBER, N)))
        integrateToOne <- TRUE
} else pars$upperdensityvalue <- rep (LARGE_POSITIVE_NUMBER, N)
if (any (names (pars) == "upperdensityvalueEndpt")) {
    if (any (pars$upperdensityvalueEndpt != rep (LARGE_POSITIVE_NUMBER, N + 1)))
        integrateToOne <- TRUE
} else pars$upperdensityvalueEndpt <- rep (LARGE_POSITIVE_NUMBER, N + 1)
#
# We also have to integrate explicity if KL or moment conditions are supplied.
#
if (any (is.element (names (pars), c("KLDensity", "upperbound1moment", "upperbound2moment"))))
    integrateToOne <- TRUE

pars$integrateToOne <- integrateToOne
#
# If unimodal is specified, turn off "continuous." You get that for free.
#
if (   (any (names (pars) == "continuous") && pars$continuous == TRUE)
    && (any (names (pars) == "unimodal")   && pars$unimodal == TRUE))
    pars$continuous <- NULL
#
# If unimodal is specified, monotonic is an error.
#
if ( (any (names (pars) == "unimodal")   && pars$unimodal == TRUE)
&& (  any (names (pars) == "monotone")   && pars$monotone != ""))
    stop ("Don't specify both unimodal and monotone.")

return (pars)


}
