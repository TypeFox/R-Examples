escalc <-
function (measure, formula, ...) 
{
    if (missing(measure) || class(measure) == "formula") 
        stop("Must specify an effect size or outcome measure.")
    if (missing(formula)) 
        formula <- NULL
    UseMethod("escalc", formula)
}
