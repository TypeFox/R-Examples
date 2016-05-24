anovaPE <-
function (object) 
{
    if (data.class(object) != "lm") 
        stop("The argument 'object' must be of class 'lm'.")
    name.x <- attr(terms(object), "term.labels")
    name.x.1 <- name.x[1]
    num.x <- length(name.x)
    if (num.x >= 2 && (length(grep(name.x.1, name.x)) < num.x)) {
        stop(paste("All predictor variables in the model must be functions of", 
            "a single variable; for example, x, x^2, etc."))
    }
    x <- model.frame(object)[, name.x.1]
    dcx <- data.class(x)
    if (!((dcx == "AsIs" & is.numeric(x)) || dcx == "numeric")) {
        if (num.x == 1) {
            stop("The single predictor variable must be numeric.")
        }
        else {
            stop("The single variable that all predictors are functions of must be numeric.")
        }
    }
    n.x <- length(x)
    n.unique.x <- length(unique(x))
    if (n.x == n.unique.x) {
        stop(paste("There must be at least two replicate values", 
            "for at least one value of the predictor variable."))
    }
    if (object$df.residual <= (n.x - n.unique.x)) {
        stop(paste("Not enough replicate values for predictors relative to", 
            "degrees of freedom for residual sums of squares."))
    }
    new.object <- update(object, formula = formula(paste(". ~ . + as.factor(", 
        name.x[1], ")", sep = "")), singular = TRUE)
    aov.object <- summary.aov(new.object)[[1]]
    rn <- row.names(aov.object)
    rn[1 + num.x] <- "Lack of Fit"
    rn[2 + num.x] <- "Pure Error"
    row.names(aov.object) <- format(rn)
    aov.object
}
