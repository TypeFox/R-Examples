##
## Extract the names of the raw terms used to fit the model
##
extractVarNames <- function(formula, mf)
{
    ans <- character()
    expanded <- logical()
    for (i in seq_len(length(formula)[2])) {
        tt <- terms(stats::formula(formula, rhs = i), data = mf,
                    simplify = TRUE)
        tt <- attr(tt, "term.labels")
        ans <- c(ans, tt)
        expanded <- c(expanded, rep(i == 1, length(tt)))
    }

    structure(ans, expanded = expanded)
}
