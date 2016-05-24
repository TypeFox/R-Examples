consequents <- function(rules) {
    if (!is.list(rules)) {
        stop("'rules' must be either a list of rules or a valid farules object")
    }
    r <- rules
    if (is.farules(r)) {
        r <- r$rules
    }
    return(lapply(r, function(rule) rule[1]))
}
