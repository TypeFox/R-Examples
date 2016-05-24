farules <- function(rules, statistics) {
    if (!is.list(rules)) {
        stop("'rules' must be a list of rules")
    }
    if (!is.matrix(statistics) || !is.numeric(statistics)) {
        stop("'statistics' must be a numeric matrix")
    }
    if (nrow(statistics) != length(rules)) {
        stop("The length of 'rules' must be the same as 'nrow(statistics)'")
    }
    res <- structure(list(rules=rules,
                          statistics=statistics), 
                     class=c('farules', 'list'))
    return(res)
}
