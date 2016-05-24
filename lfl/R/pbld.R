.defuzzEmptyRulebase <- function(values) {
    return(defuzz(rep(1, length(values)),
                    values,
                    type='dee'))
}


pbld <- function(x,
                 rules,
                 partition, 
                 values,
                 type=c('global', 'local'),
                 parallel=FALSE) {
    if (!is.fsets(x)) {
        stop("'x' must be an instance of the 'fsets' class")
    }
    if (is.farules(rules)) {
        rules <- rules$rules
    } else if (is.vector(rules) && is.character(rules)) {
        rules <- list(rules)
    }
    if (!is.list(rules)) {
        stop("'rules' must be a list of rules")
    }
    if (!is.matrix(partition) || !is.numeric(partition)) {
        stop("'partition' must be a numeric matrix")
    }
    if (nrow(partition) <= 0 || ncol(partition) <= 0) {
        stop("'partition' must not be empty matrix")
    }
    if (is.null(colnames(partition)) || !is.character(colnames(partition))) {
        stop("'partition' must have colnames")
    }
    if (!is.vector(values) || !is.numeric(values)) {
        stop("'values' must be a numeric vector")
    }
    if (nrow(partition) != length(values)) {
        stop("The length of 'values' must be equal to the number of rows of 'partition'")
    }

    diffAttr <- setdiff(unique(unlist(antecedents(rules))), colnames(x))
    if (length(diffAttr) > 0) {
        stop(paste0("'rules' contains predicates in antecedent-part that are not present in 'colnames(x)': ",
                    paste(diffAttr, collapse=', ', sep='')))
    }

    type <- match.arg(type)

    if (length(rules) <= 0) {
        return(rep(.defuzzEmptyRulebase(values), nrow(x)))
    }


    if (type == 'global') {
        innerLoop <- function(row) {
            fired <- fire(row, rules, tnorm=goedel.tnorm, onlyAnte=TRUE)
            fired <- unlist(fired)
            maxFired <- max(fired)
            most <- which(fired == maxFired)
            selected <- .perceiveGlobal(rules[most], vars(x), specs(x))
            if (length(selected) <= 0) {
                return(.defuzzEmptyRulebase(values))
            }
            conseq <- unlist(consequents(selected))
            degrees <- aggregate(conseq, rep(maxFired, length(conseq)), partition)
            res <- defuzz(degrees, values, type='dee')
            return(res)
        }
    } else {
        innerLoop <- function(row) {
            fired <- fire(row, rules, tnorm=goedel.tnorm, onlyAnte=TRUE)
            fired <- unlist(fired)
            most <- which(fired > 0)
            selected <- .perceiveLocal(rules[most], vars(x), specs(x), fired[most])
            if (length(selected) <= 0) {
                return(.defuzzEmptyRulebase(values))
            }
            #conseq <- sapply(selected, function(rule) rule[1]) # get consequents
            conseq <- unlist(consequents(rules[most][selected]))
            degrees <- aggregate(conseq, fired[most][selected], partition)
            res <- defuzz(degrees, values, type='dee')
            return(res)
        }
    }

    i <- NULL
    if (parallel) {
        result <- foreach(i=seq_len(nrow(x))) %dopar% { innerLoop(x[i, , drop=FALSE]) }
    } else {
        result <- foreach(i=seq_len(nrow(x))) %do% { innerLoop(x[i, , drop=FALSE]) }
    }

    return(unlist(result))
}
