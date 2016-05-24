#### fuzzy inference

`%is%` <- function(x, y) gset_charfun(x)(y)
## FUN <- matchfun(function(x,y) isTRUE(all.equal(x,y)))
##`%is%` <- function(x, y) cset_charfun(cset(x, matchfun = FUN))(y)

fuzzy_inference <-
function(system, values, implication = c("minimum", "product"))
{
    if (is.character(implication))
        implication <-
            switch(match.arg(implication), minimum = pmin, product = `*`)

    transform_rule <- function(expr) {
        as.call(switch(as.character(expr[[1L]]),
                       "&&" = list(as.name(".T."),
                                   Recall(expr[[2L]]), Recall(expr[[3L]])),
                       "||" = list(as.name(".S."),
                                   Recall(expr[[2L]]), Recall(expr[[3L]])),
                       "!" = list(as.name(".N."), Recall(expr[[2L]])),
                       "%is%" = list(expr[[1L]], as.call(list(as.name("$"),
                                     expr[[2L]], expr[[3L]])),
                                     values[[ as.character(expr[[2L]]) ]]),
                       "(" = list(expr[[1L]], Recall(expr[[2L]])),
                       expr))
    }

    antecedents <-
        sapply(lapply(system$rules,
                      function(i) transform_rule(i$antecedent)),
               eval, system$variables)
    antecedents[is.na(antecedents)] <- 0

    consequents <-
        lapply(system$rules,
               function(i) eval(`[[<-`(i$consequent, 1L, as.name("$")),
                                   system$variables))

    do.call(gset_union,
            Map(gset_transform_memberships,
                consequents,
                rep.int(list(implication), length(antecedents)),
                antecedents))
}

fuzzy_rule <-
function(antecedent, consequent)
{
    .structure(list(antecedent = substitute(antecedent),
                    consequent = substitute(consequent)),
               class = "fuzzy_rule")
}

print.fuzzy_rule <-
function(x, ...)
{
    writeLines(paste(format(x$antecedent), "=>", format(x$consequent)))
    invisible(x)
}

fuzzy_system <-
function(variables, rules)
    .structure(list(variables = variables, rules = rules),
               class = "fuzzy_system")


print.fuzzy_system <-
function(x, ...)
{
    writeLines(sprintf("A fuzzy system consisting of %i variables and %i rules.\n", length(x$variables), length(x$rules)))

    writeLines("Variables:\n")
    writeLines(paste(names(x$variables),
                     "(",
                     lapply(x$variables,
                            function(i) paste(names(i), collapse= ", ")),
                     ")", sep = ""))

    writeLines("\nRules:\n")
    for (i in x$rules) print(i)

    invisible(x)
}

plot.fuzzy_system <-
function(x,
         colfun = function(i) grey.colors(i, start = 0.5, end = 0), ...)
{
    vars <- as.list(x$variables)
    nams <- names(vars)

    n <- length(vars)
    nc <- ceiling(sqrt(n))
    nr <- ceiling(n / nc)
    op <- par(mfrow = c(nr, nc))
    on.exit(par(op))

    for (i in seq_len(n))
        plot(vars[[i]], main = nams[i], col = colfun(length(vars[[i]])))

    invisible(x)
}

fuzzy_partition <-
function(varnames, FUN = fuzzy_normal, universe = NULL, ...)
{
    n <- length(varnames)
    if (is.numeric(varnames)) {
        nam <- names(varnames)
        empty <- !nzchar(nam)
        if (is.null(nam) || any(is.na(empty)) || any(empty))
            stop("Position vector must be named.")
        n <- as.numeric(varnames)
        varnames <- names(varnames)
    }
    .structure(fuzzy_tuple(FUN = FUN, n = n, ...,
                           universe = universe, names = varnames),
               class = c("fuzzy_variable", "tuple"))
}

fuzzy_variable <-
function(...) {
    n <- names(l <- list(...))
    if (is.null(n) || any(!nzchar(n)))
        stop("Argument list must be named.")
    universe <- .expand()
    for (i in seq_along(l))
        if (is.function(l[[i]]))
            l[[i]] <- gset(charfun = l[[i]], universe = universe)
    .structure(as.tuple(l), class = c("fuzzy_variable", "tuple"))
}

print.fuzzy_variable <-
function(x, ...)
{
    writeLines(sprintf("A fuzzy variable with values: %s",
                       paste(names(x), collapse = ", ")))

    invisible(x)

}

plot.fuzzy_variable <-
function(x, ypos = 1, ...)
{
    plot.tuple(x, ...)
    op <- par(xpd = TRUE)
    on.exit(par(op))
    pos = sapply(x, gset_defuzzify, "meanofmax")
    text(pos, ypos, labels = names(pos), pos = 3L)

    invisible(x)
}
