`envfit.coca` <- function(ord, env,
                          which = c("response", "predictor"),
                          choices = c(1,2),
                          scaling = FALSE,
                          w,
                          na.rm = FALSE,
                          strata = NULL,
                          permutations = 999, ...) {
    ## get weights
    if(missing(w)) {
        w <- weights(ord)
    }
    vectors <- factors <- seed <- NULL
    ## what are we plotting, response or predictor?
    which <- match.arg(which)
    ## and map to X and Y for extraction
    WHICH <- ifelse(which == "response", "Y", "X")
    ## should the scores be rescaled - only for species though
    if(is.logical(scaling))
        scaling <- ifelse(scaling, 2, 1)
    X <- scores(ord, display = "sites", choices = choices,
                scaling = scaling)
    ## then extract the response or predictor scores
    X <- lapply(X, `[[`, WHICH)[[1]]
    keep <- complete.cases(X) & complete.cases(env)
    if (any(!keep)) {
        if (!na.rm)
            stop("missing values in data: consider na.rm = TRUE")
        X <- X[keep, , drop=FALSE]
        env <- env[keep, , drop=FALSE]
        na.action <- structure(seq_along(keep)[!keep], class="omit")
    }
    if (is.data.frame(env)) {
        facts <- unlist(lapply(env, is.factor))
        if (sum(facts)) {
            Pfac <- env[, facts, drop = FALSE]
            P <- env[, !facts, drop = FALSE]
            if (length(P)) {
                if (permutations) {
                    if (!exists(".Random.seed", envir = .GlobalEnv,
                                inherits = FALSE)) {
                        runif(1)
                    }
                    seed <- get(".Random.seed", envir = .GlobalEnv,
                                inherits = FALSE)
                }
                vectors <- vectorfit(X, P, permutations, strata,
                                     choices, w = w, ...)
            }
            if (!is.null(seed)) {
                assign(".Random.seed", seed, envir = .GlobalEnv)
            }
            factors <- factorfit(X, Pfac, permutations, strata,
                                 choices, w = w, ...)
            sol <- list(vector = vectors, factors = factors)
        }
        else vectors <- vectorfit(X, env, permutations, strata,
                                  choices, w = w, ...)
    } else {
        vectors <- vectorfit(X, env, permutations, strata,
                             choices, w = w, ...)
    }
    sol <- list(vectors = vectors, factors = factors)
    if (!is.null(na.action))
        sol$na.action <- na.action
    class(sol) <- "envfit"
    sol
}
