### Wrapper around etm for easier computation of cumulative incidence
### functions

etmCIF <- function(formula, data, etype, subset, na.action, failcode = 1) {

    if (missing(data)) stop("A data frame in which to interpret the formula must be supplied")
    if (missing(etype)) stop("'etype' is missing, with no default")

    Call <- match.call()
    
    ## arg.etype <- deparse(substitute(etype))
    
    mfnames <- c('formula', 'data', 'etype', 'subset', 'na.action')
    temp <- Call[c(1, match(mfnames, names(Call), nomatch=0))]
    temp[[1]] <- as.name("model.frame")
    m <- eval.parent(temp)

    n <- nrow(m)
    y <- model.extract(m, 'response')
    if (!is.Surv(y)) stop("Response must be a survival object")

    etype <- model.extract(m, "etype")
    ## cov <- model.matrix(formula, m)
    name.strata <- attr(attr(m, "terms"), "term.labels")
    if (length(name.strata) == 0) {
        cova <- rep(1, n)
    } else {
        cova <- m[[name.strata]]
    }

    ## need to deal with etype when that's a fucking factor
    if (!is.factor(etype)) etype <- factor(etype)
    levels(etype) <- c(levels(etype), "cens")
    
    ## Creating data set for using etm
    if (attr(y, "type") == "right") {
        etype[y[, 2] == 0] <- "cens"
        entry <- rep(0, n)
        exit <- y[, 1]
    } else {
        etype[y[, 3] == 0] <- "cens"
        entry <- y[, 1]
        exit <- y[, 2]
    }
    etype <- etype[, drop = TRUE]
    from <- rep(0, n)
    to <- etype
    id <- seq_len(n)
    ## cov <- cov[, ncol(cov)]
    dat.etm <- data.frame(id = id,
                          from = from,
                          to = to,
                          entry = entry,
                          exit = exit,
                          cov = cova)

    ## Now, let's use etm
    tab.cov <- sort(unique(dat.etm$cov))
    
    state.names <- as.character(c(0, as.character(sort(unique(etype[etype != "cens"])))))
    tra <- matrix(FALSE, length(state.names), length(state.names))
    tra[1, 2:length(state.names)] <- TRUE
    
    cifs <- lapply(seq_along(tab.cov), function(i) {
        etm(dat.etm[dat.etm$cov == tab.cov[i], ], state.names, tra, "cens", 0)
    })

    X <- matrix(tab.cov, nrow = 1, dimnames = list(name.strata))
    if (ncol(X) > 1)
        names(cifs) <- paste(rownames(X), X, sep = "=")
    cifs$failcode <- failcode
    cifs$call <- Call
    cifs$X <- X
    class(cifs) <- "etmCIF"

    cifs
}
