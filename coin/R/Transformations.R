### compute average scores, see Hajek, Sidak, Sen (page 131ff)
average_scores <- function(s, x) {
    for (d in unique(x))
        s[x == d] <- mean(s[x == d], na.rm = TRUE)
    return(s)
}

### identity transformation
id_trafo <- function(x) x

### rank transformation
rank_trafo <- function(x, ties.method = c("mid-ranks", "random")) {
    ties.method <- match.arg(ties.method)
    scores <- rank(x, na.last = "keep",
                   ties.method = if (ties.method == "mid-ranks") "average"
                                 else "random")
    return(scores)
}

## Klotz (1962)
klotz_trafo <- function(x, ties.method = c("mid-ranks", "average-scores")) {
    ties.method <- match.arg(ties.method)
    scores <- switch(ties.method,
        "mid-ranks" = {
            qnorm(rank_trafo(x) / (sum(!is.na(x)) + 1))^2
        },
        "average-scores" = {
            r <- rank_trafo(x, ties.method = "random")
            s <- qnorm(r / (sum(!is.na(x)) + 1))^2
            average_scores(s, x)
        }
    )
    return(scores)
}

## Mood
mood_trafo <- function(x, ties.method = c("mid-ranks", "average-scores")) {
    ties.method <- match.arg(ties.method)
    scores <- switch(ties.method,
        "mid-ranks" = {
            (rank_trafo(x) - (sum(!is.na(x)) + 1) / 2)^2
        },
        "average-scores" = {
            r <- rank_trafo(x, ties.method = "random")
            s <- (r - (sum(!is.na(x)) + 1) / 2)^2
            average_scores(s, x)
        }
    )
    return(scores)
}

### Ansari-Bradley
ansari_trafo <- function(x, ties.method = c("mid-ranks", "average-scores")) {
    ties.method <- match.arg(ties.method)
    scores <- switch(ties.method,
        "mid-ranks" = {
            r <- rank_trafo(x)
            pmin.int(r, sum(!is.na(x)) - r + 1)
        },
        "average-scores" = {
            r <- rank_trafo(x, ties.method = "random")
            s <- pmin.int(r, sum(!is.na(x)) - r + 1)
            average_scores(s, x)
        }
    )
    return(scores)
}

### Fligner
fligner_trafo <- function(x, ties.method = c("mid-ranks", "average-scores")) {
    ties.method <- match.arg(ties.method)
    scores <- switch(ties.method,
        "mid-ranks" = {
            qnorm((1 + rank_trafo(abs(x)) / (sum(!is.na(x)) + 1)) / 2)
        },
        "average-scores" = {
            s <- qnorm((1 + rank_trafo(abs(x), ties.method = "random") /
                          (sum(!is.na(x)) + 1)) / 2)
            average_scores(s, x)
        }
    )
    return(scores)
}

### Normal Scores (van der Waerden)
normal_trafo <- function(x, ties.method = c("mid-ranks", "average-scores")) {
    ties.method <- match.arg(ties.method)
    scores <- switch(ties.method,
        "mid-ranks" = {
            qnorm(rank_trafo(x) / (sum(!is.na(x)) + 1))
        },
        "average-scores" = {
            s <- qnorm(rank_trafo(x, ties.method = "random") /
                         (sum(!is.na(x)) + 1))
            average_scores(s, x)
        }
    )
    return(scores)
}

### Median Scores
median_trafo <- function(x, mid.score = c("0", "0.5", "1")) {
    ## "0.5" => symmetric median scores (Randles & Wolfe, 1979, pp. 264--266)
    x <- as.numeric(x)
    mid.score <- match.arg(mid.score)
    md <- median(x, na.rm = TRUE)
    scores <- as.numeric(x > md)
    if (mid.score != "0")
        scores[x == md] <- as.numeric(mid.score)
    return(scores)
}

### Savage scores
savage_trafo <- function(x, ties.method = c("mid-ranks", "average-scores")) {
    ties.method <- match.arg(ties.method)

    r <- function(x, t) rank(x, na.last = "keep", ties.method = t)

    scores <- switch(ties.method,
        "mid-ranks" = {
            s <- 1 / (sum(!is.na(x)) - r(x, "min") + 1)
            cumsum(s[order(x)])[r(x, "max")] - 1
        },
        "average-scores" = {
            o <- order(x)
            s <- 1 / (sum(!is.na(x)) - r(x, "first") + 1)
            average_scores(cumsum(s[o])[order(o)], x) - 1
        }
    )
    return(scores)
}

### Conover & Salsburg (1988)
consal_trafo <- function(x, ties.method = c("mid-ranks", "average-scores"),
                         a = 5) {
    ties.method <- match.arg(ties.method)

    cs <- function(a) {
        switch(ties.method,
            "mid-ranks" = {
                (rank_trafo(x) / (sum(!is.na(x)) + 1))^(a - 1)},
            "average-scores" = {
                 s <- (rank_trafo(x, ties.method = "random") /
                         (sum(!is.na(x)) + 1))^(a - 1)
                 average_scores(s, x)}
        )
    }

    scores <- if (length(a) == 1) cs(a)
              else vapply(setNames(a, paste("a =", a)), cs, as.double(x))
    return(scores)
}

## Koziol-Nemec (1979, p. 46, eq. 2.6)
koziol_trafo <- function(x, ties.method = c("mid-ranks", "average-scores"),
                         j = 1) {
    ties.method <- match.arg(ties.method)
    scores <- switch(ties.method,
        "mid-ranks" = {
            sqrt(2) * cospi(j * rank_trafo(x) / (sum(!is.na(x)) + 1))
        },
        "average-scores" = {
            s <- sqrt(2) * cospi(j * rank_trafo(x, ties.method = "random") /
                                   (sum(!is.na(x)) + 1))
            average_scores(s, x)
        }
    )
    return(scores)
}

### maximally selected (rank, chi^2, whatsoever) statistics
### ordered x
find_cutpoints <- function(x, minprob, maxprob, names) {
    qx <- quantile(x, probs = c(minprob, maxprob), na.rm = TRUE, names = FALSE,
                   type = 1)
    if (diff(qx) < .Machine$double.eps)
        return(NULL)
    ux <- sort(unique(x))
    ux <- ux[ux < max(x, na.rm = TRUE)]
    if (mean(x <= qx[2], na.rm = TRUE) <= maxprob) {
        cutpoints <- ux[ux >= qx[1] & ux <= qx[2]]
    } else {
        cutpoints <- ux[ux >= qx[1] & ux < qx[2]]
    }
    cm <- .Call("R_maxstattrafo", as.double(x), as.double(cutpoints),
                PACKAGE = "coin")
    if(names)
        dimnames(cm) <- list(1:nrow(cm), paste0("x <= ", round(cutpoints, 3)))
    cm
}

maxstat_trafo <- function(x, minprob = 0.1, maxprob = 1 - minprob) {
    cm <- find_cutpoints(x, minprob, maxprob, names = TRUE)
    cm[is.na(x)] <- NA
    cm
}

ofmaxstat_trafo <- function(x, minprob = 0.1, maxprob = 1 - minprob) {
    x <- ordered(x) # drop unused levels
    lx <- levels(x)
    x <- as.numeric(x)
    cm <- find_cutpoints(x, minprob, maxprob, names = FALSE)
    dimnames(cm) <- list(seq_len(nrow(cm)),
                         lapply(seq_len(ncol(cm)), function(i) {
                             idx <- seq_len(i)
                             paste0("{", paste0(lx[idx], collapse = ", "),
                                    "} vs. {",
                                    paste0(lx[-idx], collapse = ", "), "}")
                         }))
    cm[is.na(x)] <- NA
    cm
}

### compute index matrix of all 2^(nlevel - 1) possible splits
### code translated from package `tree'
fsplits <- function(nlevel) {

    mi <- 2^(nlevel - 1)
    index <- matrix(0, nrow = mi, ncol = nlevel)
    index[, 1] <- 1

    for (i in 0:(mi - 1)) {
        ii <- i
        for (l in 2:nlevel) {
            index[(i + 1), l] <- (ii %% 2)
            ii <- ii %/% 2
        }
    }
    storage.mode(index) <- "logical"
    index[-nrow(index),, drop = FALSE]
}

### set up transformation g(x) for all possible binary splits in an unordered x
fmaxstat_trafo <- function(x, minprob = 0.1, maxprob = 1 - minprob) {

    x <- factor(x) # drop unused levels
    sp <- fsplits(nlevels(x))
    lev <- levels(x)
    tr <- matrix(0, nrow = length(x), ncol = nrow(sp))
    cn <- vector(mode = "character", length = nrow(sp))
    for (i in 1:nrow(sp)) {
        tr[ ,i] <- x %in% lev[sp[i, ]]
        cn[i] <- paste0("{", paste0(lev[sp[i, ]], collapse = ", "), "} vs. {",
                        paste0(lev[!sp[i, ]], collapse = ", "), "}")
    }
    tr[is.na(x), ] <- NA
    dimnames(tr) <- list(1:length(x), cn)
    tr <- tr[, colMeans(tr, na.rm = TRUE) >= minprob &
               colMeans(tr, na.rm = TRUE) <= maxprob,
             drop = FALSE]
    tr
}

### weighted logrank scores; with three different methods of handling ties
logrank_trafo <-
    function(x, ties.method = c("mid-ranks", "Hothorn-Lausen", "average-scores"),
             weight = logrank_weight, ...)
{
    ## backwards compatibility
    ties.method <- match.arg(ties.method,
                             choices = c("mid-ranks", "Hothorn-Lausen",
                                         "average-scores", "logrank", "HL"),
                             several.ok = TRUE)[1]
    if (ties.method == "logrank") ties.method <- "mid-ranks"
    else if (ties.method == "HL") ties.method <- "Hothorn-Lausen"

    cc <- complete.cases(x)
    time <- x[cc, 1]
    event <- x[cc, 2]

    n <- length(time)

    if (ties.method == "average-scores") {
        noise <- runif(n, max = min(diff(sort(unique(time)))) / 2)
        time0 <- time
        time <- time - event * noise # break tied events at random
    }

    r <- rank(time, ties.method = if (ties.method != "Hothorn-Lausen") "min"
                                  else "max")
    o <- order(time, event)
    or <- r[o]
    uor <- unique(or)

    ## number at risk, number of ties and events at the ordered unique times
    n_risk <- n - uor + 1L
    n_ties <- if (ties.method != "Hothorn-Lausen") -diff(c(n_risk, 0L))
              else -diff(c(n - unique(rank(time, ties.method = "min")[o]) + 1L, 0L))
    n_event <- vapply(uor, function(i) sum(event[o][i == or]), NA_real_)

    ## index: expands ties and returns in original order
    idx <- rep.int(seq_along(n_ties), n_ties)[r] # => uor[idx] == r

    ## weights
    w <- weight(sort(unique(time)), n_risk, n_event, ...)

    ## weighted log-rank scores
    nw <- NCOL(w)
    if (nw == 1L) {
        scores <- rep.int(NA_real_, length(cc))
        scores[cc] <-
            if (ties.method != "average-scores")
                cumsum(w * n_event / n_risk)[idx] - event * w[idx]
            else # average over events only
                average_scores(cumsum(w * n_event / n_risk)[idx] - event * w[idx],
                               time0 + (1 - event) * noise)
    } else {
        scores <- matrix(NA_real_, nrow = length(cc), ncol = nw,
                         dimnames = list(NULL, colnames(w)))
        scores[cc, ] <-
            vapply(seq_len(nw), function(i) {
                if (ties.method != "average-scores")
                    cumsum(w[, i] * n_event / n_risk)[idx] - event * w[idx, i]
                else # average over events only
                    average_scores(cumsum(w[, i] * n_event / n_risk)[idx] - event * w[idx, i],
                                   time0 + (1 - event) * noise)
            }, time)
    }
    return(scores)
}

### some popular logrank weights
logrank_weight <- function(time, n.risk, n.event,
                           type = c("logrank", "Gehan-Breslow", "Tarone-Ware",
                                    "Prentice", "Prentice-Marek",
                                    "Andersen-Borgan-Gill-Keiding",
                                    "Fleming-Harrington", "Self"),
                           rho = NULL, gamma = NULL) {
    type <- match.arg(type)

    ## weight functions
    w <- function(rho, gamma) {
        switch(type,
            "logrank" = { # Mantel (1966), Peto and Peto (1972), Cox (1972)
                rep.int(1L, length(time))
            },
            "Gehan-Breslow" = { # Gehan (1965), Breslow (1970)
                n.risk
            },
            "Tarone-Ware" = { # Tarone and Ware (1977)
                n.risk^rho
            },
            "Prentice" = { # Prentice (1978), Leton and Zuluaga (2001)
                cumprod(n.risk / (n.risk + n.event)) # S(t)
            },
            "Prentice-Marek" = { # Prentice and Marek (1979)
                cumprod(1 - n.event / (n.risk + 1)) # S(t)
            },
            "Andersen-Borgan-Gill-Keiding" = { # Andersen et al (1982)
                surv <- cumprod(1 - n.event / (n.risk + 1)) # S(t)
                c(1, surv[-length(surv)]) * n.risk / (n.risk + 1) # S(t-), pred.
            },
            "Fleming-Harrington" = { # Fleming and Harrington (1991)
                surv <- cumprod(1 - n.event / n.risk) # S(t), Kaplan-Meier
                surv <- c(1, surv[-length(surv)]) # S(t-)
                surv^rho * (1 - surv)^gamma
            },
            "Self" = { # Self (1991)
                ## NOTE: this allows for arbitrary follow-up times
                v <- (time - diff(c(0, time)) / 2) / max(time[n.event > 0])
                v^rho * (1 - v)^gamma
            }
        )
    }

    ## set defaults and eliminate 'rho' and 'gamma' when redundant
    if (type == "Tarone-Ware") {
        if (is.null(rho)) rho <- 0.5
        gamma <- NULL
    } else if (type %in% c("Fleming-Harrington", "Self")) {
        if (is.null(rho)) rho <- 0
        if (is.null(gamma)) gamma <- 0
    } else rho <- gamma <- NULL

    ## find rho-gamma combinations, recycle if necessary, and re-assign
    rho_gamma <- suppressWarnings(cbind(rho, gamma)) # no warning on recycling
    if (!is.null(rho)) rho <- rho_gamma[, 1]
    if (!is.null(gamma)) gamma <- rho_gamma[, 2]

    ## weights
    wgt <- if (length(rho) < 2 && length(gamma) < 2) w(rho, gamma)
           else setColnames(vapply(seq_len(nrow(rho_gamma)),
                                   function(i) w(rho[i], gamma[i]),
                                   time),
                            ## compute names
                            paste0("rho = ", rho,
                                   if (!is.null(gamma)) ", gamma = ", gamma))
    return(wgt)
}

### factor handling
f_trafo <- function(x) {
    mf <- model.frame(~ x, na.action = na.pass, drop.unused.levels = TRUE)
    if (nlevels(mf$x) == 1)
        stop("Can't deal with factors containing only one level")
    ## construct design matrix _without_ intercept
    mm <- model.matrix(~ x - 1, data = mf)
    colnames(mm) <- levels(mf$x)
    ## the two-sample situations
    if (ncol(mm) == 2)
        mm <- mm[, -2, drop = FALSE]
    return(mm)
}

### ordered factors
of_trafo <- function(x, scores = NULL) {
    if (!is.ordered(x))
        warning(sQuote(deparse(substitute(x))), " is not an ordered factor")
    if (is.null(scores)) {
        scores <- if (nlevels(x) == 2L)
                      0:1               # must be 0:1 for exact p-values
                  else if (!is.null(s <- attr(x, "scores")))
                      s
                  else
                      seq_len(nlevels(x))
    }
    if (!is.list(scores))
        scores <- list(scores)
    if (all(lengths(scores) == nlevels(x)))
        setRownames(do.call("cbind", scores)[x, , drop = FALSE], seq_along(x))
    else
        stop(sQuote("scores"), " does not match the number of levels")
}

### transformation function
trafo <- function(data, numeric_trafo = id_trafo, factor_trafo = f_trafo,
                  ordered_trafo = of_trafo, surv_trafo = logrank_trafo,
                  var_trafo = NULL, block = NULL) {

    if (!(is.data.frame(data) || is.list(data)))
        stop(sQuote("data"), " is not a data.frame or list")

    if (!is.null(block)) {
        if (!is.factor(block) || length(block) != nrow(data))
            stop(sQuote("block"), " is not a factor with ",
                 nrow(data), " elements")

        ## need to check dimension of matrix returned by
        ## user supplied functions
        ret <- trafo(data, numeric_trafo, factor_trafo, ordered_trafo, surv_trafo)

        ## apply trafo to each block separately
        for (lev in levels(block)) {
            ret[block == lev, ] <- trafo(data[block == lev, ,drop = FALSE],
                numeric_trafo, factor_trafo, ordered_trafo, surv_trafo)
        }
        return(ret)
    }

    if (!is.null(var_trafo)) {
        if (!is.list(var_trafo)) stop(sQuote("var_trafo"), " is not a list")
        if (!all(names(var_trafo) %in% names(data)))
            stop("variable(s) ",
                 names(var_trafo)[!(names(var_trafo) %in% names(data))],
                 " not found in ", sQuote("var_trafo"))
    }

    ## compute transformations for each variable
    tr <- vector(mode = "list", length = length(data))
    names(tr) <- names(data)
    for (nm in names(data)) {
        x <- data[[nm]]
        if (nm %in% names(var_trafo)) {
            tr[[nm]] <- as.matrix(var_trafo[[nm]](x))
            next
        }
        if (class(x)[1] == "AsIs") {
            if (length(class(x)) == 1) {
                x <- as.numeric(x)
            } else {
                class(x) <- class(x)[-1]
            }
        }
        if (is.ordered(x)) {
            tr[[nm]] <- as.matrix(ordered_trafo(x))
            next
        }
        if (is.factor(x) || is.logical(x)) {
            tr[[nm]] <- as.matrix(factor_trafo(x))
            next
        }
        if (is.Surv(x)) {
            tr[[nm]] <- as.matrix(surv_trafo(x))
            next
        }
        if (is.numeric(x)) {
            tr[[nm]] <- as.matrix(numeric_trafo(x))
            next
        }
        if (is.null(tr[[nm]]))
            stop("data class ", class(x), " is not supported")
    }

    ## set up a matrix of transformations
    ## when more than one factor is in play, factor names
    ## _and_ colnames of the corresponding rows are combined by `.'
    RET <- c()
    assignvar <- c()
    cn <- c()
    for (i in 1:length(tr)) {
        if (nrow(tr[[i]]) != nrow(data))
            stop("Transformation of variable ", names(tr)[i],
                 " are not of length / nrow", nrow(data))
        RET <- cbind(RET, tr[[i]])
        if (is.null(colnames(tr[[i]]))) {
            cn <- c(cn, rep.int("", ncol(tr[[i]])))
        } else {
            cn <- c(cn, paste0(ifelse(length(tr) > 1, ".", ""), colnames(tr[[i]])))
        }
        assignvar <- c(assignvar, rep.int(i, ncol(tr[[i]])))
    }
    attr(RET, "assign") <- assignvar
    if (length(tr) > 1) {
        colnames(RET) <- paste0(rep.int(names(tr), tabulate(assignvar)), cn)
    } else {
        colnames(RET) <- cn
    }
    return(RET)
}

### multiple comparisons, cf. mcp(x = "Tukey") in multcomp
mcp_trafo <- function(...) {

    args <- list(...)
    stopifnot(length(args) == 1)

    ret <- function(data) {

        x <- data[[names(args)]]
        stopifnot(is.factor(x))
        C <- args[[1]]
        if (is.character(C)) {
            C <- contrMat(table(x), C)
        } else {
            stopifnot(is.matrix(C))
            stopifnot(ncol(C) == nlevels(x))
            if (is.null(colnames(C)))
                colnames(C) <- levels(x)
            attr(C, "type") <- "User-defined"
            class(C) <- c("contrMat", "matrix")
        }

        ret <- trafo(data, factor_trafo = function(x)
            tcrossprod(model.matrix(~ x - 1, data = model.frame(~ x, na.action = na.pass)), C))
        attr(ret, "contrast") <- C
        ret
    }
    ret
}
