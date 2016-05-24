relation_impute <-
function(x, method = NULL, control = list(), ...)
{

    .methods <- c("omit",
                  "any/G",
                  "any/L", "any/L/first", "any/L/last",
                  "any/W", "any/W/first", "any/W/last",
                  "any/O", "any/O/first", "any/O/last",
                  "average/G",
                  "average/L", "average/L/first", "average/L/last",
                  "average/W", "average/W/first", "average/W/last",
                  "average/O", "average/O/first", "average/O/last")


    n <- unique(names(l <- list(...)))
    control[n] <- l[n]

    if (!is.relation(x) && !is.relation_ensemble(x))
        stop("Argument 'x' must be a relation or a relation ensemble.")
    if (is.relation_ensemble(x)) {
        if (!is.null(control$n) && !isTRUE(control$n == 1L)) {
            control$n <- NULL
            warning("'n > 1' not implemented for relation ensembles.")
        }
        return(relation_ensemble(list =
                                 lapply(x, relation_impute, method, control)))
    }

    I <- as.array(.incidence(.get_representation(x)))
    ## could use if (!relation_has_missings(x)), but for efficiency ...
    if (!any(is.na(I))) return(x)

    NAs <- .missing_objects(I)

    method <- if (is.null(method)) {
        R <- relation(incidence = I[-NAs,-NAs])
        method <- if (isTRUE(relation_is_linear_order(R)))
            "average/L"
        else if (isTRUE(relation_is_weak_order(R)))
            "average/W"
        else if (isTRUE(relation_is_partial_order(R)))
            "average/O"
        else "average/G"
    } else match.arg(method, .methods)

    do.call(sprintf(".impute_%s", gsub("/", "_", method)),
            list(.domain(x), I, NAs, control))
}

##############################################
## imputation code

## omit

.impute_omit <-
function(D, I, NAs, control)
    .make_relation_from_domain_and_incidence(D, `[<-`(I, is.na(I), 0))

## averaged relation for G

.impute_average_G <-
function(D, I, NAs, control)
    .make_relation_from_domain_and_incidence(D, `[<-`(I, is.na(I), 0.5))

## averaged relation for L{/first, /last}

.impute_average_L <-
function(D, I, NAs, control)
{
    o <- do.call(order, split(I, col(I)))
    ro <- `[<-`(seq_along(o), o, seq_along(o))
    I <- I[o, o, drop = FALSE]

    mn <- ncol(I)
    n <- length(NAs)
    m <- mn - n

    A <- I[1:m, 1:m]
    B <- matrix(rep(seq_len(m), n), m, n) / (m + 1)
    C <- matrix(rep(rev(seq_len(m)), each = n), n, m) / (m + 1)
    D <- matrix(0.5, n, n)
    diag(D) <- 1
    I <- rbind(cbind(A, B), cbind(C, D))[ro, ro]
    .make_relation_from_domain_and_incidence(D, I)
}

.impute_average_L_first <-
function(D, I, NAs, control)
{
    I[ ,NAs] <- 1
    I[NAs, ] <- 0
    I[NAs, NAs] <- 0.5
    diag(I) <- 1
    .make_relation_from_domain_and_incidence(D, I)
}

.impute_average_L_last <-
function(D, I, NAs, control)
{
    I[ ,NAs] <- 0
    I[NAs, ] <- 1
    I[NAs, NAs] <- 0.5
    diag(I) <- 1
    .make_relation_from_domain_and_incidence(D, I)
}

## averaged relation for W{/first, /last}

.impute_average_W <-
function(D, I, NAs, control)
    .make_relation_from_domain_and_incidence(D,
                                             .impute_average_WO(I, length(NAs)))

.impute_average_W_first <-
function(D, I, NAs, control)
    .impute_any_W_first(D, I, NAs, control)

.impute_average_W_last <-
function(D, I, NAs, control)
    .impute_any_W_last(D, I, NAs, control)

## averaged relation for O{/first, /last}

.impute_average_O <-
function(D, I, NAs, control)
    .make_relation_from_domain_and_incidence(D,
                                        t(.N.(.impute_average_WO(t(.N.(I)),
                                                                 length(NAs),
                                                                 diag = 0))))

.impute_average_O_first <-
function(D, I, NAs, control)
    .impute_any_O_first(D, I, NAs, control)

.impute_average_O_last <-
function(D, I, NAs, control)
    .impute_any_O_last(D, I, NAs, control)

## averaged incidences for W and O
.impute_average_WO <-
function(I, n, diag = 1)
{
    o <- do.call(order, split(I, col(I)))
    ro <- `[<-`(seq_along(o), o, seq_along(o))
    I <- I[o, o, drop = FALSE]

    mn <- ncol(I)
    m <- mn - n
    N <- .nsol_W(m, n)

    A <- I[1:m, 1:m]

    reps <- table(cumsum(!duplicated(A)))
    c <- length(reps)

    f1 <- (.nsol_W(c, n - 1) + .nsol_W(c + 1, n - 1)) / N
    fp <- f1 * (c + 1) / 2
    fi <- rep(seq_len(c) * f1, reps)

    B <- matrix(rep(fi, n), m, n)
    C <- matrix(rep(rev(fi), each = n), n, m)
    D <- matrix(fp, n, n)
    diag(D) <- diag
    rbind(cbind(A, B), cbind(C, D))[ro, ro, drop = FALSE]
}

## individual relations for G

.impute_any_G <-
function(D, I, NAs, control)
{
    if (is.null(control$n) || isTRUE(control$n == 1L))
        return(.impute_omit(D, I, NAs, control))
    if (is.na(control$n) || control$n == "all")
        control$n <- Inf

    ind <- which(is.na(I), arr.ind = TRUE)
    I <- list(I)
    ## Recursively generate all solutions by using 0 or 1 for the
    ## NA entries of I.
    splitter <- function(x, i, j) {
        y <- x
        ## By default we use 1 for the zero entries of M.
        x[i, j] <- 0
        y[i, j] <- 1
        list(x, y)
    }
    for (k in seq_len(nrow(ind))) {
        I <- do.call("c", lapply(I, splitter, ind[k, 1L], ind[k, 2L]))
        if (length(I) > control$n) {
            I <- lapply(I[seq_len(control$n)],
                        function(i) `[<-`(i, is.na(i), 0))
            break
        }
    }

    relation_ensemble(list = lapply(I, function(i)
                      .make_relation_from_domain_and_incidence(D, i)))
}

## individual relations for L{/first,/last}

.impute_any_L <-
function(D, I, NAs, control)
    .impute_any_WL(D, I, NAs, control, by = 2)

.impute_any_L_last <-
function(D, I, NAs, control, last = TRUE)
{
    if (is.null(control$n) || isTRUE(control$n == 1L))
        .make_relation_from_domain_and_incidence(D,
                                                 .impute_L_one(seq_along(NAs),
                                                               I, NAs, last))
    else {
        perm <- .permute(seq_along(NAs))
        if (!is.na(control$n) && control$n != "all" && control$n < length(perm))
            perm <- perm[seq_len(control$n)]
        relation_ensemble(list = lapply(perm,
                        function(i) .make_relation_from_domain_and_incidence(D,
                                    .impute_L_one(i, I, NAs, last))))
    }
}

.impute_any_L_first <-
function(D, I, NAs, control)
    .impute_any_L_last(D, I, NAs, control, last = FALSE)

## individual relations for W{/first,/last}

.impute_any_W <-
function(D, I, NAs, control)
    .impute_any_WL(D, I, NAs, control)

.impute_any_W_last <-
function(D, I, NAs, control)
{
    I[ ,NAs] <- 0
    I[NAs, ] <- 1
    .make_relation_from_domain_and_incidence(D, I)
}

.impute_any_W_first <-
function(D, I, NAs, control)
{
    I[NAs, ] <- 0
    I[, NAs] <- 1
    .make_relation_from_domain_and_incidence(D, I)
}

## individual relations for O{/first,/last}

.impute_any_O <-
function(D, I, NAs, control)
    t(!.impute_any_WL(D, t(.N.(I)), NAs, control, diag = 0))

.impute_any_O_last <-
function(D, I, NAs, control)
{
    I[ ,NAs] <- 0
    I[NAs, -NAs] <- 1
    diag(I) <- 1
    .make_relation_from_domain_and_incidence(D, I)
}

.impute_any_O_first <-
function(D, I, NAs, control)
{
    I[NAs, ] <- 0
    I[-NAs, NAs] <- 1
    diag(I) <- 1
    .make_relation_from_domain_and_incidence(D, I)
}

.impute_any_WL <-
function(D, I, NAs, control, by = 1, diag = 1)
{
    ## INIT:
    ## N ... objects
    ## M ... objects already ranked
    ## L ... objects to be ranked

    ## compute scores
    S <- colSums(I, na.rm = TRUE)
    S[NAs] <- NA

    N <- LABELS(D[[1L]], quote = FALSE)
    M <- rev(tapply(names(S), S, c))
    names(M) <- dimnames(M) <- NULL
    L <- setdiff(N, unlist(M))

    add_one <- function(x, e) {
        ## prepare slots
        slots <- c(list(c()),
                   unlist(lapply(x, function(i) list(i, c())),
                          recursive = FALSE))

        ## fill in element
        lapply(seq(1, length(slots), by = by), function(i) {
            tmp <- slots
            tmp[[i]] <- c(tmp[[i]], e)
            tmp[sapply(tmp, length) > 0L]
        })
    }

    ## FIXME: limit the number of solutions to control$n
    ret <- list(M)
    for (e in L)
        ret <- unlist(lapply(ret, add_one, e), recursive = FALSE)

    FUN <- if (diag == 0)
        function(i) reflexive_reduction(as.relation(ranking(i, D[[1L]])))
    else
        function(i) as.relation(ranking(i, D[[1L]]))

    if (is.null(control$n) || isTRUE(control$n == 1L))
        FUN(ret[[1L]])
    else {
        if (!is.na(control$n) && control$n != "all" && control$n < length(ret))
            ret <- ret[seq_len(control$n)]
        relation_ensemble(list = lapply(ret, FUN))
    }
}

## helper functions

.nsol_W <-
function(c, n)
{
    if (n < 1L) return(1L)
    x <- rep.int(1, n + 1L)
    v <- seq(from = c, to = c + n)
    while((len <- length(x)) > 1L) {
        x <- v[-len] * x[-len] + v[-1L] * x[-1L]
        v <- v[-len]
    }
    x
}

.missing_objects <-
function(I)
    which(apply(I, 1, function(i) all(is.na(i))))

.impute_L_one <-
function(perm, I, NAs, last = TRUE)
{
    if (last) {
        I[,NAs] <- 0
        I[NAs,] <- lower.tri(I, diag = TRUE)[NAs,]
    } else {
        I[NAs,] <- 0
        I[,NAs] <- upper.tri(I, diag = TRUE)[,NAs]
    }
    P <- c(seq_len(ncol(I) - length(NAs)), NAs[perm])
    I[P, P]
}
