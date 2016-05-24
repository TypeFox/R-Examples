### data structure for rankings

ranking <-
function(x, domain = NULL, decreasing = TRUE, complete = FALSE)
{
    ## some checks
    elements <- unlist(x)
    if (any(duplicated(elements)))
        stop("Elements must be unique.")

    if (is.null(domain))
        domain <- elements
    domain <- as.set(domain)

    elements <- LABELS(elements, quote = FALSE)
    n <- LABELS(domain, quote = FALSE)
    if (any(duplicated(n)))
        stop("Need unique element labels for domain.")
    if (any(sapply(n, nchar) < 1L))
        stop("All domain elements need to be named.")

    if (!all(elements %in% n))
        stop("Out-of-domain elements.")

    ## put missing elements at the end
    n <- c(elements, n[! n %in% elements])

    ## create scores
    x <- as.list(x)
    SEQ <- if (!decreasing) seq_along else function(x) rev(seq_along(x))
    ret <- rep(SEQ(x), sapply(x, length))[1 : length(n)]
    names(ret) <- n

    ## complete scores, if needed
    if (complete && any(nas <- is.na(ret)))
        ret[nas] <- if (!decreasing) max(ret, na.rm = TRUE) + 1 else 0

    ## return ranking object
    .structure(list(domain = domain,
                    scores = ret[LABELS(domain, quote = FALSE)],
                    decreasing = decreasing),
               class = "ranking")
}

print.ranking <-
function(x, ...)
{
    scores <- x$scores
    nas <- is.na(scores)
    if (any(nas)) {
        na_elements <- names(scores)[nas]
        scores <- scores[!nas]
    }

    classes <- tapply(names(scores), scores, c)
    SYM <- if (x$decreasing) {
        classes <- rev(classes)
        " > "
    } else " < "

    classes <- lapply(classes, function(i) {
        s <- paste(i, collapse = " ~ ")
        if (length(i) > 1L)
            s <- paste("[", s, "]")
        s
    })

    writeLines(paste(classes, collapse = SYM))
    if (any(nas))
        writeLines(sprintf("Missing elements: %s",
                           paste(na_elements, collapse = " ")))

    invisible(x)
}

as.relation.ranking <-
function(x, ...)
{
## include decreasing information in meta data to allow
### as.ranking.relation to recover the original ranking structure.

    meta <- list(is_decreasing = x$decreasing)
    .make_relation_from_domain_and_scores(x$domain, x$scores, meta)
}

as.ranking <-
function(x, ...)
    UseMethod("as.ranking")

as.ranking.default <-
function(x, ...)
    ranking(x, ...)

as.ranking.ranking <-
function(x, ...)
    x

as.ranking.relation <-
function(x, ...)
{
    dec <- relation_property(x, "is_decreasing")
    if (is.null(dec))
        dec <- TRUE
    .as.ranking.relation(.get_representation(x), dec)
}

.as.ranking.relation <-
function(x, decreasing)
    UseMethod(".as.ranking.relation")

.as.ranking.relation.relation_by_domain_and_incidence <-
function(x, decreasing)
{
    I <- .incidence(x)
    ret <- colSums(I, na.rm = TRUE)
    ret[.missing_objects(I)] <- NA
    storage.mode(ret) <- "integer"
    labs <- .domain(x)[[1L]]
    names(ret) <- LABELS(labs)
    .make_ranking_by_domain_and_scores(labs, ret, decreasing)
}

.as.ranking.relation.relation_by_domain_and_scores <-
function(x, decreasing)
    .make_ranking_by_domain_and_scores(x$domain, x$scores, decreasing)

.make_ranking_by_domain_and_scores <-
function(domain, scores, decreasing = TRUE)
    .structure(list(domain = domain,
                    scores = scores,
                    decreasing = decreasing),
               class = "ranking")

is.ranking <-
function(x)
    inherits(x, "ranking")

rev.ranking <-
t.ranking <-
function(x)
{
    x$scores <- max(x$scores, na.rm = TRUE) + 1 - x$scores
    x
}
