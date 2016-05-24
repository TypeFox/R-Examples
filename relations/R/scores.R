### Relation scores.

relation_scores <-
function(x, method, normalize = FALSE, ...)
    UseMethod("relation_scores")

relation_scores.relation_ensemble <-
function(x,
         method = c("Borda", "Kendall", "differential", "Copeland"),
         normalize = FALSE, weights = 1, ...)
{
    method <- match.arg(method)
    if(!.is_ensemble_of_endorelations(x))
        stop("Relation scores are only available for (ensembles of) endorelations.")
    labs <- LABELS(.domain(x)[[1L]], quote = FALSE)
    weights <- rep(weights, length.out = length(x))
    I <- .weighted_sum_of_arrays(lapply(x, .incidence), weights)
    ret <- switch(method,
                  Borda =, Kendall = colSums(I),
                  Copeland =, differential = colSums(I) - rowSums(I)
                  )

    names(ret) <- labs
    if(normalize)
        ret / sum(ret)
    else
        ret
}


relation_scores.relation <-
function(x,
         method = c("ranks", "Barthelemy/Monjardet", "Borda",
                    "Kendall", "Wei", "differential", "Copeland"),
         normalize = FALSE, ...)
{
    method <- match.arg(method)
    if(!relation_is_endorelation(x))
        stop("Relation scores are only available for endorelations.")
    labs <- LABELS(.domain(x)[[1L]], quote = FALSE)

    ## <NOTE>
    ## When adding .relation_scores_FOO(x, ...) method functions, we
    ## might want to pass x as the relation itself rather than its
    ## incidence.
    I <- .incidence(x)
    ## </NOTE>

    ret <- switch(method,
                  ranks = .incidence_scores_ranks(I, ...),
                  "Barthelemy/Monjardet" = {
                      ## Use formula in Barthelemy & Monjardet, p. 258.
                      ## See also
                      ## http://mathworld.wolfram.com/ScoreSequence.html.
                      (colSums(I * (1 - t(I))) + colSums(I) - 1) / 2
                  },
                  Borda =, Kendall = colSums(I),
                  Copeland =, differential = colSums(I) - rowSums(I),
                  Wei = {
                      ## <FIXME>
                      ## Cook & Kress use "preference matrices", so to
                      ## be consistent, take complement?
                      abs(Re(eigen(1 - I)$vectors[, 1L]))
                      ## </FIXME>
                  })

    names(ret) <- labs
    if(normalize)
        ret / sum(ret)
    else
        ret
}

.incidence_scores_ranks <-
function(x, decreasing = TRUE)
{
    n <- ncol(x)
    if(decreasing)
        (n + 1 + rowSums(x) - colSums(x)) / 2
    else
        (n + 1 + colSums(x) - rowSums(x)) / 2
}

