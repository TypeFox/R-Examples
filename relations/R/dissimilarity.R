### Relation metrics (dissimilarities).

## We use the CLUE approach: have a front-end
##   relation_dissimilarity <-
##       function(x, y = NULL, method = "symdiff", ...)
## which returns dissimilarities or cross-dissimilarities, eventually
## with classes and subscript methods as in CLUE, and a mechanism for 
## registering dissimilarity methods.

relation_dissimilarity <-
function(x, y = NULL, method = "symdiff", ...)
{
    x <- as.relation_ensemble(x)

    ## Be nice.
    if(is.character(y) || is.function(y)) {
        method <- y
        y <- NULL
    }

    known_methods <-
        list(symdiff =
             c(".relation_dissimilarity_symdiff",
               "symmetric difference distance"),
             SD =
             c(".relation_dissimilarity_symdiff",
               "symmetric difference distance"),
             CKS =
             c(".relation_dissimilarity_CKS",
               "Cook-Kress-Seiford distance"),
             CS =
             c(".relation_dissimilarity_CS",
               "Cook-Seiford distance"),
             score = 
             c(".relation_dissimilarity_score",
               "score-based distance"),
             manhattan =
             c(".relation_dissimilarity_manhattan",
               "Manhattan distance"),
             euclidean =
             c(".relation_dissimilarity_euclidean",
               "Euclidean distance"),
             Jaccard =
             c(".relation_dissimilarity_Jaccard",
               "Jaccard distance")
             )
    if(is.character(method)) {
        ## Hopefully of length one, add some tests eventually ...
        if(is.na(ind <- pmatch(method, names(known_methods))))
            stop(gettextf("Method '%s' is not a valid dissimilarity method.",
                          method),
                 domain = NA)
        method <- get(known_methods[[ind]][1L])
        method_name <- known_methods[[ind]][2L]
    }
    else if(is.function(method))
        method_name <- "user-defined method"
    else
        stop("Invalid 'method' argument.")

    if(!is.null(y)) {
        y <- as.relation_ensemble(y)
        D <- relation_domain(x)
        if(!.domain_is_equal(relation_domain(y), D))
            stop("All relations must have the same domain.")
        ## Build a cross-proximity object of cross-dissimilarities.
        ## <FIXME>
        ## Not yet: if we don't want to require clue, all we can return
        ## is a matrix ...
        d <- matrix(0, length(x), length(y))
        for(j in seq_along(y))
            d[, j] <- sapply(x, method, y[[j]], ...)
        dimnames(d) <- list(names(x), names(y))
        attr(x, "description") <- method_name
        return(d)
        ## </FIXME>
    }

    ## Otherwise, build a proximity object of dissimilarities.
    ## <FIXME>
    ## Not yet: if we don't want to require clue, all we can return
    ## is a dist object ...
    n <- length(x)
    d <- vector("list", length = n - 1L)
    ind <- seq_len(n)
    while(length(ind) > 1L) {
        j <- ind[1L]
        ind <- ind[-1L]
        d[[j]] <- sapply(x[ind], method, x[[j]], ...)
    }
    ## </FIXME>
    ## Grr ... see clue:::.dist_from_vector().
    .structure(unlist(d), Size = n, Labels = names(x),
               class = "dist", description = method_name)
}

.relation_dissimilarity_symdiff <-
function(x, y, na.rm = FALSE)
{
    if(!identical(relation_is_crisp(x), TRUE) ||
       !identical(relation_is_crisp(y), TRUE))
        stop("Not implemented.")        # David said so ...
    .incidence_dissimilarity_symdiff(relation_incidence(x),
                                     relation_incidence(y),
                                     na.rm = na.rm)
}

.incidence_dissimilarity_symdiff <-
function(x, y, na.rm = FALSE)
    sum(abs(x - y), na.rm = na.rm)

.relation_dissimilarity_CKS <-
function(x, y)
{
    ## Wade D. Cook and Moshe Kress and Lawrence M. Seiford
    ## Information and preference in partial orders: a bimatrix
    ## representation.
    ## Psychometrika 51/2. 197-207.
    ## Unique paired comparison metric for partial rankings (Definition
    ## 2.1) under the assumption that indifference is the centroid
    ## between strict preferences and incomparability.
    ## <NOTE>
    ## Originally only defined between (partial) rankings, but
    ## applicable more generally.
    ## </NOTE>
    .incidence_dissimilarity_CKS(relation_incidence(x),
                                 relation_incidence(y))
}

.incidence_dissimilarity_CKS <-
function(x, y)
{
    I_x <- pmax(x, t(x)) 
    P_x <- pmin(t(x), 1 - x)
    I_y <- pmax(y, t(y))
    P_y <- pmin(t(y), 1 - y)
    sum(abs((I_x - I_y)[!upper.tri(I_x)])) + sum(abs(P_x - P_y))
}

.relation_dissimilarity_CS <-
function(x, y)
{
    ## Wade D. Cook and Lawrence M. Seiford
    ## Priority Ranking and Consensus Formation
    ## Management Science 24/16, 1721--1732
    ## Ordinal ranking in the sense of Kendall CHECK
    ## Equivalent to: complete and transitive?
    ## I.e., weak order aka preference ...
    ## <NOTE>
    ## Originally only defined between (complete) rankings, but
    ## applicable more generally.
    ## </NOTE>
    .incidence_dissimilarity_CS(relation_incidence(x),
                                relation_incidence(y))
}

.incidence_dissimilarity_CS <-
function(x, y)
    sum(abs(.incidence_scores_ranks(x) - .incidence_scores_ranks(y)))

.relation_dissimilarity_score <-
function(x, y, score = NULL, Delta = 1)
{
    ## Score-based distance Delta(score(x), score(y)).
    ## If score is NULL, the default relation_score() is used.
    ## If it is a character string, relation_scores(x, score) is used.
    ## Otherwise, it must be a function (the score function itself).
    ## If Delta is a number p (\ge 1), the p-norm is used as Delta().
    ## Otherwise, it must be a fucntion.
    ## The generalized Cook-Seiford dissimilarity is a special case of a
    ## score-based distance (corresponding to the defaults).
    if(is.null(score)) {
        s_x <- relation_scores(x)
        s_y <- relation_scores(y)
    } else if(is.character(score) && (length(score) == 1L)) {
        s_x <- relation_scores(x, score)
        s_y <- relation_scores(y, score)
    } else if(is.function(score)) {
        s_x <- score(x)
        s_y <- score(y)
    }
    else
        stop("Invalid 'score' argument.")
    
    if(is.numeric(p <- Delta) && (length(p) == 1L) && (p >= 1)) {
        Delta <- function(u, v) sum(abs(u - v) ^ p) ^ (1 / p)
    } else if(!is.function(Delta))
        stop("Invalid 'Delta' argument.")

    Delta(s_x, s_y)
}

.relation_dissimilarity_manhattan <-
function(x, y, na.rm = FALSE)
    .incidence_dissimilarity_manhattan(relation_incidence(x),
                                       relation_incidence(y),
                                       na.rm = na.rm)

.incidence_dissimilarity_manhattan <-
function(x, y, na.rm = FALSE)    
    sum(abs(x - y), na.rm = na.rm)

.relation_dissimilarity_euclidean <-
function(x, y, na.rm = FALSE)
    .incidence_dissimilarity_euclidean(relation_incidence(x),
                                       relation_incidence(y),
                                       na.rm = na.rm)

.incidence_dissimilarity_euclidean <-
function(x, y, na.rm = FALSE)
    sqrt(sum((x - y) ^ 2, na.rm = na.rm))

.relation_dissimilarity_Jaccard <-
function(x, y, na.rm = FALSE)
    .incidence_dissimilarity_Jaccard(relation_incidence(x),
                                     relation_incidence(y),
                                     na.rm = na.rm)
## One could also use
##   1 - gset_similarity(relation_graph(x), relation_graph(y))
## but proceeding directly should be more efficient.

.incidence_dissimilarity_Jaccard <-
function(x, y, na.rm = FALSE)
{
    if(identical(all(x == 0, na.rm = na.rm), TRUE) &&
       identical(all(y == 0, na.rm = na.rm), TRUE))
        return(0)
    1 - sum(.T.(x, y), na.rm = na.rm) / sum(.S.(x, y), na.rm = na.rm)
}
