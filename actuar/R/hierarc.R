### ===== actuar: An R Package for Actuarial Science =====
###
### Hierarchical credibility model calculations
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Louis-Philippe Pouliot, Tommy Ouellet.

hierarc <- function(ratios, weights, classification,
                    method = c("Buhlmann-Gisler", "Ohlsson", "iterative"),
                    tol = sqrt(.Machine$double.eps), maxit = 100,
                    echo = FALSE)
{
    ## === HANDS ON THE DATA ===
    ##
    ## Arguments 'ratios' and 'weights' must be matrices of real
    ## numbers, whereas 'classification' must be a matrix of integers
    ## giving the affiliation of each entity in the portfolio.
    nlevels <- ncol(classification)              # number of levels
    nlevels1p <- nlevels + 1L                    # frequently used

    ## To symmetrize further calculations, bind a column of ones
    ## representing the affiliation to the global portfolio.
    classification <- cbind(pf = 1L, classification)

    ## If weights are not specified, use equal weights.
    if (missing(weights))
    {
        if (any(is.na(ratios)))
            stop("missing ratios not allowed when weights are not supplied")
        array(1, dim(ratios))       # matrix of ones
    }

    ## Sanity check if weights and ratios correspond.
    if (!identical(which(is.na(ratios)), which(is.na(weights))))
        stop("missing values are not in the same positions in 'weights' and in 'ratios'")

    ## === NUMBER OF NODES AND SPLITTING FACTORS ===
    ##
    ## Future computation of per level summaries will require a set of
    ## factors based on the number of nodes in each level. An example
    ## will best explain what is achieved here: suppose there are two
    ## sectors; sector 1 has 3 units; sector 2 has 2 units. To make
    ## per sector summaries, the following factors can be used to
    ## split the unit data: 1 1 1 2 2.
    ##
    ## Generating such factors first requires to know the number of
    ## nodes at each level in a format identical to the 'nodes'
    ## argument of simul(). [In the previous example, the number of
    ## nodes would be 'list(2, c(3, 2))'.] Then, the factors are
    ## obtained by repeating a sequence the same length as the number
    ## of nodes at one level [2] according to the number of nodes at
    ## the level below [c(3, 2)].
    ##
    ## 0. Initialization
    fnodes <- nnodes <- vector("list", nlevels)

    ## 1. Calculation of the number of nodes: the main idea is to
    ## create a unique factor for each node using interaction()
    ## recursively on the columns of 'classification'. We can do
    ## something simpler for the lowest level (the entities), since we
    ## know the combinations of indexes to all be different at this
    ## level.
    fx <- vector("list", nlevels1p)
    fx[[nlevels1p]] <- factor(classification[, nlevels1p]) # entity level

    for (i in nlevels:1L)
    {
        ## Function 'interaction' expects its arguments separately or
        ## as a list, hence the lapply() below.
        fx[[i]] <- as.integer(interaction(lapply(seq.int(i),
                                                 function(j) classification[, j]),
                                          drop = TRUE))
        ## 'as.vector' below is used to get rid of names
        nnodes[[i]] <- as.vector(sapply(split(fx[[i + 1]], fx[[i]]),
                                        function(x) length(unique(x))))
    }

    ## 2. Generation of the factors. Following the rule described
    ## above, this could simply be
    ##
    ##   fnodes <- lapply(nnodes, function(x) rep(seq_along(x), x))
    ##
    ## However, this will not work if rows of data are not sorted per
    ## level. (In the example above, if the data of unit 3 of sector 1
    ## is at the bottom of the matrix of data, then the factors need
    ## to be 1 1 2 2 1.)
    ##
    ## The solution is actually simple: converting the entity level
    ## factors ('fx[[nlevels]]') to integers will assure that any
    ## summary made using these factors will be sorted. This done, it
    ## is possible to use the command above for the upper levels.
    fnodes[[nlevels]] <- as.integer(fx[[nlevels]])
    fnodes[-nlevels] <- lapply(nnodes[-nlevels],
                               function(x) rep(seq_along(x), x))

    ## === PER ENTITY SUMMARIES ===
    ##
    ## Individual weighted averages. It could happen that an entity
    ## has no observations, for example when applying the model on
    ## claim amounts. In such a situation, put the total weight of the
    ## entity and the weighted average both equal to zero. That way,
    ## the premium will be equal to the credibility weighted average,
    ## as it should, but the entity will otherwise have no
    ## contribution in the calculations.
    weights.s <- rowSums(weights, na.rm = TRUE)
    ratios.w <- ifelse(weights.s > 0, rowSums(weights * ratios, na.rm = TRUE) / weights.s, 0)

    ## === EFFECTIVE NUMBER OF NODES ===
    ##
    ## Given the possibility to have whole levels with no data, as
    ## explained above, it is necessary to count the *effective*
    ## number of nodes in each level, that is the number of nodes with
    ## data. This comes this late since it relies on 'weights.s'.
    ##
    ## Object 'eff.nnodes' is in every respect equivalent to 'nnodes'
    ## except that each element of the list is a vector of the number of
    ## non "empty" nodes for each classification of the level
    ## above.
    eff.nnodes <- vector("list", nlevels)
    w <- weights.s
    for (i in nlevels:1L)
    {
        eff.nnodes[[i]] <- tapply(w, fnodes[[i]], function(x) sum(x > 0))
        w <- tapply(w, fnodes[[i]], sum) # running totals
    }

    ## === DENOMINATORS OF VARIANCE ESTIMATORS ===
    ##
    ## The denominators for all the variance estimators never
    ## change. The denominator at one level is equal to the total
    ## number of nodes at that level minus the total number of nodes
    ## at the level above. At the lowest level (the denominator of
    ## s^2), this is
    ##
    ##   number of (non NA) ratios - (effective) number of entities.
    ##
    ## The number of (non missing) ratios is not included in
    ## 'eff.nnodes'.  For the portfolio level, the denominator is
    ##
    ##   (effective) number of "sectors" - 1
    ##
    ## The 1 neither is included in 'eff.nnodes'.
    denoms <- diff(c(1L, sapply(eff.nnodes, sum), sum(!is.na(ratios))))

    ## Final sanity checks
    if (any(!denoms))
        stop("there must be at least two nodes at every level")
    if (ncol(ratios) < 2L)
        stop("there must be at least one node with more than one period of experience")

    ## === ESTIMATION OF s^2 ===
    s2 <-  sum(weights * (ratios - ratios.w)^2, na.rm = TRUE) /
        denoms[nlevels1p]

    ## === ESTIMATION OF THE OTHER VARIANCE COMPONENTS ===
    ##
    ## Create vectors to hold values to be computed at each level
    ## (from portfolio to entity), namely: the total node weights, the
    ## node weighted averages, the between variances and the node
    ## credibility factors.
    ##
    ## Only credibility factors are not computed for the portfolio
    ## level, hence this list is one shorter than the others.
    tweights <- vector("list", nlevels1p)       # total level weights
    wmeans <- vector("list", nlevels1p)         # weighted averages
    b <- c(numeric(nlevels), s2)                # variance estimators
    cred <- vector("list", nlevels)             # credibility factors

    ## Values already computed at the entity level.
    tweights[[nlevels1p]] <- as.vector(weights.s);
    wmeans[[nlevels1p]] <- as.vector(ratios.w);

    ## The unbiased variance estimators are evaluated first as they will
    ## be used as starting values for the iterative part below.
    ##
    ## At the entity level: node weight is given by the natural
    ## weight, weighted averages use the natural weights.
    ##
    ## Level above the entity: node weight is the sum of the natural
    ## weights at the level below, weighted averages use the natural
    ## weights.
    ##
    ## All upper levels: node weight is the sum of the credibility
    ## factors at the level below, weighted averages use credibility
    ## factors from previous level.
    ##
    ## Buhlmann-Gisler estimators truncate the per node variance
    ## estimates to 0 before taking the mean, whereas the Ohlsson
    ## estimators do not make any truncation.
    method <- match.arg(method)

    if (method == "Buhlmann-Gisler")
        bexp <- expression(b[i] <- mean(pmax(ifelse(ci != 0, bi/ci, 0), 0), na.rm = TRUE))
    else                                # Ohlsson
        bexp <- expression(b[i] <- sum(bi, na.rm = TRUE) / sum(ci, na.rm = TRUE))

    for (i in nlevels:1L)
    {
        ## Total weight of the level as per the rule above.
        tweights[[i]] <- as.vector(tapply(tweights[[i + 1L]], fnodes[[i]], sum))

        ## Calculation of the weighted averages of the level. Before
        ## the between variance is estimated, these use the total
        ## weights calculated above.
        wmeans[[i]] <-
            ifelse(tweights[[i]] > 0,
                   as.vector(tapply(tweights[[i + 1L]] * wmeans[[i + 1L]],
                                    fnodes[[i]],
                                    sum) / tweights[[i]]),
                   0)

        ## Latest non-zero between variance estimate -- the one used
        ## in the estimator and in the credibility factors.
        between <- b[b != 0][1L]

        ## Calculation of the per node variance estimate.
        bi <- as.vector(tapply(tweights[[i + 1L]] *
                                (wmeans[[i + 1L]] - wmeans[[i]][fnodes[[i]]])^2,
                                fnodes[[i]],
                                sum)) -
                                    (eff.nnodes[[i]] - 1) * between
        ci <- tweights[[i]] -
            as.vector(tapply(tweights[[i + 1L]]^2, fnodes[[i]], sum)) / tweights[[i]]

        ## The final estimate is the average of all the per node estimates.
        eval(bexp)

        ## Calculation of the credibility factors. If these are
        ## non-zero, the total weights for the current level are
        ## replaced by the sum of the credibility factors and the
        ## weighted averages are recomputed with these new weights.
        #if (max(bu[i], 0))           # don't compute negative factors!
        if (b[i])
        {
            cred[[i]] <- 1/(1 + between/(b[i] * tweights[[i + 1L]]))
            tweights[[i]] <- as.vector(tapply(cred[[i]], fnodes[[i]], sum))
            wmeans[[i]] <-
                ifelse(tweights[[i]] > 0,
                       as.vector(tapply(cred[[i]] * wmeans[[i + 1L]],
                                        fnodes[[i]],
                                        sum) / tweights[[i]]),
                       0)
        }
        else
            cred[[i]] <- numeric(sum(nnodes[[i]]))
    }

    ## Iterative estimation of the structure parameters.
    ##
    ## At the entity level: total weight is the sum of the natural
    ## weights, weighted averages use the natural weights and between
    ## variance is s^2.
    ##
    ## All upper levels: total weight is the sum of the credibility
    ## factors of the level below, weighted averages use credibility
    ## factors, between variance estimated recursively and credibility
    ## factor use total weight of the level, between variance of the
    ## level below (hence the within variance) and between variance of
    ## the current level.
    if (method == "iterative")
    {
        b <- pmax(b, 0)            # truncation for starting values
        if (any(head(b, -1L) > 0)) # at least one non-zero starting value
            .External("actuar_do_hierarc", cred, tweights, wmeans, fnodes, denoms,
                      b, tol, maxit, echo)
    }

    ## Results
    structure(list(means = wmeans,
                   weights = tweights,
                   unbiased = if (method != "iterative") b,
                   iterative = if (method == "iterative") b,
                   cred = cred,
                   nodes = nnodes,
                   classification = classification[, -1L],
                   ordering = fnodes),
              class = "hierarc",
              model = "hierarchical")
}

predict.hierarc <- function(object, levels = NULL, newdata, ...)
{
    ## The credibility premium of a node at one level is equal to
    ##
    ##   p + z * (m - p)
    ##
    ## where 'p' is the credibility premium of the level above (or the
    ## collective premium for the portfolio), 'z' is the credibility
    ## factor of the node, and 'm' is the weighted average of the
    ## node.
    fnodes <- object$ordering
    cred <- object$cred
    means <- object$means
    nlevels <- length(object$nodes)
    level.names <- names(object$nodes)

    if (is.null(levels))
        levels <- seq_len(nlevels)
    if (any(is.na(levels)) || !is.numeric(levels))
        stop("invalid level number")
    n <- max(levels)

    res <- vector("list", n)

    ## First level credibility premiums
    res[[1L]] <- means[[1L]] + cred[[1L]] * (means[[2L]] - means[[1L]])

    for (i in seq(2, length.out = n - 1))
    {
        p <- res[[i - 1]][fnodes[[i]]]
        res[[i]] <- p + cred[[i]] * (means[[i + 1]] - p)
    }

    structure(res[levels], names = level.names[levels], ...)
}
