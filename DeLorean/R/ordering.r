#' Check that it is a valid ordering
#'
#' @param ordering The ordering
#' @param full.check Perform a full check
#'
ordering.is.valid <- function(ordering, full.check=FALSE) {
    if (! is.integer(ordering)
        || any(ordering <= 0)
        || any(ordering > length(ordering)))
    {
        return(FALSE)
    }
    # Check we have every item in the ordering
    if (full.check) {
        a <- rep(FALSE, length(ordering))
        a[ordering] <- TRUE
        if (! all(a)) {
            return(FALSE)
        }
    }
    return(TRUE)
}


#' Move one item in an ordering and shift the other items.
#'
#' @param ordering The ordering.
#' @param from The index of the item to move.
#' @param to The index it should be moved to.
#'
ordering.move <- function(ordering, from, to) {
    if (from != to) {
        item <- ordering[from]
        if (from < to - 1) {
            ordering[from:(to-1)] <- ordering[(from+1):to]
        } else if (from > to + 1) {
            ordering[(to+1):from] <- ordering[to:(from-1)]
        } else {
            # just swap them
            ordering[from] <- ordering[to]
        }
        ordering[to] <- item
    }
    ordering
}


#' Move a block in an ordering and shift the other items.
#'
#' @param ordering The ordering.
#' @param from The start index of the block to move.
#' @param width The width of the block to move.
#' @param to The index it should be moved to.
#' @param reverse Reverse the block?
#'
ordering.block.move <- function(ordering, from, width, to, reverse=FALSE) {
    if (from == to & ! reverse) {
        # Do nothing
        return(ordering)
    }
    # Length of whole ordering
    N <- length(ordering)
    # Check parameters
    stopifnot(width >= 1)
    stopifnot(from + width - 1 <= N)
    stopifnot(to + width - 1 <= N)
    # Block
    block <- ordering[from:(from+width-1)]
    if (reverse) {
        block <- rev(block)
    }
    # Move the other block
    if (from < to) {
        other.block <- ordering[(from+width):(to+width-1)]
        ordering[from:(to-1)] <- other.block
    } else if (from > to) {
        other.block <- ordering[to:(from-1)]
        ordering[(to+width):(from+width-1)] <- other.block
    }
    stopifnot(length(other.block) == abs(from-to))
    # Move the target block
    ordering[to:(to+width-1)] <- block
    stopifnot(ordering.is.valid(ordering))
    ordering
}


#' Randomly move one item in an ordering to another location
#'
#' @param ordering The ordering.
#'
ordering.random.move <- function(ordering) {
    .sample <- sample(length(ordering), 2)
    ordering.move(ordering, .sample[1], .sample[2])
}


#' Randomly move a block in an ordering to another location
#'
#' @param ordering The ordering
#' @param max.width The maximum width of the block
#'
ordering.random.block.move <- function(ordering, max.width=4) {
    width <- sample(max.width, 1)
    .sample <- sample(length(ordering) - width + 1, 2)
    reverse <- runif(1) > .5
    ordering.block.move(ordering, .sample[1], width, .sample[2], reverse)
}


#' Find a good ordering in the sense that some function is
#' locally maximised.
#'
#' @param fn A function to maximise.
#' @param ordering The permutation (ordering) to start from.
#'
ordering.maximise <- function(ordering, fn) {
    stopifnot(! is.null(ordering))
    while (TRUE) {
        ordering.new <- ordering.improve(fn, ordering)
        # If we didn't improve, then break the loop
        if (all(ordering.new == ordering)) {
            break
        } else {
            ordering <- ordering.new
            message('Improved ordering: f(ordering) = ', fn(ordering))
        }
    }
    stopifnot(! is.null(ordering))
    ordering
}


#' Metropolis-Hastings on orderings.
#'
#' @param ordering Initial ordering
#' @param log.likelihood Log likelihood function
#' @param proposal.fn Proposal function
#' @param iterations Number of iterations
#' @param thin Thinning parameter
#'
ordering.metropolis.hastings <- function(
    ordering,
    log.likelihood,
    proposal.fn=ordering.random.move,
    iterations=1000,
    thin=1)
{
    # Number of samples including initial
    num.samples <- floor(iterations/thin)+1
    # Create and initialise containers to store chain and likelihoods
    chain <- array(as.integer(0), dim=c(num.samples, length(ordering)))
    lls <- rep(0, num.samples)
    chain[1,] <- ordering
    last.ll <- lls[1] <- log.likelihood(ordering)
    # Set up counters and indices
    accepted <- 0
    sample.idx <- 1
    # Calculate random numbers all at once
    probs <- runif(iterations)
    # For each iteration
    for (i in 1:iterations) {
        # Get a new proposal state
        proposal <- proposal.fn(ordering)
        this.ll <- log.likelihood(proposal)
        # Accept/reject step
        if (probs[i] < exp(this.ll - last.ll)) {
            accepted <- accepted + 1
            ordering <- proposal
            last.ll <- this.ll
        }
        # Thin: only keep some samples
        if (0 == i %% thin) {
            sample.idx <- sample.idx + 1
            chain[sample.idx,] <- ordering
            lls[sample.idx] <- this.ll
        }
    }
    stopifnot(sample.idx == num.samples)
    # Return all state as a list
    list(chain = coda::mcmc(chain, start=1, end=iterations+1, thin=thin),
         log.likelihoods = lls,
         acceptance.rate = accepted / iterations)
}


#' Improve the ordering in the sense that some function is
#' maximised.
#'
#' @param fn A function to maximise.
#' @param ordering The permutation (ordering) to start from.
#'
ordering.improve <- function(fn, ordering) {
    Reduce(
        # Reduction function
        function(ordering, from) {
            scores <- sapply(1:length(ordering),
                            function(to) fn(ordering.move(ordering, from, to)))
            best.to <- which.max(scores)
            ordering.move(ordering, from, best.to)
        },
        # Choose all the elements in a random order
        sample(length(ordering)),
        # Initialise the reduction with the initial ordering
        init=ordering)
}


#' Invert the ordering
#'
#' @param ordering The permutation (ordering) to invert.
#'
ordering.invert <- function(ordering) {
    result <- rep(0, length(ordering))
    for (n in 1:length(result)) {
        result[ordering[n]] <- n
    }
    result
}


#' Test ordering score: sum every time consecutive items are in
#' order.
#'
#' @param ordering The permutation (ordering) to test.
#'
ordering.test.score <- function(ordering) {
    sum(sapply(1:(length(ordering)-1),
               function(n) as.integer(ordering[n] < ordering[n+1])))
}
