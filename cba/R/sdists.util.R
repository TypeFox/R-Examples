
###
### stuff for analyzing sequences
###
### ceeboo 2006, 2007

## Find the centroid (medoid) sequence(s), i.e. which have minimum 
## sum of distance among a collection of sequences.
##
## Alternativley, apply FUN to the distances and select the 
## distance with the minimum value of FUN (mean, median, etc.)
##
## Option 'unique' specifies to reduce the result set to a distinct
## set of sequences.
##
sdists.center <-
function(x, d = NULL, method = "ow", weight = c(1, 1, 0, 2),
    exclude = c(NA, NaN, Inf, -Inf), FUN = NULL, ..., unique = FALSE)
{
    if (is.null(d))
        d <- sdists(x, method = method, weight = weight,
            exclude = exclude)
    r <- 
    if (is.null(FUN)) 
        rowSums.dist(d) 
    else 
        apply(as.matrix(d), 1, FUN, ...)
    k <- which(r == min(r))
    r <- x[k]
    if (unique && length(r) > 1)
        r <- r[!duplicated(sapply(r, paste, collapse = ""))]
    r
}

## Compute a global alignment of a collection of sequences using
## the center star tree heuristic, i.e. 'c' is assumed to be the
## center and each remaining sequence 'x' is aligned in turn 
## against the (aligned) center. Spaces are inserted as needed.
##
## If transitive = TRUE the space symbols in the center are
## replaced by the (non-space) symbols in the current sequence
## the center sequence was aligned with. This results in a
## 'transitive' global alignment, i.e. each pair of sequences
## is implicitly aligned, too. However, this usually  results
## in spreading out the alignments (considerably).
##
## NOTE unfortunately, there may not exist a unique alignment
##      of a pair of sequences, so that the global 
##      alignment may not be unique either. In this case we
##      make a first or random choice.
##
sdists.center.align <-
function(x, center, method = "ow", weight = c(1, 1, 0, 2),
    exclude = c(NA, NaN, Inf, -Inf), break.ties = TRUE, transitive = FALSE,
    to.data.frame = FALSE)
{
    if (!is.list(x) && !is.character(x))
        stop("'x' not a list")
    if (missing(center)) {
        center <- sdists.center(x, method = method, weight = weight,
            exclude = exclude, unique = TRUE)
        k <- length(center)
        if (break.ties && k > 1)
            k <- sample(k, 1)
        else
            k <- 1
        center <- center[[k]]
    }
    n <- 0                          ## total number of ties
    r <- list()
    for (s in x) {
        a <- sdists.trace(center, s, method, weight, exclude)
        k <- length(a)
        if (break.ties && k > 1) {
            n <- n + k
            k <- sample(k, 1)       ## random choice
        } else
            k <- 1                  ## first
        center <- a[[k]][[1]]
        na <- is.na(center)
        if (any(na)) {
            center[na] <- ""
            if (length(r) > 0) {    ## update
                t <- gsub("[Dd]","?", names(a)[k])
                if (regexpr("[Ii]", t) > -1) 
                    r <- lapply(r, function(x) {
                        x <- .Call(R_sdists_align, x, center, t)[[1]]
                        x[is.na(x)] <- ""
                        x
                    })
            }
            if (transitive)
                center[na] <- a[[k]][[2]][na]
        }
        s <- a[[k]][[2]]
        s[is.na(s)] <- ""
        r <- c(r, list(s))
    }
    names(r) <- names(x)
    if (to.data.frame) {
        if (is.null(names(r)))
            names(r) <- seq_len(length(r))
        r <- data.frame(center = center, r, check.names = FALSE)
    } else {
        is.na(center) <- center == ""           ## recode space symbol
        names(center) <- seq(length(center))    ## add positional index
        r <- lapply(r, function(x) {
            is.na(x) <- x == ""
            names(x) <- names(center)
            x
        })
        attr(r, "center") <- center
        attr(r, "ties")   <- n
    }
    r
}

###
