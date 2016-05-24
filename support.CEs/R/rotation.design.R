rotation.design <- function(candidate.array = NULL, attribute.names, nalternatives,
                            nblocks, row.renames = TRUE, randomize = FALSE, seed = NULL) 
{
# Name: rotation.design
# Title: Creating a choice experiment design using the rotation or mix-and-match method
# Arguments:
#  candidate.array   A data frame containing an array created by the user.
#  attribute.names   A list of the names of attributes and levels.
#  nalternatives     An integer value describing the number of alternatives per choice set,
#                      excluding an opt-out alternative or common base alternative.
#  nblocks           An integer value describing the number of blocks into which
#                      a choice experiment design is divided.
#  row.renames       A logical variable describing whether or not the row names of a choice 
#                      experiment design created by this function are changed
#  randomize         If TRUE, the function executes the mix-and-match method;
#                      if FALSE, the function executes the rotation method.
#  seed              Seed for a random number generator.



# initial setting

    attribute.levels <- sapply(attribute.names, length)
    nattributes <- length(attribute.levels)

# search an array corresponding to the argument

    if (is.null(candidate.array) == TRUE) {
        OA <- oa.design(nlevels = attribute.levels, seed = seed)
    }
    else {
        OA <- candidate.array
    }

    nquestions <- nrow(OA)

# check number of blocks against number of rows of the design

    if ((nquestions%%nblocks) != 0) {
        cat("\n")
        cat("Your nblocks =", nblocks, "\n")
        cat("Number of rows of the design =", nrow(OA), "\n")
        stop("Number of blocks should be divisors of number of rows of the design.", call. = FALSE)
    }

# rotation design

    ALTS <- vector("list", nalternatives)
    ALTS[[1]] <- matrix(as.numeric(as.matrix(OA)), nrow = nquestions, ncol = nattributes)
    for (i in 2:nalternatives) {
        ALTS[[i]] <- matrix(0, nrow = nquestions, ncol = ncol(OA))
    }
    for (i in 2:nalternatives) {
        ALTS[[i]] <- sweep(ALTS[[(i - 1)]], 2, apply(ALTS[[(i - 1)]], 2, max), "%%") + 1
    }

# rename row

    if (row.renames != TRUE) { 
        OriginalRowNames <- row.names(OA)
        for (i in 1:nalternatives) {
            rownames(ALTS[[i]]) <- OriginalRowNames
        }
    }

# mix-and-Match design

    if (is.null(seed) == FALSE) {
        set.seed(seed)
    }

    if (randomize == TRUE) {
        repeat {
            for (i in 1:nalternatives) {
                ALTS[[i]] <- cbind(ALTS[[i]], sample(1:nquestions, nquestions))
                ALTS[[i]] <- ALTS[[i]][order(ALTS[[i]][, (nattributes + 1)]), ]
                ALTS[[i]] <- ALTS[[i]][, 1:nattributes]
            }
            chk <- FALSE
            for (i in 1:(nalternatives - 1)) {
                for (j in (i + 1):nalternatives) {
                    for (k in 1:nquestions) {
                        chk <- chk + all(ALTS[[i]][k, ] == ALTS[[j]][k, ])
                    }
                }
            }
            if (sum(chk) == 0) {
                break
            }
        }
    }

# randomize order of questions

    nquestions_nblocks <- nquestions / nblocks
    RndUnif <- runif(nquestions)
    for (i in 1:nalternatives) {
        ALTS[[i]] <- cbind(ALTS[[i]], RndUnif)
        ALTS[[i]] <- ALTS[[i]][order(ALTS[[i]][, (nattributes + 1)]), ]
        ALTS[[i]] <- ALTS[[i]][, 1:nattributes]
    }

# add BLOCK, QES, and ALT variables to ALTS[[i]]

    for (i in 1:nalternatives) {
        ALTS[[i]] <- cbind(rep(1:nblocks, each = nquestions_nblocks),
                           rep(1:nquestions_nblocks, nblocks),
                           rep(i, nquestions),
                           ALTS[[i]])
    }

# convert ALTS[[i]] into data.frame

    for (i in 1:nalternatives) {
        ALTS[[i]] <- data.frame(ALTS[[i]])
    }

# treat attribute columns as factors

    for (i in 1:nalternatives) {
        for (j in 1:nattributes) {
            ALTS[[i]][, (j + 3)] <- factor(ALTS[[i]][, (j + 3)])
        }
    }

# assign names to attributes/levels

    for (i in 1:nalternatives) {
        colnames(ALTS[[i]]) <- c("BLOCK", "QES", "ALT", names(attribute.names))
        for (j in 1:nattributes){
            levels(ALTS[[i]][, (j + 3)]) <- attribute.names[[j]]
        }
    }

# format output

    names(ALTS) <- paste("alt.", 1:nalternatives, sep = "")
    desinf <- list(nalternatives = nalternatives,
                   nblocks = nblocks,
                   nquestions = nquestions_nblocks,
                   nattributes = nattributes
                   )
    my.choice.experiment.design <- list(alternatives = ALTS,
                                        candidate = OA,
                                        design.information = desinf)
    class(my.choice.experiment.design) <- "cedes"

    return(my.choice.experiment.design)
}

