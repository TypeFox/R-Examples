Lma.design <- function(candidate.array = NULL, attribute.names, nalternatives,
                       nblocks, row.renames = TRUE, seed = NULL) 
{
# Name: Lma.design
# Title: Creating a choice experiment design using the L^MA method
# Arguments:
#  candidate.array   A data frame containing an array created by the user
#  attribute.names   A list of the names of attributes and levels
#  nalternatives     An integer value describing the number of alternatives per choice set,
#                      excluding an opt-out alternative or a common base alternative
#  nblocks           An integer value describing the number of blocks into which
#                      a choice experiment design is divided
#  row.renames       A logical variable describing whether or not the row names of a choice 
#                      experiment design created by this function are changed
#  seed              Seed for a random number generator.



# set variable

    # calculate number of levels for each attribute 
    attribute.levels <- rep(sapply(attribute.names, length), nalternatives)
    if (nblocks >= 2) {
        attribute.levels <- c(attribute.levels, nblocks)
    }

# create orthogonal array

    if (is.null(candidate.array) == TRUE) { # generate it using oa.design()
        OA <- oa.design(nlevels = attribute.levels, seed = seed)
    }
    else { # use arbitrary user-defined array
        OA <- candidate.array
    }

# set variables, matrix, and list

    nattributes <- length(attribute.names)           # number of attributes per alternative
    total.nattributes <- nattributes * nalternatives # total number of attribute
    nquestions <- nrow(OA)                           # number of questions
    nquestions_nblocks <- nquestions / nblocks       # number of question per block
    alt <- vector("list", nalternatives)             # list object "alt"

    if (is.null(candidate.array) == TRUE) {          # assign dummy column to orthogonal array
        if (nblocks == 1) {                          #  if number of blocks is one
            backupOA <- OA
            OA <- cbind(OA, DUMMY=rep(1, nquestions))
        }
    }

# randomize order of questions

    if (is.null(seed) == FALSE) { # set seed for random number generator
        set.seed(seed)
    }
    ALTS <- transform(OA, r = runif(nquestions)) # assign column containing uniform random number
    ALTS <- ALTS[order(ALTS[,(total.nattributes + 1)], ALTS$r), ] # randomize order of rows (questions)

# create alternatives

    if (nattributes == 1) { # number of attributes == 1
        for (i in 1:nalternatives){
            temp <- cbind(BLOCK = as.numeric(ALTS[, (total.nattributes + 1)]),
                          QES = rep(1:nquestions_nblocks, nblocks),
                          ALT = rep(i, nquestions),
                          ALTS[, i:(i + 1)])[, 1:4]
            colnames(temp)[4] <- names(attribute.names)
            levels(temp[, 4]) <- attribute.names[[1]]          
            alt[[i]] <- temp
        }
    }
    else { # number of attributes >= 2
        for (i in 1:nalternatives) {
            temp <- ALTS[, (1 + (i - 1) * nattributes):(i * nattributes)]
            colnames(temp) <- names(attribute.names)
            for (j in 1:nattributes){
                levels(temp[, j]) <- attribute.names[[j]]
            }
            alt[[i]] <- cbind(BLOCK = as.numeric(ALTS[, (total.nattributes + 1)]),
                              QES = rep(1:nquestions_nblocks, nblocks),
                              ALT = rep(i, nquestions),
                              temp)
        }
    }

# format output
    if (row.renames == TRUE) {
        for (i in 1:nalternatives) {
            rownames(alt[[i]]) <- c(1:nquestions)
        }
    }
    if (is.null(candidate.array) == TRUE) {
        if (nblocks == 1){
            OA <- backupOA
        }
    }
    ALTS <- list(NULL)
    for (i in 1:nalternatives) {
        ALTS[i] <- list(data.frame(alt[[i]]))
    }
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

