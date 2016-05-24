# A constructor function for a bdsmatrix with blocks of ones
#   It is used for nested effects
# group: the grouping variable
# id   : the eventual dimnames
#
bdsBlock <- function(id, group) {
    if (any(is.na(group))) stop ("Missing group indicator not allowed")
    blocksize <- as.vector(table(group))
    id <- id[order(group)]  # resort the data in group order
    temp <- sum(blocksize * (blocksize+1)/2)
    bdsmatrix(blocksize=blocksize, blocks=rep(1.0, temp),
              dimnames=list(id, id))
    }
