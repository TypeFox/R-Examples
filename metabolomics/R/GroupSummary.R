GroupSummary <- function(inputdata)
{
    # Get groups information
    groups <- factor(inputdata[, 1], levels=unique(inputdata[, 1]))
    # Get levels for groups
    grp_levs <- levels(groups)

    # Initialise an empty vector and create a vector of group lengths
    grp_len <- c()
    for (ii in 1:length(grp_levs)) {
        grp_len <- c(
            # Use the initial vector,
            grp_len,
            # and then append the number of samples in the group
            length(which(groups == levels(groups)[ii]))
        )
    }

    #
    #    Split the matrix by group
    #
    new_mats <- c()
    for (ii in 1:length(grp_levs)) {
        new_mats[ii] <- list(inputdata[which(groups == levels(groups)[ii]), ])
    }

    #
    #    Calculate the means
    #
    # For each matrix, calculate the averages and std per column
    submeans <- c()
    substd <- c()
    # Preallocate a matrix for the means
    std <- means <- matrix(
        nrow=length(grp_levs),
        ncol=length(colnames(inputdata[, -1])),
        dimnames=list(grp_levs, colnames(inputdata[, -1]))
    )
    # Calculate the means for each variable per group
    for (ii in 1:length(new_mats)) {
        submeans[ii] <- list(apply(new_mats[[ii]][, -1], 2, mean, na.rm=TRUE))
        means[ii, ] <- submeans[[ii]]
    }
    # Calculate the std for each variable per group
    for (ii in 1:length(new_mats)) {
        substd[ii] <- list(apply(new_mats[[ii]][, -1], 2, sd, na.rm=TRUE))
        std[ii, ] <- substd[[ii]]
    }

    cv <- matrix(
        nrow=length(means[, 1]),
        ncol=length(means[1, ]),
        dimnames=list(rownames(means), colnames(means))
    )
    for (ii in 1:length(means[, 1])) {
        for (jj in 1:length(means[1, ])) {
            cv[ii, jj] <- (std[ii, jj] / means[ii, jj])
        }
    }

    return(list(cv=cv, means=means, std=std))
}
