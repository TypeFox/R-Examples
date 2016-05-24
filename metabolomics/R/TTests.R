TTests <- function(inputdata, alternative="two.sided", paired=FALSE, 
    padjmethod="BH", saveoutput=FALSE, outputname="ttest.results", ...)
{
    # Get groups information
    groups <- factor(inputdata[, 1], levels=unique(inputdata[, 1]))
    # Remove groups for data processing
    ttest_data <- inputdata[, -1]
    ttest_data_t <- t(ttest_data)
    grp_levs <- levels(groups)
    
    if (length(grp_levs) > 2) {
        stop(
            paste("The number of groups is greater than 2.",
                "Use LinearModelFit() instead."
            )
        )
    }
    
    # Edit the column names if necessary
    colnames(ttest_data) <- if (
        length(
            grep("^X[\\d]", colnames(ttest_data), perl=TRUE)
        ) != 0
    ) {
        gsub("^X([\\d].*)", "\\1", colnames(ttest_data), perl=TRUE)
    } else {
        colnames(ttest_data)
    }
    row.names(ttest_data_t) <- colnames(ttest_data)
    
    #
    #    Prepare the supporting data for the test

    # Separate the groups
    x <- t(ttest_data_t[, which(groups==grp_levs[1])])
    y <- t(ttest_data_t[, which(groups==grp_levs[2])])
    # Create an empty matrix
    vals <- matrix(NA, nrow(ttest_data_t), ncol=3)
    
    #
    #    Perform the t-Test
    #
    for(ii in 1:nrow(ttest_data_t)) {
        if (min(length(na.omit(x[, ii])), length(na.omit(y[, ii]))) > 1) {
            if (var(c(x[, ii], y[, ii]), na.rm=TRUE) != 0) {
                test <- t.test(y[, ii], x[, ii], alternative=alternative,
                    paired=paired, ...
                )
                vals[ii, 1] <- test$statistic
                vals[ii, 2] <- test$p.value
            }
        }
    }
    vals[, 3] <- p.adjust(vals[, 2], method=padjmethod, 
        n=length(na.omit(vals[, 2]))
    )
    
    # Prepare row and column labels for output
    row.names(vals) <- row.names(ttest_data_t)
    colnames(vals) <- c("t-statistic", "p-value", "Adjusted p-value")
    
    if (saveoutput) {
        write.csv(vals, paste(c(outputname, ".csv"), collapse=""))
    }
    
    return(vals)
}
