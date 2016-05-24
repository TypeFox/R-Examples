FoldChange <- function(inputdata, paired=FALSE, plot.hist=TRUE, 
    saveoutput=FALSE, outputname="fc.results") 
{
    # Read in data, collect groups information
    inputdata <- editcolnames(inputdata)
    groups <- factor(inputdata[, 1], levels=unique(inputdata[, 1]))
    grp_levs <- levels(groups)
    
    if (length(grp_levs) > 2) {
        stop(
            paste("The number of groups is greater than 2. Use", 
                " LinearModelFit() instead."
            )
        )
    }
    # Prepare empty matrices and populate
    folds <- matrix(NA, nrow=length(grp_levs), 
        ncol=ncol(inputdata) - 1, 
        dimnames=list(grp_levs, colnames(inputdata)[2:ncol(inputdata)])
    )
    grp_len <- c()
    for (ii in 1:length(grp_levs)) {
        grp_len <- c(grp_len, length(which(groups == levels(groups)[ii])))
    }
    
    new_mats <- c()
    for (ii in 1:length(grp_levs)) {
        new_mats[ii] <- list(inputdata[which(groups == levels(groups)[ii]), ])
    }
    
    # Perform fold change calculations
    if (!paired) {
        submeans <- c()
        means <- matrix(nrow=length(grp_levs), 
            ncol=length(colnames(inputdata[, -1])), 
            dimnames=list(grp_levs, colnames(inputdata[, -1]))
        )
        for (ii in 1:length(new_mats)) {
            submeans[ii] <- list(
                apply(new_mats[[ii]][, -1], 2, mean, na.rm=TRUE)
            )
            means[ii, ] <- submeans[[ii]]
        }
        for (ii in 1:length(means[, 1])) {
            for (jj in 1:length(means[1, ])) {
                folds[ii, jj] <- means[ii, jj] - means[1, jj]
            }
        }
    } else {
        folds[1, ] <- 1
        folds[2, ] <- apply(
            (new_mats[[2]][, -1] - new_mats[[1]][, -1]), 2, mean, na.rm=TRUE
        )
    }
    
    # Plot and/or save
    if (plot.hist) {
        hist(folds[2, ], breaks=50, xlab="Fold change", main="Histogram")
    }
    if (saveoutput) {
        write.csv(t(folds), paste(c(outputname, ".csv"), collapse=""))
    }
    
    return(t(folds))
}
