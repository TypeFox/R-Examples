# across, within group
# across - remove med for met, boxplot sample
# within - all groups, remove median for each group, then plot
RlaPlots <- function(inputdata, type=c("ag", "wg"), cols=NULL,
    cex.axis=0.8, las=2, ylim=c(-2, 2), oma=c(7, 4, 4, 2) + 0.1, ...)
{
    type <- match.arg(type)
    groups <- factor(inputdata[, 1], levels = unique(inputdata[, 1]))
    unique.groups <- levels(groups)
    if (is.null(cols)) 
        cols <- ColList(length(unique.groups))
    box_cols <- c(rep(NA, length(rownames(inputdata))))
    for (ii in 1:length(inputdata[, 1])) 
        box_cols[ii] <- cols[which(unique.groups == inputdata[, 1][ii])]
    
    
    
    # Within groups
    if(type == "wg") {
        out_data<-data.frame()
        for (grp in unique.groups) {
            submat <- inputdata[which(inputdata[, 1] == grp), -1]
            med_vals <- apply(submat, 2, median)
            swept_mat <- sweep(submat, 2, med_vals, "-")
            out_data <- rbind(out_data, swept_mat)
        }
    # Across groups (i.e. type == "ag")
    } else  {
        med_vals <- apply(inputdata[, -1], 2, median)
        out_data <- sweep(inputdata[, -1], 2, med_vals, "-")
    }
        
    boxplot(t(out_data),
        cex.axis=cex.axis,                 # font size
        las=las,                           # label orientation
        col=box_cols,                      # colours
        ylim=ylim,                         # y-axis range
        oma=oma,                           # outer margin size
        ...
    )
    
    abline(h=0)
}
