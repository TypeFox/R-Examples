MetBoxPlots <- function(inputdata, metname, cols=NULL, main=NULL, 
    cex.main=NULL, ...)
{       ## comments!
    groups <- factor(inputdata[, 1], levels = unique(inputdata[, 1]))
    unique.groups <- levels(groups)
    if(is.null(cols))  {
        cols <- ColList(length(unique.groups))
        box_cols <- c(rep(NA, length(rownames(inputdata))))
        for (ii in 1:length(inputdata[, 1])) {
          box_cols[ii] <- cols[which(unique.groups == inputdata[, 1][ii])]
        }
    }    
    plot(groups, inputdata[, metname], col=cols, cex.main=cex.main, 
        main=main, ...
    )
}
