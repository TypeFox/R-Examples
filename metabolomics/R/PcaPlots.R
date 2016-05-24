
PcaPlots <- function(inputdata, y.axis=1, x.axis=2, center=TRUE, scale=TRUE, main=NULL,
    varplot=FALSE, multiplot=FALSE,n=5,cols=NULL, ...)
{
    # Get groups information
    group_list <- factor(inputdata[, 1], levels=unique(inputdata[, 1]))
    # Remove groups for data processing
    pca_data <- inputdata[, -1]

    write(' -> Performing PCA...', '')

    const_rows <- which(apply(pca_data, 2, var) == 0)
    if (length(const_rows) != 0) {
        pca_data <- pca_data[, -const_rows]
    }

    pca <- prcomp(pca_data, scale.=scale, center=center, ...)
    # Get the eigenvectors (loadings)
    eigenvecs <- pca$rotation

    # Get summary information
    summ <- summary(pca)
    importance <- summ$importance[2, ]

    if (varplot) {
        # Plot the explained variance
        barplot(summ$importance[2,c(1:10) ],
            col="#ee3333",                  # colour to plot bars
            main="Variance",               # plot title
            xlab="Principal Component",    # x-axis title
            ylab="Explained Variance",     # y-axis title
            cex.axis=1, cex.names=1,       # font size
            las=1                          # horizontal labels
        )
    }
    # Prepare a list of colours to use
    if (is.null(cols)) {
        col_list <- ColList(nlevels(group_list))
        cols_used <- col_list[as.numeric(group_list)]
    } else
        cols_used<-cols[as.numeric(group_list)]
    pch_list <- PchList(nlevels(group_list))
    pch_used <- pch_list[as.numeric(group_list)]

    # Store PCA data in a meaningful namespace
    pca_scores <- pca$x
    rownames(pca_scores) <- rownames(pca_data)

    if (multiplot) {
        pairs(pca_scores[,1:n], pch=pch_used, col=cols_used,
        labels=paste(
            "PC", c(1:n), "(", round(importance[c(1:n)] * 100, 2), "%)",
            sep="")
        )
    }

    # Plot PCA scores with sample names
    pca_mat <- cbind(pca_scores[, x.axis],pca_scores[, y.axis])
    x_percent <- sprintf("%.2f", importance[x.axis] * 100)
    y_percent <- sprintf("%.2f", importance[y.axis] * 100)
    pic_gen(pca_mat,
        plot_title=if (!is.null(main)) main else "PCA Score Plot\nSamples",
        x_label= paste("PC", x.axis, " (", x_percent, "%)", sep=""),
        y_label=paste("PC", y.axis, " (", y_percent, "%)", sep=""),
        cols_used=cols_used,
        pch_used=pch_used
    )
    # Plot PCA scores with group names
    pic_gen(pca_mat,
        plot_title=if (!is.null(main)) main else "PCA Score Plot\nGroups",
        plot_labels=group_list,
        x_label=paste("PC", x.axis, " (", x_percent, "%)", sep=""),
        y_label=paste("PC", y.axis," (", y_percent, "%)",  sep=""),
        cols_used=cols_used,
        pch_used=pch_used
    )
    eigen_mat <- cbind(eigenvecs[, x.axis],eigenvecs[, y.axis])
    # Loadings plot
    pic_gen(eigen_mat,
        "PCA Loading Plot",
        x_label=paste("PC", x.axis," (", x_percent, "%)", sep=""),
        y_label=paste("PC", y.axis," (", y_percent, "%)", sep=""),
        cols_used="black"
        #text_on=FALSE
    )

    write(' -> Done!', '')
}

