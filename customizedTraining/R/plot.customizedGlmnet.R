plot.customizedGlmnet <-
function(x, lambda, ...)
{
    groups = as.character(sort(unique(x$groupid)))
    nonzeroVariables = matrix(0, nrow = length(x$CTset) + 1,
        ncol = ncol(x$x$train))
    rownames(nonzeroVariables) = c("Standard", paste("Group", groups))

    selected = unique(unlist(predict(x$standard, s = lambda,
        type = 'nonzero')))
    nonzeroVariables['Standard', selected] = 1
    for (group in groups) {
        selected = NULL
        if (class(x$fit[[group]])[1] != "singleton") {
            selected = unique(unlist(predict(x$fit[[group]],
                s = lambda/x$fit[[group]]$nobs, type = "nonzero")))
        }
        nonzeroVariables[paste("Group", group), selected] = 1
    }
    
    selected = apply(t(nonzeroVariables), 1, rev)
    selected = t(selected[, colSums(selected) > 0])

    par(mar = c(5.1, 6.1, 4.1, 2.1))
    image(selected, col = c("white", "forestgreen"),
        axes = FALSE, xlab = "Variable", main = "Variables selected")
    axis(2, at = 0:length(groups)/length(groups),
        labels = colnames(selected), las = 1, lwd = 0)
}
