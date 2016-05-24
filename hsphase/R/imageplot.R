# This is a modified version of function which was written by Chris Seidel, available via this link:
# http://www.phaget4.org/R/image_matrix.html

imageplot <- function(x, title = c(), rv = FALSE, ...)
{
    x[x == 1] <- 3
    x[x == 2] <- 4
    min <- min(x)
    max <- max(x)
    yLabels <- rownames(x)
    xLabels <- colnames(x)
    if (length(list(...)))
    {
        Lst <- list(...)
        if (!is.null(Lst$zlim))
        {
            min <- Lst$zlim[1]
            max <- Lst$zlim[2]
        }
        if (!is.null(Lst$yLabels))
        {
            yLabels <- c(Lst$yLabels)
        }
        if (!is.null(Lst$xLabels))
        {
            xLabels <- c(Lst$xLabels)
        }
        if (!is.null(Lst$title))
        {
            title <- Lst$title
        }
    }
    if (is.null(xLabels))
    {
        xLabels <- c(1:ncol(x))
    }
    if (is.null(yLabels))
    {
        yLabels <- c(1:nrow(x))
    }
    if (rv == TRUE)
    {
        if (length(x[x == 0]) == 0)
        {
            ColorRamp <- c("#386CB0", "#FB8072")
        }
        else
        {
            ColorRamp <- c("#FFFFFF", "#FFFFFF", "#FFFFFF", "#386CB0", "#FB8072")
        }
    }
    else
    {
        if (length(x[x == 0]) == 0)
        {
            ColorRamp <- c("#FB8072", "#386CB0")
        }
        else
        {
            ColorRamp <- c("#FFFFFF", "#FFFFFF", "#FFFFFF", "#FB8072", "#386CB0")
        }
    }
    ColorLevels <- seq(min, max, length = length(ColorRamp))
    reverse <- nrow(x):1
    yLabels <- yLabels[reverse]
    x <- x[reverse, ]
    image(1:length(xLabels), 1:length(yLabels), t(x), col = ColorRamp, xlab = "", ylab = "", axes = FALSE, zlim = c(min, 
        max))
    if (!is.null(title))
    {
        title(main = title)
    }
    axis(BELOW <- 1, at = 1:length(xLabels), labels = xLabels, cex.axis = 0.7)
    axis(LEFT <- 2, at = 1:length(yLabels), labels = yLabels, las = HORIZONTAL <- 1, cex.axis = 0.7)
    mtext(side = 4, paste("Number of half-sibs: ", nrow(x)))
    mtext(side = 3, paste("Number of markers:", ncol(x)))
    layout(1)
} 
