parallelplot <-
function (x, name = "", comb = -1, class = 0, obs = rep(0, 0), 
    col = 2, lty = 1, ...) 
{
    classes = as.numeric(factor(x[, ncol(x)]))
    len.class = table(classes)
    numclass = length(len.class)
#    print(numclass)
    x = x[, -ncol(x)]
    r = nrow(x)
    c = ncol(x)
    numgraphs = combinations(c)
    class.list = as.integer(rownames(len.class))
    xrnms = rownames(x)
    x <- apply(x, 2, function(x) (x - min(x))/(max(x) - min(x)))
    rownames(x) = xrnms
    if (class != 0) {
        same = (classes == class)
        x = x[same, ]
        col = class + 1
    }
    graphtitle = paste("Parallel Coordinate Plot for ", name)
    if (class == 0) 
        col = classes + 1
    if (comb == 0) {
        j = 1
        def.par <- par(mfrow = c(2, 2), font.lab = 2, font.sub = 2, 
            cex = 0.75, las = 2)
        for (k in 1:ncol(numgraphs)) {
            if (k%%4 == 1) {
                par(mfrow = c(2, 2))
            }
            varorder = numgraphs[, j]
            subtitle = paste("Combination #", j)
            matplot(1:c, t(x[, varorder]), type = "l", col = col, 
                lty = lty, xlab = "", ylab = "", main = graphtitle, 
                sub = subtitle, axes = FALSE, ...)
            axis(1, at = 1:c, labels = colnames(x[, varorder]))
            for (i in 1:c) lines(c(i, i), c(0, 1), col = "grey70")
            j = j + 1
        }
        par(def.par)
    }
    else {
        if (comb == -1) 
            varorder = 1:c
        else varorder = numgraphs[, comb]
        def.par <- par(font.lab = 2, font.sub = 2, cex = 0.75, 
            las = 2, bg = gray(0.8))
        subtitle = ("Original Attribute Order")
        matplot(1:c, t(x[, varorder]), type = "l", col = col, 
            lty = lty, xlab = "", ylab = "", main = graphtitle, 
            sub = subtitle, axes = FALSE, ...)
        axis(1, at = 1:c, labels = colnames(x[, varorder]))
        for (i in 1:c) lines(c(i, i), c(0, 1), col = "grey70")
        if (length(obs) != 0) {
            obsers = rep(0, 0)
            for (i in 1:length(obs)) obsers = c(obsers, which(rownames(x) == 
                obs[i]))
            colors = palette()[(numclass + 2):8]
            if (length(obsers) == 1) 
                matlines(1:c, x[obsers, varorder], lty = 1, 
                  lwd = 3, col = colors)
            else matlines(1:c, t(x[obsers, varorder]), lty = 1, 
                lwd = 3, col = colors)
            par(cex = 0.75)
            text(1, x[obsers, varorder[1]], rownames(x[obsers, 
                ]), pos = 3)
            text(c, x[obsers, varorder[c]], rownames(x[obsers, 
                ]), pos = 3)
            palette("default")
        }
    }
    invisible()
    par(def.par)
}
