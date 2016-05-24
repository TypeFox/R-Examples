"anovaPlot" <-
function (obj, stacked = TRUE, base = TRUE, axes = TRUE, faclab = TRUE, 
    labels = FALSE, cex = par("cex"), cex.lab = par("cex.lab"), 
    ...) 
{
    if (!any(class(obj) == "lm")) 
        stop(paste("Object", deparse(substitute(obj)), "should be of class 'lm'"))
    if (!any(class(obj) == "aov")) 
        obj <- do.call("aov", as.list(obj$call[-1]))
    tables <- model.tables(obj, type = "effects")[[1]]
    Residuals <- resid(obj)
    factors <- names(tables)
    k <- length(factors)
    Df <- anova(obj)[, "Df"]
    Scale.Factor <- sqrt(Df[length(Df)]/Df[-length(Df)])
    lst <- lst.dev <- lst.names <- list()
    for (i in 1:k) {
        nf <- length(unlist(strsplit(names(tables[i]), ":")))
        if (nf == 1) {
            label <- dimnames(tables[i][[1]])[[1]]
            eff <- as.numeric(tables[i][[1]])
            names(eff) <- label
        }
        else {
            label <- dimnames(tables[i][[1]])[[1]]
            eff <- as.numeric(tables[i][[1]])
            for (j in 2:nf) {
                lab <- dimnames(tables[i][[1]])[[j]]
                label <- paste(rep(label, length(lab)), rep(lab, 
                  each = length(label)), sep = ":")
            }
            names(eff) <- label
        }
        lst.dev[[factors[i]]] <- Scale.Factor[i] * eff
    }
    xmax <- max(abs(range(c(unlist(lst.dev), Residuals))))
    xlim <- c(-xmax, +xmax)
    xpd <- par("xpd")
    on.exit(par(xpd = xpd))
    par(xpd = TRUE)
    plot(c(0, 1), c(0, 1), xlim = xlim, type = "n", xlab = "", 
        ylab = "", frame = FALSE, axes = FALSE, ...)
    if (stacked) 
        h <- 1/(k + 2)
    else h <- 1/(k + 1)
    hinc <- h/(k + 1)
    for (i in 1:k) {
        y <- (k - i + 2) * h
        if (labels) 
            lab <- names(lst.dev[[i]])
        else lab <- NULL
        lst[[factors[i]]] <- dots(lst.dev[[i]], y = y + hinc, 
            xlim = xlim, hmax = y + hinc + h, stacked = stacked, 
            base = base, axes = FALSE, labels = lab, cex = cex)
        if (faclab) 
            text(0, y, labels = paste("scaled <", factors[i], 
                "> deviations"), adj = 0.5, cex = cex.lab)
    }
    lst[["Rediduals"]] <- dots(Residuals, y = 0, xlim = xlim, 
        hmax = 2 * h, stacked = stacked, base = FALSE, axes = FALSE, 
        labels = NULL)
    if (axes) 
        axis(1, at = pretty(xlim), mgp = c(1.5, 0.5, 0), line = -1 + 
            2 * hinc)
    if (faclab & axes) 
        mtext("Residuals", side = 1, cex = cex.lab, line = 1 + 
            hinc)
    if (faclab & !axes) 
        mtext("Residuals", side = 1, cex = cex.lab, line = 0)
    invisible(lst)
}
