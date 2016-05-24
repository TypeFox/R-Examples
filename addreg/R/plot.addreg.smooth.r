plot.addreg.smooth <- function(x, type = c("response","link"), at = data.frame(), 
                        knotlines = TRUE, nobs = 1000, ...)
{
  type <- match.arg(type)
  gp <- interpret.addreg.smooth(x$full.formula)
  allvars <- names(get_all_vars(delete.response(gp$terms), data = x$data))
  smthterms <- names(gp$smooth.spec)
  for (i in 1:length(smthterms)) {
    allvars2 <- allvars[allvars != smthterms[i]]
    smthtype <- class(gp$smooth.spec[[smthterms[i]]])
	ltype = "l"
	if (smthtype == "Iso.smooth")
        ltype = "s"
    smthknots <- x$knots[[smthterms[i]]]
    smthx <- seq(smthknots[1], smthknots[length(smthknots)], len = nobs)
    dat <- list()
    at.new <- at
    at.new[[smthterms[i]]] <- NULL
    if (!setequal(intersect(allvars2, names(at.new)), allvars2))
        stop(gettextf("error predicting for %s: 'at' must contain %s",
                smthterms[i], paste(allvars2, collapse = ", ")), domain = NA)
    if (nrow(at.new) == 0)           
        dat[[1]] <- setNames(data.frame(tmp = smthx), smthterms[i])
    else {
        for (j in 1:nrow(at.new))
		      dat[[j]] <- cbind(setNames(data.frame(tmp = smthx), smthterms[i]), 
                                setNames(data.frame(tmp = at.new[j,]), names(at.new)), row.names = NULL)
    }
    pred <- list()
    for (j in 1:length(dat))
        pred[[j]] <- predict(x, newdata = dat[[j]], type = type)
    ylim <- c(min(sapply(pred, min)), max(sapply(pred, max)))
    plot(smthx, pred[[1]], type = "n", ylim = ylim, xlab = smthterms[i], ylab = type)
		if(knotlines & smthtype != "Iso.smooth")
        abline(v = smthknots, col = "gray", lwd = 0.5)
    plotcols <- rainbow(length(pred))
    for (j in 1:length(pred))
        lines(smthx, pred[[j]], type = ltype, col = plotcols[j], ...)
  }
}