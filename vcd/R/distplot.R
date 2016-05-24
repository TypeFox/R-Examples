# added lwd arg, changed default point sizes

distplot <-
    function(x, type = c("poisson", "binomial", "nbinomial"),
             size = NULL, lambda = NULL, legend = TRUE, xlim = NULL, ylim = NULL,
             conf_int = TRUE, conf_level = 0.95, main = NULL,
             xlab = "Number of occurrences", ylab = "Distribution metameter",
             gp = gpar(cex = 0.8), lwd=2, name = "distplot", newpage = TRUE,
             pop = TRUE, return_grob = FALSE, ...)
{
    if(is.vector(x)) {
        x <- table(x)
    }
    if(is.table(x)) {
        if(length(dim(x)) > 1) stop ("x must be a 1-way table")
        freq <- as.vector(x)
        count <- as.numeric(names(x))
    } else {
        if(!(!is.null(ncol(x)) && ncol(x) == 2))
            stop("x must be a 2-column matrix or data.frame")
        freq <- as.vector(x[,1])
        count <- as.vector(x[,2])
    }

    myindex <- (1:length(freq))[freq > 0]
    mycount <- count[myindex]
    myfreq <- freq[myindex]

    switch(match.arg(type),

           "poisson" = {
               par.ml <- suppressWarnings(goodfit(x, type = type)$par$lambda)

               phi <- function(nk, k, N, size = NULL)
                   ifelse(nk > 0, lgamma(k + 1) + log(nk/N), NA)
               y <- phi(myfreq, mycount, sum(freq))
               if(!is.null(lambda)) y <- y + lambda - mycount * log(lambda)
               fm <- lm(y ~ mycount)
               par.estim <- exp(coef(fm)[2])
               names(par.estim) <- "lambda"
               txt <- "exp(slope)"
               if(!is.null(lambda)) {
                   par.estim <- par.estim * lambda
                   txt <- paste(txt, "x lambda")
               }
               legend.text <- paste(txt, "=", round(par.estim, digits = 3))
               if(is.null(main)) main <- "Poissoness plot"
  },

           "binomial" = {
               if(is.null(size)) {
                   size <- max(count)
                   warning("size was not given, taken as maximum count")
               }
               par.ml <- suppressWarnings(goodfit(x, type = type, par = list(size = size))$par$prob)

               phi <- function(nk, k, N, size = NULL)
                   log(nk) - log(N * choose(size, k))
               y <- phi(myfreq, mycount, sum(freq), size = size)
               fm <- lm(y ~ mycount)
               par.estim <- exp(coef(fm)[2])
               par.estim <- par.estim / (1 + par.estim)
               names(par.estim) <- "prob"
               legend.text <- paste("inv.logit(slope) =", round(par.estim, digits = 3))
               if(is.null(main)) main <- "Binomialness plot"
           },

           "nbinomial" = {
               if(is.null(size)) {
                   par.ml <- suppressWarnings(goodfit(x, type = type)$par)
                   size <- par.ml$size
                   par.ml <- par.ml$prob
               }else{
                   xbar <- weighted.mean(mycount, myfreq)
                   par.ml <- size / (size+xbar)
               }
               phi <- function(nk, k, N, size = NULL)
                   log(nk) - log(N * choose(size + k - 1, k))
               y <- phi(myfreq, mycount, sum(freq), size = size)
               fm <- lm(y ~ mycount)
               par.estim <- 1 - exp(coef(fm)[2])
               names(par.estim) <- "prob"
               legend.text <- paste("1-exp(slope) =", round(par.estim, digits = 3))
               if(is.null(main)) main <- "Negative binomialness plot"
           })

    yhat <- ifelse(myfreq > 1.5, myfreq - 0.67, 1/exp(1))
    yhat <- phi(yhat, mycount, sum(freq), size = size)
    if(!is.null(lambda)) yhat <- yhat + lambda - mycount * log(lambda)

    phat <- myfreq / sum(myfreq)
    ci.width <- qnorm(1-(1 - conf_level)/2) *
        sqrt(1-phat)/sqrt(myfreq - (0.25 * phat + 0.47)*sqrt(myfreq))

    RVAL <- cbind(count, freq, NA, NA, NA, NA, NA)
    RVAL[myindex,3:7] <- cbind(y,yhat,ci.width, yhat-ci.width, yhat + ci.width)
    RVAL <- as.data.frame(RVAL)
    names(RVAL) <- c("Counts", "Freq", "Metameter", "CI.center",
                     "CI.width", "CI.lower", "CI.upper")

    if(is.null(xlim)) xlim <- range(RVAL[,1])
    if(is.null(ylim)) ylim <- range(RVAL[,c(3,6,7)], na.rm = TRUE)
    xlim <- xlim + c(-1, 1) * diff(xlim) * 0.04
    ylim <- ylim + c(-1, 1) * diff(ylim) * 0.04

    if(newpage) grid.newpage()
    pushViewport(plotViewport(xscale = xlim, yscale = ylim, default.units = "native", name = name))
    grid.points(x = RVAL[,1], y = RVAL[,3], default.units = "native", gp = gp, ...)
    grid.lines(x = xlim, y = predict(fm, newdata = data.frame(mycount = xlim)),
               default.units = "native", gp = gpar(lwd=lwd, col = 2))
    grid.rect(gp = gpar(fill = "transparent"))
    grid.xaxis()
    grid.yaxis()
    grid.text(xlab, y = unit(-3.5, "lines"))
    grid.text(ylab, x = unit(-3, "lines"), rot = 90)
    grid.text(main, y = unit(1, "npc") + unit(2, "lines"), gp = gpar(fontface = "bold"))

    if(conf_int) {
        grid.points(x = RVAL[,1], y = RVAL[,4], pch = 19, gp = gpar(cex = 0.5))
        grid.segments(RVAL[,1], RVAL[,6], RVAL[,1], RVAL[,7], default.units = "native", gp = gpar(lty = 3))
    }

    if(legend) {
        mymin <- which.min(RVAL[,5])
        leg.x <- RVAL[mymin,1]
        if(RVAL[mymin,6] - ylim[1] > ylim[2] - RVAL[mymin,7])
            leg.y <- ylim[1] + 0.7 * (RVAL[mymin,6] - ylim[1])
        else leg.y <- ylim[2]

        legend.text <- c(paste("slope =", round(coef(fm)[2], digits = 3)),
                         paste("intercept =", round(coef(fm)[1], digits = 3)),
                         "", paste(names(par.estim),": ML =", round(par.ml, digits=3)),
                         legend.text)
        legend.text <- paste(legend.text, collapse = "\n")
        grid.text(legend.text, leg.x, leg.y - 0.05 * abs(leg.y), default.units = "native", just = c("left", "top"))
    }

    if(pop) popViewport() else upViewport()

    if (return_grob)
        structure(invisible(RVAL),
                  grob = grid.grab())
    else
        invisible(RVAL)
}
