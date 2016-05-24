sampdist <-
function (sampsize = c(3, 9, 30), seed = NULL, nsamp = 1000, FUN = mean,
            sampvals = function(n) exp(rnorm(n, mean = 0.5, sd = 0.3)),
            tck = NULL, plot.type = c("density", "qq"), layout = c(3,
                                                          1))
{
  if (!is.null(seed))
    set.seed(seed)
  ncases <- length(sampsize)
  y <- sampvals(nsamp)
  xlim = quantile(y, c(0.01, 0.99))
  xlim <- xlim + c(-1, 1) * diff(xlim) * 0.1
  samplingDist <- function(sampsize=3, nsamp=1000, FUN=mean)
    apply(matrix(sampvals(sampsize*nsamp), ncol=sampsize), 1, FUN)
  df <- data.frame(sapply(sampsize, function(x)samplingDist(x, nsamp=nsamp)))
  names(df) <- paste("y", sampsize, sep="")
  form <- formula(paste("~", paste(names(df), collapse="+")))
  lab <- lapply(sampsize, function(x) substitute(A, list(A = paste(x))))
  if (plot.type[1] == "density")
    gph <- lattice::densityplot(form, data=df, layout = layout, outer=TRUE,
                       plot.points = FALSE, panel = function(x, ...) {
                         lattice::panel.densityplot(x, ..., col = "black")
                         lattice::panel.densityplot(y, col = "gray40", lty = 2,
                                           ...)
                       }, xlim = xlim, xlab = "", scales = list(tck = tck),
                       between = list(x = 0.5), strip = lattice::strip.custom(strip.names = TRUE,
                                                  factor.levels = as.expression(lab), var.name = "Sample size",
                                                  sep = expression(" = ")))
  else if (plot.type[1] == "qq")
    gph <- lattice::qqmath(form, data = df, layout = layout, plot.points = FALSE,
                  outer=TRUE,
                  panel = function(x, ...) {
                    lattice::panel.qqmath(x, ..., col = "black", alpha=0.5)
                    lattice::panel.qqmath(y, col = "gray40", lty = 2, type = "l",
                                 ...)
                  }, xlab = "", xlim = c(-3, 3), ylab = "", scales = list(tck = tck),
                  between = list(x = 0.5), strip = lattice::strip.custom(strip.names = TRUE,
                                             factor.levels = as.expression(lab), var.name = "Sample size",
                                             sep = expression(" = ")))
  if (plot.type[1] %in% c("density", "qq"))
    print(gph)
  invisible(df)
}

