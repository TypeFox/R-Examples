seasonTrend <-
function(x, plot = FALSE, type = c("slope", "relative"), pval = .05, ...) {

  # validate args
  if (!is(x, "ts"))
    stop("x must be a 'ts'")
  type <- match.arg(type)

  # extend to full years
  first = start(x)[1]
  last = end(x)[1]
  fr <- frequency(x)
  x <- window(x, start = first, end = c(last, fr), extend = TRUE)

  # function for single vector
  st <- function(x) {
    x1 <- matrix(x, ncol = fr, byrow = TRUE)
    mannKen(x1)[, c(1:3, 6)]
  }

  # construct a data frame of all trends
  if (!is.matrix(x)) {
    ans <- data.frame(season = as.factor(1:fr), st(x), row.names = 1:fr)
  } else {
    nc <- ncol(x)
    colx <- colnames(x)
    series <- factor(rep(colx, each = fr), levels = colx, ordered = TRUE)
    season <- as.factor(rep(1:fr, times = nc))
    ans0 <- do.call(rbind, lapply(1:nc, function(i) st(x[, i])))
    ans <- data.frame(series, season, ans0, row.names = 1:nrow(ans0))
  }

  if (!plot) return(ans)

  ans[ans$miss >= 0.5, c("sen.slope", "sen.slope.rel")] <- NA
  ans$sig <- ifelse(ans$p.value < pval, TRUE, FALSE)
  v1 <- switch(type, slope = "sen.slope", relative = "sen.slope.rel")
  ylb <- switch(type,
                slope = expression(paste("Trend, units ", yr^{-1})),
                relative = expression(paste("Relative trend, ", yr^{-1})))
  names(ans)[match(v1, names(ans))] <- "trend"
  plt <- ggplot(ans, aes_string(x="season", y="trend", fill="sig")) +
           geom_bar(stat = "identity") +
           scale_fill_manual(name = "", values = c(`FALSE` = "grey65",
             `TRUE` = "dodgerblue"),
             labels = c(bquote(italic(p)>=.(pval)),bquote(italic(p)<.(pval)))) +
           labs(x = "Season", y = ylb) +
           theme(panel.grid.minor = element_blank())
  if (!is.matrix(x)) return(plt)
  plt + facet_wrap(~ series, ...)
}
