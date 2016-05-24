eofNum <-
function(x, n = nrow(x), scale. = TRUE) {

  # eigenvectors
  eigs <- prcomp(x, scale.=scale.)[["sdev"]]^2
  eigs.pct <- 100 * eigs/sum(eigs)

  # 0.95 confidence limits
  eigs.lo <- eigs * (1 - sqrt(2/n))
  eigs.hi <- eigs * (1 + sqrt(2/n))

  # cum. variance
  cumvar <- round(cumsum(eigs.pct), 1)

  # plot
  p <- ncol(x)
  d <- data.frame(rank = factor(1:p), eigs, eigs.lo, eigs.hi, cumvar)
  d <- within(d, cumvar.line <- eigs.hi + 0.02 * max(eigs.hi))
  d <- d[1:min(p, 10), ]
  ggplot(data = d, aes(x = rank, y = eigs)) +
    geom_errorbar(aes(x = rank, ymin = eigs.lo, ymax = eigs.hi),
                  width = 0.3) +
    geom_point(size = 3) +
    geom_text(aes(x = rank, y = cumvar.line, label = cumvar),
              size = 3, vjust = 0) +
    labs(list(x = "Rank", y = "Eigenvalue")) +
    theme(panel.grid.minor = element_blank())
}
