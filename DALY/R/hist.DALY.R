## 'hist' S3 Method for class 'DALY'
## based on 'grid' package

hist.DALY <-
function(x, xval = c("DALY", "YLD", "YLL", "cases", "deaths"),
         prob = 0.95, central = c("mean", "median"),
         breaks = 25, fill = "grey95", ...){

  ## Check 'xval'
  xval <- match.arg(xval)

  ## Check 'prob'
  if (is.null(prob) || is.na(prob) || !is.numeric(prob))
    stop("'prob' must be a numeric value between 0 and 1")
  if (prob < 0 | prob > 1)
    stop("'prob' must be a numeric value between 0 and 1")

  ## Check 'central'
  central <- match.arg(central)

  ## Obtain aggregated 'xval'
  y <- aggregate(x, by = "total")[xval][[1]]

  ## Plot settings
  h <- hist(y, breaks = breaks, plot = FALSE)
  h2 <- hist(h$counts, plot = FALSE)

  x_val <- (h$breaks - min(h$breaks)) / (max(h$breaks) - min(h$breaks))
  y_val <- (h$counts / max(h$counts)) * .95
  y_val <- c(t(cbind(y_val, y_val, 0)))

  x_lb <- seq(min(h$breaks), max(h$breaks), length.out = 5)
  x_at <- (x_lb - min(h$breaks)) / (max(h$breaks) - min(h$breaks))

  y_lb <- round(seq(0, max(h2$breaks), length.out = 5), 0)
  y_at <- y_lb / max(h$counts)
  if (max(y_at) > 1.05){
    y_lb <- y_lb[-5]
    y_at <- y_at[-5]
  }

  ## Obtain summary statistics of 'xval'
  s <- summarize(y, .prob = prob)
  q <- approx(h$breaks, seq(0, 1, length.out = length(h$breaks)), s)$y
  at <- c(q[3] / 2, mean(q[3:4]), q[4] + (1 - q[4]) / 2)

  ## Plot histogram
  grid.newpage()

  ## top strip
  pushViewport(viewport(y = unit(1, "npc") - unit(3.5, "lines"),
                        x = unit(0, "npc") + unit(4, "lines"),
                        width = unit(.95, "npc") - unit(4, "lines"),
                        height = unit(1, "lines"),
                        just = c("left", "top"), name = "strip"))
  grid.rect(gp = gpar(col = "black", fill = "grey95"))
  grid.rect(x = unit(q[3], "npc"), width = unit(q[4] - q[3], "npc"),
            just = "left", gp = gpar(fill = "red"))
  p_lwr <- paste(100 * (1 - prob) / 2, "%", sep = "")
  p_mid <- paste(100 * prob, "%", sep = "")
  p_upr <- paste(100 * (1 - (1 - prob) / 2), "%", sep = "")
  grid.text(p_lwr, x = unit(at[1], "npc"), y = unit(.5, "npc"),
            gp = gpar(fontsize = 7))
  grid.text(p_mid, x = unit(at[2], "npc"), y = unit(.5, "npc"),
            gp = gpar(fontsize = 7))
  grid.text(p_upr, x = unit(at[3], "npc"), y = unit(.5, "npc"),
            gp = gpar(fontsize = 7))
  grid.points(pch = 25, size = unit(.5, "lines"), gp = gpar(fill = "black"),
              x = unit(q[3], "npc"), y = unit(1, "npc"))
  grid.points(pch = 25, size = unit(.5, "lines"), gp = gpar(fill = "black"),
              x = unit(q[4], "npc"), y = unit(1, "npc"))
  grid.points(pch = 25, size = unit(.5, "lines"), gp = gpar(fill = "red"),
              x = unit(ifelse(central == "mean", q[1], q[2]), "npc"),
              y = unit(1, "npc"))
  grid.text(round(ifelse(central == "mean", s[1], s[2]), 0),
            gp = gpar(fontsize = 7),
            x = unit(ifelse(central == "mean", q[1], q[2]), "npc"),
            y = unit(1.5, "npc"))
  grid.text(round(s[3], 0), gp = gpar(fontsize = 7),
            x = unit(q[3], "npc"), y = unit(1.5, "npc"))
  grid.text(round(s[4], 0), gp = gpar(fontsize = 7),
            x = unit(q[4], "npc"), y = unit(1.5, "npc"))

  ## histogram
  upViewport(1)
  pushViewport(viewport(y = unit(1, "npc") - unit(4.5, "lines"),
                        x = unit(0, "npc") + unit(4, "lines"),
                        width = unit(.95, "npc") - unit(4, "lines"),
                        height = unit(1, "npc") - unit(8.5, "lines"),
                        just = c("left", "top"), name = "hist"))
  grid.rect(gp = gpar(col = "black"))
  grid.xaxis(at = x_at, label = x_lb,
             gp = gpar(fontsize = 8))
  grid.yaxis(at = .95 * y_at, label = y_lb,
             gp = gpar(fontsize = 8))
  grid.text("Frequency", rot = 90,
            y = unit(0.5, "npc"), x = unit(-3, "lines"))
  grid.text(x$name, y = unit(1, "npc") + unit(2.5, "lines"),
            gp = gpar(fontsize = 15))
  grid.text(xval, y = unit(0, "npc") - unit(3, "lines"))
  grid.polygon(x = rep(x_val, each = 3)[-1],
               y = c(0, y_val, 0),
               gp = gpar(fill = fill, ...))
  grid.lines(x = unit(q[3], "npc"), y = unit(c(0, 1), "npc"))
  grid.lines(x = unit(q[4], "npc"), y = unit(c(0, 1), "npc"))
}