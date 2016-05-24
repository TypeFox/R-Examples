plotTsTile <-
function(x, plot.title = NULL, legend.title = NULL, four = TRUE,
         loganom = TRUE, square = TRUE, legend = TRUE,
         trim = TRUE, overall = TRUE, stat = c("median", "mean")) {

  # Validate args
  if (!is(x, "ts") || is(x, "mts") || !identical(frequency(x), 12))
    stop("x must be a vector of class 'ts' with frequency = 12")
  stat <- match.arg(stat)

  # Define center function
  center <- function(x, type=stat) {
    switch(type,
           mean = mean(x, na.rm=TRUE),
           median = median(x, na.rm=TRUE)
    )
  }

  # trim leading and trailing NAS
  if (trim) {
    x <- as.zoo(x)
    x <- na.trim(x)
    x <- as.ts(x)
  }

  # Complete partial years by padding with NAs
  sx <- start(x)[1]
  ex <- end(x)[1]
  x1 <- window(x, start = c(sx, 1), end = c(ex, 12), extend = TRUE)

  # Transform to log-anomalies
  if (loganom) {
    if (any(x1 <= 0, na.rm = TRUE)) {
      stop("All values must be positive if loganom=TRUE")
    }
    else {
      if (overall) {
        x1 <- x1/center(x1)
      }
      else {
        x1 <- as.matrix(ts2df(x1))
        x1 <- sweep(x1, 2, apply(x1, 2, center), "/")
        x1 <- ts(as.vector(t(x1)), start = c(sx, 1), frequency=12)
      }
    }
    x1 <- log10(x1)
  }

  # Break data.
  if (four) {
    mmin <- min(x1, na.rm = TRUE)
    mlo <- center(x1[x1 < 0])
    mhi <- center(x1[x1 > 0])
    mmax <- max(x1, na.rm = TRUE)
    the.breaks <- c(mmin, mlo, 0, mhi, mmax)
  }
  else {
    the.breaks <- quantile(x1, probs = seq(0, 1, 0.1),
                           na.rm = TRUE)
  }
  len <- length(the.breaks)
  if (length(unique(the.breaks)) < len)
    stop("Breaks between groups are\nnot unique: insufficient unique data.")
  x2 <- cut(x1, breaks = the.breaks, include.lowest = TRUE,
            dig.lab = 2)
  x3 <- data.frame(yr = floor(time(x1)), mon = ordered(month.abb[cycle(x1)],
                   levels = month.abb), value = x2)

  # Plot it.
  mypalette <- colorRampPalette(c("darkblue", "lightblue",
                                  "pink", "red"))
  cols <- mypalette(len - 1)
  p1 <- ggplot(x3, aes_string(x="yr", y="mon", fill="value")) +
    geom_tile(colour = "white", size = 0.25) +
    scale_x_continuous(name = "", expand = c(0, 0)) +
    scale_y_discrete(name = "", expand = c(0, 0)) +
    scale_fill_manual(name = legend.title, values = cols,
                      breaks = levels(x3$value), labels = levels(x3$value)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    labs(title = plot.title)
  if (!legend)
    p1 <- p1 + theme(legend.position = "none")
  if (square)
    p1 + coord_equal()
  else p1
}
