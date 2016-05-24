plotTsAnom <-
function(x, xlab = NULL, ylab = NULL,
         strip.labels = colnames(x), ...) {

  # Validate arguments
  if (!is.ts(x)) stop("x must be of class 'ts'")
  if (missing(xlab)) xlab = ""
  if (missing(ylab)) ylab = ""

  if (is.matrix(x)) {  # a matrix time series

    # fill possible spaces in column names so melt+merge will work
    strip.labels <- strip.labels
    colnames(x) <- gsub(' ', '.', colnames(x))

    # Create data frame
    x.mean = apply(x, 2, mean, na.rm=TRUE)
    x.mean.df <- data.frame(variable = factor(names(x.mean)), x.mean)
    d <- data.frame(time=as.Date(time(x)), x)
    d1 <- melt(d, id = 'time')
    d2 <- merge(d1, x.mean.df)
    d3 <- within(d2, variable <- factor(variable, levels = levels(variable),
                                        labels = strip.labels))
    d3 <- na.omit(d3)
    d3$ymin. <- with(d3, ifelse(value >= x.mean, x.mean, value))
    d3$ymax. <- with(d3, ifelse(value >= x.mean, value, x.mean))
    d3$colour. <- with(d3, value >= x.mean)
    # Plot
    ggplot(d3, aes_string(x="time", y="value", ymin="ymin.", ymax="ymax.",
                          colour="colour.")) +
      geom_linerange() +
      geom_hline(aes(yintercept = x.mean), size = 0.25) +
      labs(x = xlab, y = ylab) +
      facet_wrap(~ variable, ...) +
      theme(legend.position='none', panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=45, colour="grey50"))

  } else {  # a vector time series

    # Create data frame
    x.mean <- mean(x, na.rm = TRUE)
    d1 <- data.frame(time = as.Date(time(x)), x = as.numeric(x), x.mean)
    d1 <- na.omit(d1)
    d1$ymin. <- with(d1, ifelse(x >= x.mean, x.mean, x))
    d1$ymax. <- with(d1, ifelse(x >= x.mean, x, x.mean))
    d1$colour. <- with(d1, x >= x.mean)
    # Plot
    ggplot(d1, aes_string(x="time", y="x", ymin="ymin.", ymax="ymax.",
                          colour="colour.")) +
      geom_linerange() +
      geom_hline(aes(yintercept = x.mean), size = 0.25) +
      labs(x = xlab, y = ylab) +
      theme(legend.position='none', panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle=45, colour="grey50"))
  }
}
