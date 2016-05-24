plotSeason <-
function(x, type = c('by.era', 'by.month'), num.era = 4,
  same.plot = TRUE, ylab = NULL, num.col = 3) {

  # Validate args
  if (!is(x, 'ts') || is(x, 'mts'))
    stop("x must be a single 'ts'")
  type <- match.arg(type)

  # Turn time series into data.frame
  sx <- start(x)[1]
  ex <- end(x)[1]
  x <- window(x, start = sx, end = c(ex, 12), extend = TRUE)
  d <- data.frame(x = as.numeric(x), mon = ordered(month.abb[cycle(x)],
      levels = month.abb), yr = as.numeric(floor(time(x))))

  # Take care of case where num.era is a scalar
  if (length(num.era)==1) {
    if (num.era<1 || round(num.era)!=num.era) {
      stop("num.era must be a whole number > 0")
    } else {
      num.era <- round((0:num.era) * (ex-sx)/num.era + sx, 0)
    }
  }

  if (type == 'by.era') {
    # Break data into eras
    d$era <- cut(d$yr, breaks = num.era, include.lowest = TRUE, dig.lab = 4,
                 ordered_result = TRUE)
    colnames(d)[1] <- 'value'
    d <- na.omit(d)

    # Find missing fraction by month and era
    t0 <- table(d$mon, d$era)
    t1 <- sweep(t0, 2, diff(num.era), '/')
    t2 <- t1 < 0.5
    t3 <- melt(t2)
    colnames(t3) <- c('mon', 'era', 'too.few')
    t4 <- within(t3, {
      mon <- ordered(mon, levels = levels(d$mon))
      if (length(unique(era))>1)
        era <- ordered(era, levels = levels(d$era))
      }
    )
    d1 <- merge(d, t4)

    if (same.plot) {
       # Nest eras within months
       ggplot(d1, aes_string(x="mon", y="value", fill="era")) +
          geom_boxplot(size=.2, position='dodge') +
          labs(x="", y=ylab, fill="Era")
    } else {
       # Nest months within eras
       cols <- c(`TRUE` = "red", `FALSE` = "blue")
       p1 <- ggplot(d1, aes_string(x="mon", y="value", colour="too.few")) +
          geom_boxplot(size = .2) +
          scale_x_discrete('', breaks = month.abb,
                           labels = c('Jan', ' ', ' ', 'Apr', ' ', ' ', 'Jul',
                                      ' ', ' ', 'Oct', ' ', ' ')) +
          scale_y_continuous(ylab) +
          scale_colour_manual("", values=cols, guide="none") +
          theme(panel.grid.minor = element_blank(),
                axis.text.x = element_text(angle=45, colour="grey50"))
       if (length(num.era) > 2)
          p1 <- p1 + facet_wrap(~ era, nrow = 1)
       p1
    }

    } else {
      # Plot standardized anomalies for each month
      x1 <- ts2df(x)
      x2 <- ts(x1, start = start(x))
      plotTsAnom(x2, ylab = ylab, scales = "free_y")
    }
}
