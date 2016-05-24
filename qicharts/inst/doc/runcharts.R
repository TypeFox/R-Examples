## ----fig1, fig.cap='Figure 1'--------------------------------------------
library(qicharts)
set.seed(9)         # Lock random number generator
y <- rpois(24, 16)  # Random values from to plot
qic(y)              # Plot run chart of y

## ----fig2, fig.cap='Figure 2'--------------------------------------------
y[13:24] <- rpois(12, 24)  # Introduce a shift in process mean
qic(y)                     # Plot run chart of y

## ----limits table--------------------------------------------------------
n <- 10:30
data.frame(
  n.useful      = n,
  longest.run   = round(log2(n) + 3),
  min.crossings = qbinom(0.05, n - 1, 0.5))

## ----fig3, fig.cap='Figure 3'--------------------------------------------
qic(y, freeze = 12)

## ----fig4, fig.cap='Figure 4'--------------------------------------------
qic(y, breaks = 12)

## ----fig5, fig.cap='Figure 5'--------------------------------------------
y <- rbinom(24, 20, 0.5)                # Numerator
n <- sample(16:20, 24, replace = TRUE)  # Denominator
qic(y, n)                               # Plot run chart of y/n

## ----fig6, fig.cap='Figure 6'--------------------------------------------
startdate <- as.Date('2014-1-6')
date      <- seq.Date(startdate,         # Dates for x axis labels
                      by = 'day',
                      length.out = 24)
notes     <- NA
notes[18] <- 'This is a note'            # Character vector of annotations
qic(y, n,
    x     = date,
    main  = 'Run Chart', 
    ylab  = 'Proportion',
    xlab  = 'Date',
    notes = notes)

## ----data frame----------------------------------------------------------
date      <- seq.Date(startdate, by = 'day',       # 20 week long day sequence
                      length.out = 7 * 20)
n         <- sample(3:5, 7 * 20, replace = TRUE)   # Denominator vector
y         <- rbinom(7 * 20, n, 0.5)                # Numerator vector
week      <- as.Date(cut(date, 'week'))            # Subgrouping vector
d         <- data.frame(date, y, n, week)          # Data frame
head(d, 10)

## ----fig7, fig.cap='Figure 7'--------------------------------------------
qic(y, n, x = week, data = d)

