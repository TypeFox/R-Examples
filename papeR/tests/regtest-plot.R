#####
## Tests for plot interface

require("papeR")

## make sure other plot functions work
example("plot.data.frame", package = "graphics", ask = FALSE)
example("plot.formula", package = "graphics", ask = FALSE)

####
## check plot.ldf

data <- data.frame(a = 1:10, b = 10:1, c = rep(1:2, 5))
data$z <- as.factor(rep(2:3, each = 5))
data <- as.ldf(data)

## plot the data auto"magically"; numerics as boxplot, factors as barplots
par(mfrow = c(2,2))
plot(data)

## a single plot
plot(data, variables = "a")
## grouped plot
plot(data, variables = "a", by = "z")
## make "c" a factor and plot "c" vs. "z"
data$c <- as.factor(data$c)
plot(data, variables = "c", by = "z")
## the same
plot(data, variables = 3, by = 4)

## plot everithing against "b"
## (grouped boxplots, stacked barplots or scatterplots)
plot(data, with = "b")
