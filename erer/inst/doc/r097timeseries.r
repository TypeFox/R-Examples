# Creation and properties of a time series object
ts1 <- ts(data = 1:8, start = c(1990, 1), frequency = 12)
ts2 <- ts(data = matrix(data = 1:24, ncol = 4, byrow = FALSE),
  frequency = 12, start = c(2000, 1),
  names = paste("p", letters[1:4], sep = ""))
ts1; ts2
class(ts1); class(ts2)

start(ts2); end(ts2); frequency(ts2); deltat(ts2); cycle(ts2)
time(ts2); tsp(ts2); colnames(ts2); names(ts2); rownames(ts2)
colnames(ts2) <- paste("p", 1:4, sep = "")

# One time series object: less rows or columns
ma <- ts2[1:2, ]  # like a matrix in base R
mb <- window(x = ts2, start = c(2000, 3), end = c(2000, 5))
mc <- ts2[, 1:2]
md <- ts2[, c("p1", "p2")]

# Two time series objects: combining by row or column
ha <- ts.union(ts2, mb, dframe = FALSE)
hb <- cbind(ts2, mb); identical(ha, hb)
hc <- ts.intersect(ts2, mb, dframe = FALSE)
hd <- ts(data = rbind(ts2, ts2), start = c(1990, 1), frequency = 12)

# Convesion from time series to data frame
ff <- as.Date(paste(c(start(ts2), 1), collapse = "-"))
tt <- as.Date(paste(c(end(ts2), 1), collapse = "-"))
my.date <- seq(from = ff, to = tt, by = "month")
df <- data.frame(date = my.date, ts2)

# Lead, lag, and difference
ga <- lag(x = ts2, k = 1)   # a leading series
gb <- lag(x = ts2, k = -1)  # a lagged series
library(erer); gc <- bsLag(h = ts2, lag = 1)  # a lagged series with names
gd <- diff(x = ts2, lag = 1)                  # taking difference

# Head and tail
library(erer); head(ts2)                         # as a time series in erer
detach("package:erer", unload = TRUE); head(ts2) # as a matrix in base R