## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, 
  comment = '>',
  fig.align = 'center',
  fig.show = 'hold'
)

## ---- echo = FALSE-------------------------------------------------------
library(linbin)

## ---- fig.width = 5, fig.height = 5--------------------------------------
# Load event data
e <- elwha
e.filled <- fill_event_gaps(e, max.length = 1) # fill small gaps for the variable length bins (d)

# Design bins using different strategies
bins.a <- event_overlaps(e)[1:2]
bins.b <- seq_events(event_range(e), length.out = 33, adaptive = FALSE)
bins.c <- seq_events(event_coverage(e), length.out = 20, adaptive = FALSE)
bins.d <- seq_events(event_coverage(e.filled), length.out = 20, adaptive = TRUE)
bins <- rbind(cbind(bins.a, g = 1), cbind(bins.b, g = 2), cbind(bins.c, g = 3), cbind(bins.d, g = 4))

# Sample events at bins
e.bins <- sample_events(e, bins, list(weighted.mean, "mean.width", "unit.length"), 
                        scaled.cols = "unit.length")

# Plot binned data
plot_events(e.bins, group.col = "g", data.cols = "mean.width", col = "grey", border = "#666666", 
            main = c("(a) Flattened original data", "(b) Equal length bins", 
                     "(c) Equal coverage bins", "(d) Variable length bins"),
            xlabs = "Distance upstream (km)", ylabs = "Wetted width (m)",
            dim = c(4, 1), ylim = c(0, 56), xpd = NA)

## ---- fig.width = 6, fig.height = 6--------------------------------------
# Load event data
e <- quinault

# Design bins
bin.lengths <- c(100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600) # m
bins <- seq_events(event_range(e), by = bin.lengths / 1000) # km

# Sample events at bins
e.bins <- sample_events(e, bins, list(sum, "ONXX.*"), scaled.cols = "ONXX.*")

# Plot binned data
plot_events(e.bins, group.col = "group", data.cols = "ONXX.total", 
            main = paste0("Bin length = ", prettyNum(bin.lengths, ","), " m"), 
            xlabs = "Distance upstream (km)", ylabs = "Trout abundance",
            dim = c(3, 3), byrow = TRUE, oma = c(3, 3, 2, 2))

## ---- fig.width = 6, fig.height = 6--------------------------------------
plot_events(e.bins, group.col = "group", data.cols = "ONXX.[0-9]+", 
            main = paste0("Bin length = ", prettyNum(bin.lengths, ","), " m"), 
            xlabs = "Distance upstream (km)", ylabs = "Trout abundance",
            dim = c(3, 3), byrow = TRUE, oma = c(3, 3, 2, 2), col = heat.colors(3), border = NA)

## ---- fig.width = 5, fig.height = 5--------------------------------------
# Load NetMap data
d <- netmap

# Convert to event table
# (compute from and to endpoints from Netmap variables)
# (OUT_DIST = distance from outlet in km, LENGTH_M = length of unit in m)
d$from <- d$OUT_DIST
d$to <- d$from + (d$LENGTH_M / 1000)

# Seperate into mainstem and network
e.main <- d[d$CHAN_ID == 1, ]
e.net <- d

# Design bins
bins = seq_events(event_range(e.net), length.out = 10)

# Sample events at bins
fields = c("IP_CHINOOK", "IP_COHO", "IP_STEELHD", "BeavHab", "DEPTH_M")
e.bins.main = sample_events(e.main, bins, list(weighted.mean, fields, "LENGTH_M"), 
                            scaled.cols = "LENGTH_M")
e.bins.net = sample_events(e.net, bins, list(weighted.mean, fields, "LENGTH_M"), 
                           scaled.cols = "LENGTH_M")
e.bins = rbind(cbind(e.bins.main, group = 1), cbind(e.bins.net, group = 2))

# Plot binned data
plot_events(e.bins, group.col = 'group', data.cols = fields, sigfigs = c(3, 2),
            xlabs = c('Distance upstream (km)\nmainstem', 'Distance upstream (km)\nnetwork'),
            ylabs = c('Depth (m)', 'Proportion', 'IP', 'IP', 'IP'),
            oma = c(4, 3, 2, 2), mar = c(2, 4, 1.5, 0.5))

## ---- fig.width = 6, fig.height = 4--------------------------------------
# Load event data
d = fishmotion
e.motion = d[[1]]
e.origin = d[[2]]

# Design hourly bins
# (endpoints are in seconds since 1970-01-01 UTC)
bins = seq_events(event_range(e.motion), by = 3600)

# Sample events at bins
e.motion.bins = sample_events(e.motion, bins, list(length, 'region', by = 'region'))
e.origin.bins = sample_events(e.origin, bins, list(length, 'region', by = 'region'))

# Normalize by total fish present tagged in region 1
e.motion.bins$fish.1.norm = e.motion.bins$region.1 / e.origin.bins$region.1

# Prepare weekly data labels
bins[c("from.date", "to.date")] = lapply(bins[c("from", "to")], 
                                         as.POSIXct, origin = '1970-01-01', tz = "US/Alaska")
week.ticks = seq(trunc(min(bins$from.date), "day"), trunc(max(bins$from.date), "day"), by = "week")
week.labels = format(week.ticks, '%b-%d')

# Plot binned data
plot_events(e.motion.bins, data.cols = "fish.1.norm", yticks = c(0, 1, 2), 
            col = par("fg"), ylim = c(0, 2), plot.grid = TRUE, xpd = FALSE, 
            main = NA, xlabs = "Date (2008)", ylabs = "Relative abundance",
            xticks = week.ticks, xtick.labels = week.labels, oma = c(3, 2, 1, 2))

# Add daily vertical lines
days = seq(trunc(min(bins$from.date), "day"), trunc(max(bins$from.date), "day"), by = "day")
abline(v = days, col = 'grey')

