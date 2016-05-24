## ---- include=FALSE, echo=FALSE------------------------------------------
knitr::opts_chunk$set(fig.width = 7, fig.height = 5)
options(digits = 2)

## ---- echo=FALSE---------------------------------------------------------
library(cycleRtools)

## ---- fig.height = 10----------------------------------------------------
plot(x = intervaldata,     # "x" is the data, for consistency with other methods.
     y = 1:3,              # Which plots should be created? see below.
     xvar = "timer.min",   # What should be plotted on the x axis?
     xlab = "Time (min)",  # x axis label.
     laps = TRUE,          # Should different laps be coloured?
     breaks = TRUE)        # Should stoppages in the ride be shown?

## ------------------------------------------------------------------------
## Zoom to 0-50 minutes.
plot(intervaldata, y = 3, xvar = "timer.min", xlim = c(0, 50))

## ------------------------------------------------------------------------
zone_time(data = intervaldata,
          column = power.W,           # What are we interested in?
          zbounds = c(100, 200, 300), # Zone boundaries.
          pct = FALSE) / 60           # Output in minutes.

## How about time above and below CP, as a percentage?
## NB: column = power.W is the default.
zone_time(intervaldata, zbounds = 310, pct = TRUE) 

## ------------------------------------------------------------------------
zdist_plot(data = intervaldata,
           binwidth = 10,               # 10 Watt bins.
           zbounds = c(100, 200, 300),  # Zone boundaries.
           xlim = c(50, 400))           # Zoom to 50-400 Watts.

## ------------------------------------------------------------------------
summary(intervaldata)

## ------------------------------------------------------------------------
times_sec <- 2:20 * 60   # 2-20 minutes.
prof <- mmv(data = intervaldata, 
            column = power.W,    # Could also use speed.kmh.
            windows = times_sec)
print(prof)

## ------------------------------------------------------------------------
hypm <- lm(prof[1, ] ~ {1 / times_sec})  # Hyperbolic model.

## Critical Power (Watts) and W' (Joules) estimates
hypm <- setNames(coef(hypm), c("CP", "W'"))
print(hypm)

## Plot with the inverse model overlaid.
plot(times_sec, prof[1, ], ylim = c(hypm["CP"], max(prof[1, ])),
     xlab = "Time (sec)", ylab = "Power (Watts)")
curve((hypm["W'"] / x) + hypm["CP"], add = TRUE, col = "red")
abline(h = hypm["CP"], lty = 2)
legend("topright", legend = c("Model", "CP"), bty = "n",
       lty = c(1, 2), col = c("red", "black"))

## ------------------------------------------------------------------------
ms <- Pt_model(prof[1, ], times_sec)
print(ms)

plot(times_sec, prof[1, ], ylim = c(hypm["CP"], max(prof[1, ])),
     xlab = "Time (sec)", ylab = "Power (Watts)")
## Showing an exponential model, as it best fits these data.
curve(ms$Pfn$exp(x), add = TRUE, col = "red")

## ---- fig.height=3-------------------------------------------------------
library(leaflet)
leaflet(intervaldata) %>% addTiles() %>% addPolylines(~lon, ~lat)

