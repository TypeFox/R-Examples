## ---- echo = FALSE, message = FALSE--------------------------------------
library(knitr)
library(ggvis)
opts_chunk$set(comment = "#>", error = FALSE, tidy = FALSE)
opts_chunk$set(fig.width = 2.5, fig.height = 1.5, dpi = 100)

## ---- echo = FALSE, fig.width = 4----------------------------------------
# Histogram
faithful %>% ggvis(~eruptions, fill := "#ffffdd", fill.hover := "#eebbbb") %>%
  layer_histograms(width = 0.2) %>%
  add_axis("x", title = "eruptions") %>%
  add_axis("y", title = "count")

## ---- echo = FALSE, fig.width = 3, fig.height = 3------------------------
# Scatter plot with model fit line
mtcars %>%
  ggvis(~wt, ~mpg) %>%
  layer_points() %>%
  layer_smooths(span = input_slider(0.2, 1, 0.75, step = 0.05,
    label = "Smoothing span"))

