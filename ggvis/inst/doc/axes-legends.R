## ---- echo = FALSE, message = FALSE--------------------------------------
library(knitr)
library(ggvis)
opts_chunk$set(comment = "#>", error = FALSE, tidy = FALSE)
opts_chunk$set(fig.width = 3.5, fig.height = 2.5, dpi = 100)

## ---- eval = FALSE-------------------------------------------------------
#  add_axis(vis, "x")
#  add_axis(vis, "y")
#  add_legend(vis, "stroke")
#  add_legend(vis, "size")
#  # Display multiple scales in one legend:
#  add_legend(vis, "stroke", "size")

## ---- eval = FALSE-------------------------------------------------------
#  add_axis(vis, "x", title = "My x variable")
#  add_legend(vis, "fill", title = "Some interesting colours")

## ---- results = 'asis'---------------------------------------------------
mtcars %>% ggvis(~wt, ~mpg) %>%
  layer_points() %>%
  add_axis("x", properties = axis_props(
    axis = list(stroke = "red", strokeWidth = 5),
    grid = list(stroke = "blue"),
    ticks = list(stroke = "blue", strokeWidth = 2),
    labels = list(angle = 45, align = "left", fontSize = 20)
  ))

## ---- results = 'asis'---------------------------------------------------
mtcars %>% ggvis(~wt, ~mpg) %>% layer_points()

mtcars %>% ggvis(~wt, ~mpg) %>% layer_points() %>%
  add_axis("x", title = "Weight") %>%
  add_axis("y", title = "Miles per gallon")

# Use title offset to push the titles further away
mtcars %>% ggvis(~wt, ~mpg) %>%
  layer_points() %>%
  add_axis("x", title = "Weight", title_offset = 50) %>%
  add_axis("y", title = "Miles per gallon", title_offset = 50)

## ---- results = 'asis'---------------------------------------------------
# Change ticks and subdivide with minor ticks
mtcars %>% ggvis(~wt, ~mpg) %>%
  layer_points() %>%
  add_axis("x", subdivide = 9, values = 1:6) %>%
  add_axis("y", subdivide = 1, values = seq(10, 34, by = 2))

# Make the minor ticks smaller and the end ticks longer
mtcars %>% ggvis(~wt, ~mpg) %>%
  layer_points() %>%
  add_axis("x", subdivide = 9, values = 1:6, tick_size_major = 10,
    tick_size_minor = 5, tick_size_end = 15, tick_padding = 20)

## ---- results = 'asis'---------------------------------------------------
mtcars %>% ggvis(~wt, ~mpg) %>%
  layer_points() %>%
  add_axis("x", orient = "top") %>%
  add_axis("y", orient = "right")

## ---- results = 'asis'---------------------------------------------------
mtcars %>% ggvis(~wt, ~mpg) %>%
  layer_points() %>%
  add_axis("x", orient = "bottom") %>%
  add_axis("x", orient = "top")

## ---- results = 'asis'---------------------------------------------------
mtcars %>% ggvis(~wt, ~mpg) %>%
  layer_points() %>%
  add_axis("x") %>%
  add_axis("x", offset = 40, grid = FALSE)

## ---- results = 'asis'---------------------------------------------------
mtcars2 <- mtcars %>% dplyr::mutate(cyl = ordered(mtcars$cyl))
mtcars2 %>% ggvis(~mpg, ~wt, size = ~cyl, fill = ~cyl) %>% layer_points()
mtcars2 %>% ggvis(~mpg, ~wt, size = ~cyl, fill = ~cyl) %>% layer_points() %>%
  add_legend(c("size", "fill"))

