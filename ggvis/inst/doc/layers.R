## ---- echo = FALSE, message = FALSE--------------------------------------
library(knitr)
library(ggvis)
opts_chunk$set(comment = "#>", error = FALSE, tidy = FALSE)
opts_chunk$set(fig.width = 3.5, fig.height = 2.5, dpi = 100)

## ------------------------------------------------------------------------
mtcars %>% ggvis(x = ~wt, y = ~mpg, stroke := "red") %>% layer_points()
mtcars %>% ggvis() %>% layer_points(x = ~wt, y = ~mpg, stroke := "red")

## ------------------------------------------------------------------------
hec <- as.data.frame(xtabs(Freq ~ Hair + Eye, HairEyeColor))

hec %>% 
  ggvis(~Hair, ~Eye, fill = ~Freq) %>% 
  layer_rects(width = band(), height = band()) %>%
  scale_nominal("x", padding = 0, points = FALSE) %>%
  scale_nominal("y", padding = 0, points = FALSE)

## ------------------------------------------------------------------------
df <- data.frame(x = c(1, 1, 2, 2), y = c(2, 1, 1, 2))
df %>% ggvis(~x, ~y, stroke := "red") %>% layer_paths()
# Add a fill colour to make it a polygon
df %>% ggvis(~x, ~y, fill := "red") %>% layer_paths()

## ---- eval = FALSE-------------------------------------------------------
#  ggvis() %>%
#    layer_points(x = ~disp, y = ~wt, data = mtcars) %>%
#    layer_paths(x := 0, y = ~mean(mtcars$wt, x2 := prop_group())) %>%
#    layer_paths(x = ~mean(mtcars$disp), y := 0, y2 := prop_group())

## ------------------------------------------------------------------------
df <- data.frame(x = 1:10, y = (1:10) ^ 2)
df %>% ggvis(~x, ~y, y2 := 0) %>% layer_ribbons()

# Set height in pixels
df %>% ggvis(~x, ~y, height := 20) %>% layer_ribbons()

## ------------------------------------------------------------------------
df %>% ggvis(~x, ~y, prop("height", 80, scale = "y")) %>% layer_ribbons()

df <- data.frame(x = 1:10, y = (1:10) ^ 2)
df %>% ggvis(~x, ~y) %>%
  layer_ribbons(prop("height", input_slider(0, 100), scale = "y")) %>%
  layer_paths(stroke := "red", strokeWidth := 10)

## ------------------------------------------------------------------------
df %>% ggvis(~x, y = ~y - 2, y2 = ~y + 2) %>% layer_ribbons()

## ------------------------------------------------------------------------
mtcars %>% ggvis(~wt, ~mpg) %>%
  layer_points() %>%
  group_by(cyl) %>%
  layer_paths()

mtcars %>% ggvis(~wt, ~mpg) %>%
  layer_points() %>%
  auto_group() %>%
  layer_paths()

## ------------------------------------------------------------------------

mtcars %>% 
  dplyr::mutate(cyl2 = factor(cyl)) %>% 
  ggvis(~wt, ~mpg, stroke = ~cyl2) %>% 
  layer_lines()

## ------------------------------------------------------------------------
layer_smooths
layer_histograms

