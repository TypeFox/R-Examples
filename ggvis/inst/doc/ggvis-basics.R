## ---- echo = FALSE, message = FALSE--------------------------------------
library(knitr)
library(ggvis)
opts_chunk$set(comment = "#>", error = FALSE, tidy = FALSE)
opts_chunk$set(fig.width = 3.5, fig.height = 2.5, dpi = 100)

## ------------------------------------------------------------------------
p <- ggvis(mtcars, x = ~wt, y = ~mpg)

## ------------------------------------------------------------------------
layer_points(p)

## ------------------------------------------------------------------------
layer_points(ggvis(mtcars, x = ~wt, y = ~mpg))

## ------------------------------------------------------------------------
mtcars %>%
  ggvis(x = ~wt, y = ~mpg) %>%
  layer_points()

## ---- message = FALSE----------------------------------------------------
library(dplyr)
mtcars %>%
  ggvis(x = ~mpg, y = ~disp) %>%
  mutate(disp = disp / 61.0237) %>% # convert engine displacment to litres
  layer_points()

## ------------------------------------------------------------------------
mtcars %>%
  ggvis(~mpg, ~disp) %>%
  layer_points()

## ------------------------------------------------------------------------
mtcars %>% ggvis(~mpg, ~disp, stroke = ~vs) %>% layer_points()
mtcars %>% ggvis(~mpg, ~disp, fill = ~vs) %>% layer_points()
mtcars %>% ggvis(~mpg, ~disp, size = ~vs) %>% layer_points()
mtcars %>% ggvis(~mpg, ~disp, shape = ~factor(cyl)) %>% layer_points()

## ------------------------------------------------------------------------
mtcars %>% ggvis(~wt, ~mpg, fill := "red", stroke := "black") %>% layer_points()
mtcars %>% ggvis(~wt, ~mpg, size := 300, opacity := 0.4) %>% layer_points()
mtcars %>% ggvis(~wt, ~mpg, shape := "cross") %>% layer_points()

## ------------------------------------------------------------------------
mtcars %>% 
  ggvis(~wt, ~mpg, 
    size := input_slider(10, 100),
    opacity := input_slider(0, 1)
  ) %>% 
  layer_points()

## ------------------------------------------------------------------------
mtcars %>% 
  ggvis(~wt) %>% 
  layer_histograms(width =  input_slider(0, 2, step = 0.10, label = "width"),
                   center = input_slider(0, 2, step = 0.05, label = "center"))

## ------------------------------------------------------------------------
keys_s <- left_right(10, 1000, step = 50)
mtcars %>% ggvis(~wt, ~mpg, size := keys_s, opacity := 0.5) %>% layer_points()

## ------------------------------------------------------------------------
mtcars %>% ggvis(~wt, ~mpg) %>% 
  layer_points() %>% 
  add_tooltip(function(df) df$wt)

## ------------------------------------------------------------------------
mtcars %>% ggvis(~wt, ~mpg) %>% layer_points()

## ------------------------------------------------------------------------
df <- data.frame(x = 1:10, y = runif(10))
df %>% ggvis(~x, ~y) %>% layer_paths()

## ------------------------------------------------------------------------
t <- seq(0, 2 * pi, length = 100)
df <- data.frame(x = sin(t), y = cos(t))
df %>% ggvis(~x, ~y) %>% layer_paths(fill := "red")

## ------------------------------------------------------------------------
df <- data.frame(x = 1:10, y = runif(10))
df %>% ggvis(~x, ~y) %>% layer_ribbons()
df %>% ggvis(~x, ~y + 0.1, y2 = ~y - 0.1) %>% layer_ribbons()

## ------------------------------------------------------------------------
set.seed(1014)
df <- data.frame(x1 = runif(5), x2 = runif(5), y1 = runif(5), y2 = runif(5))
df %>% ggvis(~x1, ~y1, x2 = ~x2, y2 = ~y2, fillOpacity := 0.1) %>% layer_rects()

## ------------------------------------------------------------------------
df <- data.frame(x = 3:1, y = c(1, 3, 2), label = c("a", "b", "c"))
df %>% ggvis(~x, ~y, text := ~label) %>% layer_text()
df %>% ggvis(~x, ~y, text := ~label) %>% layer_text(fontSize := 50)
df %>% ggvis(~x, ~y, text := ~label) %>% layer_text(angle := 45)

## ------------------------------------------------------------------------
t <- seq(0, 2 * pi, length = 20)
df <- data.frame(x = sin(t), y = cos(t))
df %>% ggvis(~x, ~y) %>% layer_paths()
df %>% ggvis(~x, ~y) %>% layer_lines()

## ------------------------------------------------------------------------
df %>% ggvis(~x, ~y) %>% arrange(x) %>% layer_paths()

## ------------------------------------------------------------------------
mtcars %>% ggvis(~mpg) %>% layer_histograms()

# Or equivalently
binned <- mtcars %>% compute_bin(~mpg) 
binned %>% 
  ggvis(x = ~xmin_, x2 = ~xmax_, y2 = 0, y = ~count_, fill := "black") %>%
  layer_rects()

## ------------------------------------------------------------------------
mtcars %>% ggvis(~wt, ~mpg) %>% layer_smooths()

# Or equivalently
smoothed <- mtcars %>% compute_smooth(mpg ~ wt)
smoothed %>% ggvis(~pred_, ~resp_) %>% layer_paths()

## ------------------------------------------------------------------------
span <- input_slider(0.2, 1, value = 0.75)
mtcars %>% ggvis(~wt, ~mpg) %>% layer_smooths(span = span)

## ------------------------------------------------------------------------
mtcars %>% 
  ggvis(~wt, ~mpg) %>% 
  layer_smooths() %>% 
  layer_points()

## ------------------------------------------------------------------------
mtcars %>% ggvis(~wt, ~mpg) %>%
  layer_smooths(span = 1) %>%
  layer_smooths(span = 0.3, stroke := "red")

