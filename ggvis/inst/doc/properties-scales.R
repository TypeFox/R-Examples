## ---- echo = FALSE, message = FALSE--------------------------------------
library(knitr)
library(ggvis)
opts_chunk$set(comment = "#>", error = FALSE, tidy = FALSE)
opts_chunk$set(fig.width = 3.5, fig.height = 2.5, dpi = 100)

## ----eval = FALSE--------------------------------------------------------
#  geom_point(aes(x = wt, y = mpg), colour = "red", size = 5)

## ---- eval = FALSE-------------------------------------------------------
#  layer_paths(x = ~wt, y = ~mpg, stroke := "red", strokeWidth := 5)

## ---- eval = FALSE-------------------------------------------------------
#  mtcars %>% ggvis(x = ~wt, y = ~mpg) %>%
#    layer_points() %>%
#    layer_model_predictions(model = "lm", stroke = "lm") %>%
#    layer_smooths(stroke = "loess")

## ---- eval = FALSE-------------------------------------------------------
#  mtcars %>% ggvis() %>% layer_lines(strke = ~cyl)
#  #> Error: Unknown properties: strke. Did you mean: stroke?
#  mtcars %>% ggvis(strke = ~cyl) %>% layer_lines()
#  #> Error: Unknown properties: strke. Did you mean: stroke?

## ------------------------------------------------------------------------
df <- data.frame(x = 1:10)
f <- function(n) {
  df %>% ggvis(x = ~x, y = ~x ^ n) %>% layer_paths()
}
f(1)
f(2)
f(4)

## ------------------------------------------------------------------------
prop("x", quote(mpg))
prop("y", ~cyl)

## ------------------------------------------------------------------------
var <- "mpg"
prop("x", as.name(var))

## ------------------------------------------------------------------------
expr <- "mpg / wt"
prop("x", parse(text = expr)[[1]])

## ------------------------------------------------------------------------
# Override the default data limits:
mtcars %>% ggvis(~disp, ~wt) %>%
  layer_points() %>%
  scale_numeric("x", domain = c(50, 500), nice = FALSE) %>%
  scale_numeric("y", domain = c(0, 6), nice = FALSE)

## ---- eval = FALSE-------------------------------------------------------
#  scale_numeric(vis, "x")
#  scale_numeric(vis, "y")
#  scale_nominal(vis, "shape")

## ---- eval = FALSE-------------------------------------------------------
#  scale_numeric(vis, "x", domain = c(10, 100))
#  scale_numeric(vis, "x", trans = "log")

## ------------------------------------------------------------------------
df <- data.frame(x = runif(5), y = runif(5),
  labels = c("a", "b", "b", "a", "b"))
df %>% ggvis(~x, ~y, text := ~labels, font = ~labels, fontSize := 40) %>%
  layer_text() %>%
  scale_ordinal("font", range = c("Helvetica Neue", "Times New Roman"))

## ------------------------------------------------------------------------
mtcars %>% ggvis(y = ~mpg) %>%
  layer_points(prop("x", ~disp, scale = "xdisp")) %>%
  layer_points(prop("x", ~wt, scale = "xwt"), fill := "blue") %>%
  add_axis("x", "xdisp", orient = "bottom") %>%
  add_axis("x", "xwt", orient = "bottom", offset = 20,
    properties = axis_props(labels = list(fill = "blue")))

## ------------------------------------------------------------------------
df <- data.frame(x = 1:5, y = 1:5, a = runif(5), b = -runif(5))

df %>% 
  ggvis(x = ~x, y = ~y, stroke = ~a, fill = ~b, 
    strokeWidth := 5, size := 1000) %>%
  layer_points() %>%
  add_legend("stroke", properties = legend_props(legend = list(y = 50)))

df %>% 
  ggvis(x = ~x, y = ~y, stroke = ~a, prop("fill", ~b, scale = "stroke"),
    strokeWidth := 5, size := 1000) %>%
  layer_points() %>%
  add_legend("stroke", properties = legend_props(legend = list(y = 50)))

