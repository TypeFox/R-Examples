## ---- echo = FALSE, message = FALSE--------------------------------------
library(knitr)
opts_chunk$set(comment = "#>", error = FALSE, tidy = FALSE)
opts_chunk$set(fig.width = 2, fig.height = 1.25, dpi = 100)

## ---- message = FALSE----------------------------------------------------
library(ggvis)
library(dplyr)

## ------------------------------------------------------------------------
# The first few rows of mtcars
head(mtcars)
mtcars %>% ggvis(~wt, ~mpg) %>% layer_points()

## ------------------------------------------------------------------------
mtcars %>% 
  ggvis(~wt, ~mpg) %>% 
  layer_points(size := 25, shape := "diamond", stroke := "red", fill := NA)

## ---- message = FALSE----------------------------------------------------
mtcars %>% 
  ggvis(~wt, ~mpg) %>%
  layer_points() %>%
  layer_smooths()

## ---- message = FALSE----------------------------------------------------
mtcars %>% 
  ggvis(~wt, ~mpg) %>%
  layer_points() %>%
  layer_model_predictions(model = "lm", se = TRUE)

## ---- message = FALSE----------------------------------------------------
mtcars %>% 
  ggvis(~wt, ~mpg) %>% 
  layer_points(fill = ~factor(cyl))

## ---- message = FALSE----------------------------------------------------
mtcars %>% 
  ggvis(~wt, ~mpg, fill = ~factor(cyl)) %>% 
  layer_points() %>% 
  group_by(cyl) %>% 
  layer_model_predictions(model = "lm")

## ------------------------------------------------------------------------
# The first few rows
head(pressure)

## ---- message = FALSE, fig.width = 4-------------------------------------
pressure %>% 
  ggvis(~temperature, ~pressure) %>%
  layer_bars()

## ---- message = FALSE, fig.width = 4-------------------------------------
pressure %>% 
  ggvis(~temperature, ~pressure) %>%
  layer_bars(width = 10)

## ---- message = FALSE, fig.width = 4-------------------------------------
# First, modify the pressure data set so that the x variable is a factor
pressure2 <- pressure %>% mutate(temperature = factor(temperature))

pressure2 %>% ggvis(~temperature, ~pressure) %>%
  layer_bars()

## ---- message = FALSE----------------------------------------------------
pressure %>% ggvis(~temperature, ~pressure) %>% layer_lines()

## ---- message = FALSE----------------------------------------------------
pressure %>% ggvis(~temperature, ~pressure) %>%
  layer_points() %>% 
  layer_lines()

## ------------------------------------------------------------------------
# The first few rows
head(faithful)

## ---- message = FALSE----------------------------------------------------
faithful %>% ggvis(~eruptions) %>% layer_histograms()

## ---- message = FALSE----------------------------------------------------
faithful %>% ggvis(~eruptions) %>% layer_histograms(width=0.5, boundary=0)
faithful %>% ggvis(~eruptions) %>% layer_histograms(width=0.5, center=0)

## ---- message = FALSE----------------------------------------------------
faithful %>% ggvis(~eruptions, fill := "#fff8dc") %>%
  layer_histograms(width = 0.25)

## ---- message = FALSE----------------------------------------------------
cocaine %>% ggvis(~month, fill := "#fff8dc") %>%
  layer_histograms() %>%
  add_axis("x", title = "month") %>%
  add_axis("y", title = "count")

## ---- message = FALSE----------------------------------------------------
cocaine %>% ggvis(~month, fill := "#fff8dc") %>%
  layer_histograms(width = 1, center = 0) %>%
  add_axis("x", title = "month") %>%
  add_axis("y", title = "count")

## ----message=FALSE, warning=FALSE----------------------------------------
mtcars %>% ggvis(~factor(cyl), ~mpg) %>% layer_boxplots()

