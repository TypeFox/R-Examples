## ---- echo = FALSE, message = FALSE--------------------------------------
library(knitr)
library(ggvis)
opts_chunk$set(comment = "#>", error = FALSE, tidy = FALSE)
opts_chunk$set(fig.width = 3.5, fig.height = 2.5, dpi = 100)

## ------------------------------------------------------------------------
mtcars %>%
  ggvis(~wt, ~mpg) %>%
  layer_smooths(span = input_slider(0.5, 1, value = 1)) %>%
  layer_points(size := input_slider(100, 1000, value = 100))

## ---- eval=FALSE---------------------------------------------------------
#  mtcars %>% ggvis(x = ~wt) %>%
#      layer_densities(
#        adjust = input_slider(.1, 2, value = 1, step = .1, label = "Bandwidth adjustment"),
#        kernel = input_select(
#          c("Gaussian" = "gaussian",
#            "Epanechnikov" = "epanechnikov",
#            "Rectangular" = "rectangular",
#            "Triangular" = "triangular",
#            "Biweight" = "biweight",
#            "Cosine" = "cosine",
#            "Optcosine" = "optcosine"),
#          label = "Kernel")
#      )

## ---- eval = FALSE-------------------------------------------------------
#  input_slider(-5, 5, label = "lambda", map = function(x) 10 ^ x)

## ------------------------------------------------------------------------
mtcars %>%
  ggvis(~wt, ~mpg, size := input_slider(10, 1000)) %>%
  layer_points(fill := "red") %>%
  layer_points(stroke := "black", fill := NA)

## ------------------------------------------------------------------------
slider <- input_slider(10, 1000)
mtcars %>% ggvis(~wt, ~mpg) %>%
  layer_points(fill := "red", size := slider) %>%
  layer_points(stroke := "black", fill := NA, size := slider)

## ------------------------------------------------------------------------
slider <- input_slider(100, 1000)
mtcars %>% ggvis(~wt, ~mpg) %>%
  layer_points(size := slider) %>% 
  layer_points(fill := "red", opacity := 0.5, size := slider)

mtcars %>% ggvis(~wt, ~mpg) %>%
  layer_points(size := input_slider(100, 1000, label = "black")) %>%
  layer_points(fill := "red", size := input_slider(100, 1000, label = "red"))

## ---- eval = FALSE-------------------------------------------------------
#  # A reactive expression wrapper for input$size
#  input_size <- reactive(input$size)
#  
#  mtcars %>%
#    ggvis(~disp, ~mpg, size := input_size) %>%
#    layer_points() %>%
#    bind_shiny("ggvis", "ggvis_ui")

## ---- eval = FALSE-------------------------------------------------------
#  shinyUI(pageWithSidebar(
#    sidebarPanel(
#      sliderInput("size", "Area", 10, 1000)
#    ),
#    mainPanel(
#      uiOutput("ggvis_ui"),
#      ggvisOutput("ggvis")
#    )
#  ))

