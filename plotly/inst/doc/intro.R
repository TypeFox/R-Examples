## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  message = FALSE, 
  warning = FALSE, 
  fig.width = 7, 
  fig.height = 3
)

## ------------------------------------------------------------------------
library(plotly)
plot_ly(economics, x = date, y = unemploy / pop)

## ------------------------------------------------------------------------
library(plotly)
plot_ly(economics, x = date, y = unemploy / pop, 
        type = "scatter", mode = "markers+lines")

## ------------------------------------------------------------------------
m <- loess(unemploy / pop ~ as.numeric(date), data = economics)
p <- plot_ly(economics, x = date, y = unemploy / pop, name = "raw") 
add_trace(p, y = fitted(m), name = "loess")

## ------------------------------------------------------------------------
economics %>%
  plot_ly(x = date, y = unemploy / pop) %>% 
  add_trace(y = fitted(m)) %>%
  layout(showlegend = F)

## ------------------------------------------------------------------------
economics %>%
  transform(rate = unemploy / pop) %>%
  plot_ly(x = date, y = rate) %>% 
  subset(rate == max(rate)) %>%
  layout(
    showlegend = F, 
    annotations = list(x = date, y = rate, text = "Peak", showarrow = T)
  )

## ---- eval = FALSE-------------------------------------------------------
#  s <- plot_ly(z = volcano, type = "surface")

## ---- eval = FALSE-------------------------------------------------------
#  plotly_POST(s)

## ------------------------------------------------------------------------
plot_ly(iris, x = Petal.Length, y = Petal.Width, 
        color = Species, mode = "markers")

## ------------------------------------------------------------------------
plot_ly(iris, x = Petal.Length, y = Petal.Width, 
        color = Species, colors = "Set1", mode = "markers")

## ------------------------------------------------------------------------
cols <- RColorBrewer::brewer.pal(9, "Set1")
scales::show_col(cols)

## ------------------------------------------------------------------------
cols <- RColorBrewer::brewer.pal(nlevels(iris$Species), "Set1")
plot_ly(iris, x = Petal.Length, y = Petal.Width, 
        color = Species, colors = cols, mode = "markers")

## ------------------------------------------------------------------------
plot_ly(iris, x = Petal.Length, y = Petal.Width, 
        color = as.ordered(Species), mode = "markers")

## ------------------------------------------------------------------------
plot_ly(iris, x = Petal.Length, y = Petal.Width, 
        color = Sepal.Length, mode = "markers")

## ------------------------------------------------------------------------
plot_ly(iris, x = Petal.Length, y = Petal.Width, 
        color = Sepal.Length, colors = c("#132B43", "#56B1F7"), 
        mode = "markers")

## ------------------------------------------------------------------------
plot_ly(iris, x = Petal.Length, y = Petal.Width, 
        color = Sepal.Length, colors = "PuOr", mode = "markers")

## ------------------------------------------------------------------------
plot_ly(iris, x = Petal.Length, y = Petal.Width, 
        symbol = Species, mode = "markers")

## ------------------------------------------------------------------------
plot_ly(iris, x = Petal.Length, y = Petal.Width, mode = "markers",
        symbol = Species, symbols = c("cross", "square", "triangle-down"))

## ------------------------------------------------------------------------
plot_ly(iris, x = Petal.Length, y = Petal.Width, 
        group = Species, mode = "markers")

## ------------------------------------------------------------------------
iris$id <- as.integer(iris$Species)
p <- plot_ly(iris, x = Petal.Length, y = Petal.Width, group = Species,
             xaxis = paste0("x", id), mode = "markers")
subplot(p)

## ------------------------------------------------------------------------
p2 <- layout(
  p, 
  xaxis = list(range = range(Petal.Length) + c(-0.1, 0.1)),
  yaxis = list(range = range(Petal.Width) + c(-0.1, 0.1))
)
subplot(p2)

## ------------------------------------------------------------------------
layout(
    subplot(p2),
    yaxis2 = list(title = ""), 
    yaxis3 = list(title = "")
)

