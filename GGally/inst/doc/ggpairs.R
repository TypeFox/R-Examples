## ----echo=FALSE, message=FALSE-------------------------------------------
ignore <- suppressMessages(library(ggplot2))
ignore <- suppressMessages(library(grid))
ignore <- suppressMessages(library(gtable))
ignore <- lapply(dir(file.path("..", "R"), full.names = TRUE), source)
knitr::opts_chunk$set(fig.width = 9, fig.height = 7, fig.retina = 1)

## ----ggpairs_columns-----------------------------------------------------
data(tips, package = "reshape")
pm <- ggpairs(tips)
pm
## too many plots for this example.  

## reduce the columns being displayed
## these two lines of code produce the same plot matrix
pm <- ggpairs(tips, columns = c(1, 6, 2))
pm <- ggpairs(tips, columns = c("total_bill", "time", "tip"), columnLabels = c("Total Bill", "Time of Day", "Tip"))
pm

## ----ggpairs_mapping-----------------------------------------------------
library(ggplot2)
pm <- ggpairs(tips, mapping = aes(color = sex), columns = c("total_bill", "time", "tip"))
pm

## ----ggpairs_section-----------------------------------------------------
library(ggplot2)
pm <- ggpairs(
  tips, columns = c("total_bill", "time", "tip"),
  lower = list(
    continuous = "smooth",
    combo = "facetdensity",
    mapping = aes(color = time)
  )
)
pm

## ----ggpairs_blank-------------------------------------------------------
pm <- ggpairs(
  tips, columns = c("total_bill", "time", "tip"),
  upper = "blank",
  diag = NULL
)
pm

## ------------------------------------------------------------------------
custom_function <- function(data, mapping, ...){
  # produce ggplot2 object here
}

## ----ggpairs_custom_function---------------------------------------------
my_bin <- function(data, mapping, ..., low = "#132B43", high = "#56B1F7") {
  ggplot(data = data, mapping = mapping) +
    geom_bin2d(...) +
    scale_fill_gradient(low = low, high = high)
}
pm <- ggpairs(
  tips, columns = c("total_bill", "time", "tip"),
  lower = list(
    continuous = my_bin
  )
)
pm

## ----ggpairs_wrap--------------------------------------------------------
pm <- ggpairs(
  tips, columns = c("total_bill", "time", "tip"),
  lower = list(
    combo = wrap("facethist", binwidth = 1),
    continuous = wrap(my_bin, binwidth = c(5, 0.5), high = "red")
  )
)
pm

## ----ggpairs_matrix------------------------------------------------------
pm <- ggpairs(tips, columns = c("total_bill", "time", "tip"))
# retrieve the third row, first column plot
p <- pm[3,1]
p <- p + aes(color = time)
p
pm[3,1] <- p
pm

## ----ggpairs_theme-------------------------------------------------------
pmBW <- pm + theme_bw()
pmBW

