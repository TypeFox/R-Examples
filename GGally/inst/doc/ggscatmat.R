## ----echo=FALSE, message=FALSE-------------------------------------------
ignore <- suppressMessages(library(ggplot2))
ignore <- suppressMessages(library(plyr))
ignore <- lapply(dir(file.path("..", "R"), full.names = TRUE), source)
knitr::opts_chunk$set(fig.width = 9, fig.height = 7, fig.retina = 1)

## ----basic-usage, fig.height=7, fig.width=7------------------------------
data(flea)
ggscatmat(flea, columns = 2:4, color="species", alpha=0.8)

