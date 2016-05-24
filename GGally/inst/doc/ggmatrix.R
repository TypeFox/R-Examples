## ----echo=FALSE, message=FALSE-------------------------------------------
ignore <- suppressMessages(library(ggplot2))
ignore <- suppressMessages(library(grid))
ignore <- suppressMessages(library(gtable))
ignore <- lapply(dir(file.path("..", "R"), full.names = TRUE), source)
knitr::opts_chunk$set(fig.width = 9, fig.height = 7, fig.retina = 1)

## ----ggmatrix_genExample-------------------------------------------------
plotList <- list()
for (i in 1:6) {
  plotList[[i]] <- ggally_text(paste("Plot #", i, sep = ""))
}

# bare minimum of plotList, nrow, and ncol
pm <- ggmatrix(plotList, 2, 3)
pm

# provide more information
pm <- ggmatrix(
  plotList,
  nrow = 2, ncol = 3,
  xAxisLabels = c("A", "B", "C"),
  yAxisLabels = c("D", "E"),
  title = "Matrix Title"
)
pm

# display plots in column order
pm <- ggmatrix(
  plotList,
  nrow = 2, ncol = 3,
  xAxisLabels = c("A", "B", "C"),
  yAxisLabels = c("D", "E"),
  title = "Matrix Title",
  byrow = FALSE
)
pm

## ----ggmatrix_place------------------------------------------------------
pm <- ggmatrix(
  plotList,
  nrow = 2, ncol = 3,
  xAxisLabels = c("A", "B", "C"),
  yAxisLabels = c("D", "E"),
  title = "Matrix Title"
)
pm
p2 <- pm[1,2]
p3 <- pm[1,3]
p2
p3
pm[1,2] <- p3
pm[1,3] <- p2
pm


## ----ggmatrix_theme------------------------------------------------------
pm <- ggmatrix(
  plotList,
  nrow = 2, ncol = 3,
  xAxisLabels = c("A", "B", "C"),
  yAxisLabels = c("D", "E"),
  title = "Matrix Title",
  byrow = FALSE
)
pm <- pm + theme_bw()
pm

## ----ggmatrix_axisControl------------------------------------------------
pm <- ggmatrix(
  plotList, nrow = 2, ncol = 3,
  xAxisLabels = c("A", "B", "C"),
  yAxisLabels = c("D", "E"),
  title = "No Left Plot Axis",
  showYAxisPlotLabels = FALSE
)
pm
pm <- ggmatrix(
  plotList, nrow = 2, ncol = 3,
  xAxisLabels = c("A", "B", "C"),
  yAxisLabels = c("D", "E"),
  title = "No Bottom Plot Axis",
  showXAxisPlotLabels = FALSE
)
pm
pm <- ggmatrix(
  plotList, nrow = 2, ncol = 3,
  xAxisLabels = c("A", "B", "C"),
  yAxisLabels = c("D", "E"),
  title = "No Plot Axes",
  showAxisPlotLabels = FALSE
)
pm

## ----ggmatrix_stripControl-----------------------------------------------
data(tips, package = "reshape")
plotList <- list(
  qplot(total_bill, tip, data = subset(tips, smoker == "No" & sex == "Female")) +
    facet_grid(time ~ day),
  qplot(total_bill, tip, data = subset(tips, smoker == "Yes" & sex == "Female")) +
    facet_grid(time ~ day),
  qplot(total_bill, tip, data = subset(tips, smoker == "No" & sex == "Male")) +
    facet_grid(time ~ day),
  qplot(total_bill, tip, data = subset(tips, smoker == "Yes" & sex == "Male")) +
    facet_grid(time ~ day)
)


pm <- ggmatrix(
  plotList, nrow = 2, ncol = 2,
  yAxisLabels = c("Female", "Male"),
  xAxisLabels = c("Non Smoker", "Smoker"),
  title = "Total Bill vs Tip",
  showStrips = NULL # default
)
pm
pm <- ggmatrix(
  plotList, nrow = 2, ncol = 2,
  yAxisLabels = c("Female", "Male"),
  xAxisLabels = c("Non Smoker", "Smoker"),
  title = "Total Bill vs Tip",
  showStrips = TRUE
)
pm
pm <- ggmatrix(
  plotList, nrow = 2, ncol = 2,
  yAxisLabels = c("Female", "Male"),
  xAxisLabels = c("Non Smoker", "Smoker"),
  title = "Total Bill vs Tip",
  showStrips = FALSE
)
pm

## ----ggmatrix_customPrinting---------------------------------------------
data(presidential, package = "ggplot2")
plotList <- list(
  qplot(name, data = presidential, geom = "bar", fill = party) + ylim(c(0,6)) + coord_flip(),
  qplot(party, data = presidential, geom = "bar", fill = party) + coord_flip()
)
pm <- ggmatrix(
  plotList,
  2, 1,
  yAxisLabels = c("Presidents", "Party"),
  xAxisLabels = c("Count"),
  title = "President Term Length"
)
pm # default spacing
# use print.ggmatrix parameters
print(pm, leftWidthProportion = 0.12, spacingProportion = 0.05)

