trellis.par.set(theme = col.mosaic(bw = TRUE))
histogram( ~ Sepal.Length, data = iris, breaks = c(4,5,5.5,5.75,6,6.5,7,8,9),type = "count")
trellis.focus('panel', 1, 1)
grid::grid.text(y = .7, 'Never do this!', gp = grid::gpar(col = 'red', cex = 2, alpha = .6))
trellis.unfocus()
trellis.par.set(theme = col.mosaic())

