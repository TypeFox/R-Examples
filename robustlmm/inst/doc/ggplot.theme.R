require(grid)
## set ggplot theme
theme <- theme_bw(base_size = 10)
theme$legend.key.size <- unit(1, "lines") 
theme$panel.grid.minor <- element_line(colour = "grey95", size = 0.5)
theme$plot.margin <- unit(c(1/2, 1/10, 0, 0), "lines")
theme$axis.ticks.margin <- unit(0.1, "lines")
theme_set(theme)
