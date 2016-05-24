## testing ggplot guides
ord = c(2,3,1)

qplot(data = mpg, x = displ, y = cty, size = hwy, colour = cyl, shape = drv) +
  guides(colour = guide_colourbar(order = ord[1]),
  alpha = guide_legend(order = ord[2]),
   size = guide_legend(order = ord[3])) +
  ggtitle("ord " %&% paste(c("colour","alpha","size")[ord], collapse=","))
