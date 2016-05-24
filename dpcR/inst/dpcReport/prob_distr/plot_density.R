p <- ggplot(dens, aes(x = x, y = y)) + geom_line(colour = "lightskyblue1", size = 1.2) + 
  geom_area(aes(fill = conf_up)) + 
  geom_area(aes(fill = conf_low)) +
  scale_fill_manual(values = c("FALSE" = NA, "TRUE" = adjustcolor("cyan4", alpha.f = 0.5)), guide = FALSE) +
  cool_theme + 
  scale_y_continuous("Density")

p <- if(input[["density_plot_avg"]]) {
  p + scale_x_continuous(expression(lambda))
} else {
  p + scale_x_continuous("k")
}

if(input[["density_plot_bars"]])
  p <- p + geom_bar(stat = "identity", fill = adjustcolor("lightskyblue1", alpha.f = 0.5))