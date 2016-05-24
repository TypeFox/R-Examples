summ <- summary_plot_dat()
dat <- cbind(summ, selected = rep(FALSE, nrow(summary_plot_dat())))
dat[as.numeric(summary_point[["selected"]]), "selected"] <- TRUE

p <- ggplot(dat, aes(x = experiment, y = lambda, shape = selected)) +
  geom_point(size = 4, alpha = 0.6, lty = 2, colour = "blue") + cool_theme +
  geom_boxplot(outlier.colour = NA, fill = adjustcolor("lightgrey", alpha.f = 0.25), shape = 15) + 
  ggtitle("Experiment boxplot") +
  scale_x_discrete("Experiment name") +
  scale_y_continuous(expression(lambda)) + 
  scale_shape_manual(guide = FALSE, values = c(15, 18))
