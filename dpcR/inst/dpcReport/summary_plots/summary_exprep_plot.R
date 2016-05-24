summ <- summary_exprep_plot_dat()
dat <- cbind(summ, selected = rep(FALSE, nrow(summary_exprep_plot_dat())))
dat[as.numeric(summary_exprep_point[["selected"]]), "selected"] <- TRUE

p <- ggplot(dat, aes(y = exprep, x = lambda, shape = selected, colour = experiment,
                     ymin = exprep, ymax = exprep, linetype = selected)) +
  geom_point(size = 4) + cool_theme +
  ggtitle(paste0("Experiment/replicate scatter chart\nCI method: ", cap1L(input[["CI_method"]]))) +
  scale_y_discrete("Replicate id", labels = dat[["replicate"]] ) +
  scale_x_continuous(expression(lambda)) + 
  scale_color_discrete("Experiment name") +
  scale_linetype_manual(guide = FALSE, values = c("solid", "dashed")) + 
  scale_shape_manual(guide = FALSE, values = c(15, 18)) + 
  geom_errorbarh(aes(x = lambda, xmin = lambda.low, xmax = lambda.up), 
                 size = 1.2, heigth = nlevels(dat[["exprep"]])/160)