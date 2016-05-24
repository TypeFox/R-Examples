dat <- test_counts_groups_summary()
dat[["selected"]] <- rep(FALSE, nrow(dat))
dat[as.numeric(test_count_point[["selected"]]), "selected"] <- TRUE

p <- ggplot(dat, aes(y = run, x = lambda, shape = selected, colour = experiment,
                     linetype = selected, label = group)) +
  geom_point(size = 4) + cool_theme +
  geom_text(aes(x = lambda.up, y = run, size = selected), 
            show_guide = FALSE, hjust = -0.25, vjust = 0) +
  ggtitle("Grouped experiments") +
  scale_y_discrete("Replicate id", labels = dat[["replicate"]] ) +
  scale_x_continuous(expression(lambda)) + 
  coord_cartesian(xlim = c(ifelse(min(dat[["lambda.low"]]) > 0,
                                  min(dat[["lambda.low"]]) * 0.9, 
                                  min(dat[["lambda.low"]]) * 1.1),
                           ifelse(max(dat[["lambda.up"]]) < 0,
                                  max(dat[["lambda.up"]]) * 0.9, 
                                  max(dat[["lambda.up"]]) * 1.1))) +
  scale_size_discrete(guide = FALSE, range = c(5, 7)) + 
  scale_color_discrete("Experiment name") +
  scale_linetype_manual(guide = FALSE, values = c("solid", "dashed")) + 
  scale_shape_manual(guide = FALSE, values = c(15, 18)) + 
  geom_errorbarh(aes(x = lambda, xmin = lambda.low, xmax = lambda.up), 
                size = 1.2, heigth = nlevels(dat[["run"]])/160)
