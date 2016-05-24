if(!is.factor(df[["value"]]))
  df[["value"]] <- as.factor(df[["value"]])

p <- ggplot(df, aes(x = col, y = row , fill = value, shape = selected)) +
  geom_tile(colour = "black", linetype = 2) + cool_theme  +
  geom_point(size = 6) +
  scale_x_discrete("Column") + scale_y_discrete("Row") +
  scale_fill_discrete("Value") +
  scale_shape_manual(guide = FALSE, values = c(NA, 18)) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  theme(panel.border = element_blank(),
        panel.background = element_blank())
