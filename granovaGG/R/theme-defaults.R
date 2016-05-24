theme_granova_1w <- function(base_size = 12) {
  theme_bw() %+replace%
    theme(
      axis.text.x      = element_text(size = 8, lineheight = 0.9, vjust = 0.5, angle = 90, colour = "grey50"),
      axis.text.y      = element_text(size = 8, lineheight = 0.9, hjust = 1, colour = "grey50"),
      axis.title.x     = element_text(size = 8, vjust = 0, colour = "grey20"),
      axis.title.y     = element_text(size = 8, angle = 90, vjust = 0.3, hjust = 0.5, colour = "grey20"),
      axis.line        = element_blank(),

      legend.text      = element_text(size = 8, lineheight = 8),
      # legend.key.size  = unit(0.5, "lines"),

      panel.background = element_rect(colour = NA),
      panel.border     = element_rect(fill = NA, colour = NA),
      panel.grid.minor = element_line(colour = NA, size = 0.25),
      panel.grid.major = element_line(colour = "grey90", size = 0.1),
      plot.title       = element_text(face = "bold", size = 10, vjust = 1)
    )
}

theme_granova_ds <- function(base_size = 12) {
  theme_granova_ds <- theme_bw()

  theme_granova_ds$axis.title.x    <- element_text(size = 10)
  theme_granova_ds$axis.title.y    <- element_text(size = 10, angle = 90)

  # theme_granova_ds$legend.key.size <- unit(1, "lines")
  theme_granova_ds$legend.text     <- element_text(size = 8, lineheight = 8)
  theme_granova_ds$plot.title      <- element_text(face = "bold", size = base_size)

  return(theme_granova_ds)
}

theme_granova_contr <- function(base_size = 12) {
  theme_granova_contr <- theme_bw()
  theme_granova_contr$axis.text.x      <- element_text()
  theme_granova_contr$axis.text.y      <- element_text(hjust = 1)
  theme_granova_contr$axis.title.x     <- element_text(size = 10)
  theme_granova_contr$axis.title.y     <- element_text(size = 10, angle = 90, vjust = 0.3, hjust = 0.5)

  theme_granova_contr$panel.grid.major <- element_blank()
  theme_granova_contr$panel.grid.minor <- element_blank()

  theme_granova_contr$plot.title       <- element_text(face = "bold", size = base_size, vjust = 1)

  return(theme_granova_contr)
}

# theme_gray <- function(base_size = 12) {
#   structure(list(
#     axis.line         = element_blank(),
#     axis.text.x       = element_text(size = base_size * 0.8 , lineheight = 0.9, colour = "grey50", vjust = 1),
#     axis.text.y       = element_text(size = base_size * 0.8, lineheight = 0.9, colour = "grey50", hjust = 1),
#     axis.ticks        = theme_segment(colour = "grey50"),
#     axis.title.x      = element_text(size = base_size, vjust = 0.5),
#     axis.title.y      = element_text(size = base_size, angle = 90, vjust = 0.5),
#     axis.ticks.length = unit(0.15, "cm"),
#     axis.ticks.margin = unit(0.1, "cm"),
#
#     legend.background = theme_rect(colour="white"),
#     legend.key        = theme_rect(fill = "grey95", colour = "white"),
#     legend.key.size   = unit(1.2, "lines"),
#     legend.text       = element_text(size = base_size * 0.8),
#     legend.title      = element_text(size = base_size * 0.8, face = "bold", hjust = 0),
#     legend.position   = "right",
#
#     panel.background = theme_rect(fill = "grey90", colour = NA),
#     panel.border     = element_blank(),
#     panel.grid.major = theme_line(colour = "white"),
#     panel.grid.minor = theme_line(colour = "grey95", size = 0.25),
#     panel.margin     = unit(0.25, "lines"),
#
#     strip.background = theme_rect(fill = "grey80", colour = NA),
#     strip.text.x     = element_text(size = base_size * 0.8),
#     strip.text.y     = element_text(size = base_size * 0.8, angle = -90),
#
#     plot.background = theme_rect(colour = NA, fill = "white"),
#     plot.title      = element_text(size = base_size * 1.2),
#     plot.margin     = unit(c(1, 1, 0.5, 0.5), "lines")
#   ), class          = "options")
# }
