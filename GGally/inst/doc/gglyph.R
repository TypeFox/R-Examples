## ----echo=FALSE, message=FALSE-------------------------------------------
ignore <- suppressMessages(library(ggplot2))
ignore <- suppressMessages(library(plyr))
ignore <- suppressMessages(library(GGally))
ignore <- lapply(dir(file.path("..", "R"), full.names = TRUE), source)
knitr::opts_chunk$set(fig.width = 9, fig.height = 7, fig.retina = 1)

## ----basic-usage, fig.height=7, fig.width=7------------------------------
data(nasa)
temp.gly <- glyphs(nasa, "long", "day", "lat", "surftemp", height=2.5)
ggplot(temp.gly, ggplot2::aes(gx, gy, group = gid)) +
  add_ref_lines(temp.gly, color = "grey90") +
  add_ref_boxes(temp.gly, color = "grey90") +
  geom_path() +
  theme_bw() +
  labs(x = "", y = "")

