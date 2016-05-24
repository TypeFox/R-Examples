## ---- echo = TRUE--------------------------------------------------------
require("ggplot2")
huron <- data.frame(year = 1875:1972, level = as.vector(LakeHuron))
h <- ggplot(huron, aes(year))
h + RcmdrPlugin.KMggplot2::geom_stepribbon(
      aes(ymin = level - 1, ymax = level + 1),
      fill = "grey70"
    ) +
    geom_step(aes(y = level))
h + geom_ribbon(
      aes(ymin = level - 1, ymax = level + 1),
      fill = "grey70"
    ) +
    geom_line(aes(y = level))

