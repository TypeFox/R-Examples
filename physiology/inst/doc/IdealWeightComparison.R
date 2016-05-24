## ----, echo=FALSE, fig.width=7, fig.height=7-----------------------------
library(physiology)
library(ggplot2)
library(reshape2)
hts <- seq(from = 1.5, to = 2.0, by = 0.01)
suppressWarnings({
  ideal_weights <- data.frame(
    "Height" = hts,
    Devine = ideal_weight_Devine(heightm = hts, male = rep(T, times=length(hts))),
    Broca = ideal_weight_Broca(heightm = hts, male = rep(T, times=length(hts))),
    Robinson = ideal_weight_Robinson(heightm = hts, male = rep(T, times=length(hts))),
    Miller = ideal_weight_Miller(heightm = hts, male = rep(T, times=length(hts))),
    Lemmens = ideal_weight_Lemmens(heightm = hts)
    )
  })
p <- ggplot(melt(ideal_weights, id.vars = "Height"),
            aes(x = Height, y = value, group = variable)) +
  geom_line(aes(colour = variable)) +
  #  scale_x_continuous(breaks = seq(0.5, 3, 0.5), limits = c(0.5, 3)) +
  #scale_y_continuous(breaks = seq(0, 5, 1), limits = c(0, 5)) +
  labs(x = "height, m", y = "ideal weight, kg") +
  #  ggtitle("Comparison of ideal weight formulae") +
  theme_bw() +
  theme(legend.position = c(0.25, 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = rel(1.1), face = "bold"))
p

