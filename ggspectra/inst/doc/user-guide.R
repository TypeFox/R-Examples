## ------------------------------------------------------------------------
library(ggplot2)
library(photobiology)
library(photobiologyWavebands)
library(ggspectra)

## ----setup, include = FALSE, eval = FALSE--------------------------------
#  library(svglite)
#  knitr::opts_chunk$set(
#    dev = "svglite",
#    fig.ext = ".svg"
#  )

## ---- include=FALSE, echo=FALSE------------------------------------------
library(knitr)
opts_chunk$set(fig.path = 'figure/pos-', fig.align = 'center', fig.show = 'hold',
               fig.width = 7, fig.height = 4)
options(warnPartialMatchArgs = FALSE)

## ------------------------------------------------------------------------
two_suns.spct <- rbindspct(list(sun1 = sun.spct, sun2 = sun.spct * 2))

## ------------------------------------------------------------------------
theme_set(theme_bw())

## ------------------------------------------------------------------------
ggplot(sun.spct) + geom_line()

## ------------------------------------------------------------------------
ggplot(two_suns.spct) + aes(color = spct.idx) + geom_line()

## ------------------------------------------------------------------------
ggplot(two_suns.spct, aes(w.length, s.e.irrad, color = spct.idx)) + geom_line()

## ------------------------------------------------------------------------
ggplot(sun.spct, unit.out = "photon") + geom_line()

## ------------------------------------------------------------------------
ggplot(yellow_gel.spct) + geom_line()

## ------------------------------------------------------------------------
ggplot(yellow_gel.spct, plot.qty = "absorbance") + geom_line()

## ------------------------------------------------------------------------
ggplot(sun.spct) + geom_line() + stat_peaks()

## ------------------------------------------------------------------------
ggplot(sun.spct, unit.out = "photon") + geom_line() + stat_peaks()

## ------------------------------------------------------------------------
ggplot(sun.spct) + geom_line() + stat_valleys()

## ------------------------------------------------------------------------
ggplot(sun.spct) + geom_line() + 
  stat_peaks(shape = 21) + scale_fill_identity()

## ------------------------------------------------------------------------
ggplot(sun.spct) + geom_line() + 
  stat_peaks(span = 21, shape = 4, color = "red", size = 2) +
  stat_peaks(span = 21, color = "red", geom = "rug", sides = "b")

## ------------------------------------------------------------------------
ggplot(sun.spct) + geom_line() + 
  stat_peaks(span = 21, geom = "text")

## ------------------------------------------------------------------------
ggplot(sun.spct) + geom_line() + 
  stat_peaks(shape = 21, span = 25, size = 2) + scale_fill_identity() +
  stat_peaks(geom = "label", span = 25, vjust = "bottom", size = 3, 
             color = "white")

## ------------------------------------------------------------------------
ggplot(sun.spct) + geom_line() + 
  stat_peaks(span = NULL, geom = "vline", linetype = "dotted") +
  stat_peaks(span = NULL, geom = "hline", linetype = "dotted")

## ------------------------------------------------------------------------
ggplot(sun.spct) + geom_line() + 
  stat_peaks(span = 21, geom = "text", aes(label = ..y.label..))

## ------------------------------------------------------------------------
ggplot(sun.spct) + geom_line() + 
  stat_peaks(span = NULL) +
  stat_peaks(span = NULL, geom = "text", vjust = -0.5, size = 2.5,
             aes(label = paste(..y.label.., "at", ..x.label.. , "nm")))

## ------------------------------------------------------------------------
ggplot(sun.spct) + geom_line() + 
  stat_peaks(span = 21, geom = "point", colour = "red") +
  stat_valleys(span = 21, geom = "point", colour = "blue") +
  stat_peaks(span = 51, geom = "text", colour = "red", 
             vjust = -0.3, label.fmt = "%3.0f nm") +
  stat_valleys(span = 51, geom = "text", colour = "blue", 
               vjust = 1.2, label.fmt = "%3.0f nm")

## ------------------------------------------------------------------------
ggplot(two_suns.spct) + aes(color = spct.idx) +
  geom_line() + ylim(NA, 1.8) +
  stat_peaks(span = NULL, color = "black") +
  stat_peaks(span = NULL, geom = "text", vjust = -0.5, size = 3, 
             color = "black", 
             aes(label = paste(..y.label.., "at", ..x.label.. , "nm"))) +
  facet_grid(spct.idx~.)

## ------------------------------------------------------------------------
ggplot(sun.spct) + 
  stat_color() + scale_color_identity()

## ------------------------------------------------------------------------
ggplot(sun.spct) +
  geom_line() +
  stat_color(shape = 21, color = "black") + 
  scale_fill_identity()

## ------------------------------------------------------------------------
ggplot(sun.spct) + 
  stat_color(geom = "bar") + 
  geom_point(shape = 21, color = "black", stroke = 1.2, fill = "white") +
  scale_fill_identity() + 
  scale_color_identity() + 
  theme_bw()

## ------------------------------------------------------------------------
ggplot(two_suns.spct) + aes(shape = spct.idx) +
  stat_color() + scale_color_identity() +
  geom_line() + 
  facet_grid(spct.idx~., scales = "free_y")

## ------------------------------------------------------------------------
ggplot(sun.spct) + geom_line() + stat_wl_summary()

## ------------------------------------------------------------------------
ggplot(sun.spct) + 
  stat_wl_summary(range = c(300,350), geom = "rect") +
  geom_line()

## ------------------------------------------------------------------------
ggplot(sun.spct) +
  geom_line() + 
  stat_wl_summary(geom = "hline", color = "red") +
  stat_wl_summary(label.fmt = "Mean = %.3g", color = "red", vjust = -0.3)

## ------------------------------------------------------------------------
ggplot(sun.spct) +
  stat_wl_summary(range = c(400,500), geom = "rect", alpha = 0.2, fill = color(450)) +
  stat_wl_summary(range = c(400,500), label.fmt = "Mean = %.3g", vjust = -0.3, geom = "text") + 
  geom_line()

## ------------------------------------------------------------------------
ggplot(two_suns.spct) + aes(color = spct.idx) +
  geom_line() + 
  stat_wl_summary(geom = "hline") +
  stat_wl_summary(label.fmt = "Mean = %.3g", vjust = 1.2, show.legend = FALSE) +
  facet_grid(spct.idx~.)

## ------------------------------------------------------------------------
ggplot(two_suns.spct) + aes(color = spct.idx) +
  geom_line() + 
  stat_wl_summary(geom = "hline") +
  stat_wl_summary(label.fmt = "Mean = %.3g", vjust = 1.2, show.legend = FALSE) +
  facet_grid(spct.idx~., scales = "free_y")

## ------------------------------------------------------------------------
ggplot(sun.spct) +
  stat_wb_mean(w.band = PAR(), geom = "rect", alpha = 0.5) +
  stat_wb_mean(w.band = PAR(), geom = "text") +
  geom_line() + 
  scale_color_identity() + 
  scale_fill_identity() +  
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.spct) +
  stat_wb_mean(w.band = c(400,500), geom = "rect", alpha = 0.5) +
  stat_wb_mean(w.band = c(400,500), geom = "text") +
  geom_line() + 
  scale_color_identity() + 
  scale_fill_identity() +
  theme_bw()

## ------------------------------------------------------------------------
ggplot(yellow_gel.spct) +
  stat_wb_mean(w.band = Plant_bands(), geom = "rect", alpha = 0.7) +
  stat_wb_mean(w.band = Plant_bands(), angle = 90, hjust = 0, geom = "text", ypos.fixed = 0.1, 
                 aes(label = paste(..wb.name.., " = ", ..y.label.., sep = ""))) +
  geom_line() + 
  scale_fill_identity() +
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.spct) +
  stat_wb_total(w.band = PAR(), geom = "rect", alpha = 0.5) +
  stat_wb_total(w.band = PAR(), geom = "text") +
  geom_line() + 
  scale_color_identity() + 
  scale_fill_identity() +  
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.spct) +
  stat_wb_total(w.band = c(400,500), geom = "rect", alpha = 0.5) +
  stat_wb_total(w.band = c(400,500), geom = "text") +
  geom_line() + 
  scale_color_identity() + 
  scale_fill_identity() +
  theme_bw()

## ------------------------------------------------------------------------
ggplot(yellow_gel.spct) +
  stat_wb_total(w.band = Plant_bands(), geom = "rect", alpha = 0.4, color = "white",
                ypos.fixed = 1) +
  stat_wb_mean(w.band = Plant_bands(), geom = "text", label.fmt = "%1.2f",
               ypos.fixed = 1) +
  geom_line() + 
  scale_fill_identity() +
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.spct) +
  stat_wb_irrad(w.band = PAR(), geom = "rect", alpha = 0.5,
                unit.in = "energy", time.unit = "second") +
  stat_wb_irrad(w.band = PAR(), geom = "text",
                unit.in = "energy", time.unit = "second") +
  geom_line() + 
  scale_color_identity() + 
  scale_fill_identity() +  
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.spct, unit.out = "photon") +
  stat_wb_irrad(w.band = PAR(), geom = "rect", alpha = 0.5,
                unit.in = "photon", time.unit = "second") +
  stat_wb_irrad(w.band = PAR(), geom = "text",
                unit.in = "photon", time.unit = "second", 
                aes(label = sprintf("%s = %.3g", ..wb.name.., ..yint.. * 1e6))) +
  geom_line() + 
  scale_color_identity() + 
  scale_fill_identity() +  
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.spct) +
  stat_wb_irrad(w.band = CIE(), geom = "rect", alpha = 0.5,
                unit.in = "energy", time.unit = "second") +
  stat_wb_irrad(w.band = CIE(), geom = "text",
                unit.in = "energy", time.unit = "second") +
  geom_line() + 
  scale_color_identity() + 
  scale_fill_identity() +  
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.daily.spct) +
  stat_wb_irrad(w.band = CIE(), geom = "rect", alpha = 0.5,
                unit.in = "energy", time.unit = "day") +
  stat_wb_irrad(w.band = CIE(), geom = "text",
                unit.in = "energy", time.unit = "day", label.mult = 1e-3) +
  geom_line() + 
  scale_color_identity() + 
  scale_fill_identity() +  
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.spct, unit.out = "photon") +
  stat_wb_irrad(w.band = VIS_bands(), geom = "rect", alpha = 0.5,
                unit.in = "photon", time.unit = "second") +
  stat_wb_irrad(w.band = VIS_bands(), geom = "text",
                unit.in = "photon", time.unit = "second", 
                label.mult = 1e6, size = rel(2)) +
  geom_line() + 
  scale_color_identity() + 
  scale_fill_identity() +  
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.spct) +
  stat_wb_sirrad(w.band = PAR(), geom = "rect", alpha = 0.5,
                unit.in = "energy", time.unit = "second") +
  stat_wb_sirrad(w.band = PAR(), geom = "text",
                unit.in = "energy", time.unit = "second") +
  geom_line() + 
  scale_color_identity() + 
  scale_fill_identity() +  
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.spct, unit.out = "photon") +
  stat_wb_sirrad(w.band = PAR(), geom = "rect", alpha = 0.5,
                unit.in = "photon", time.unit = "second") +
  stat_wb_sirrad(w.band = PAR(), geom = "text",
                unit.in = "photon", time.unit = "second", 
                aes(label = sprintf("Total %s = %.3g", ..wb.name.., ..yint.. * 1e6))) +
  stat_wb_sirrad(w.band = PAR(), geom = "text",
                unit.in = "photon", time.unit = "second", 
                mapping = aes(label = sprintf("Mean %s = %.3g", ..wb.name.., ..ymean.. * 1e6)), 
                ypos.mult = 0.45) +
  geom_line() + 
  scale_color_identity() + 
  scale_fill_identity() +  
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.spct) +
  stat_wb_sirrad(w.band = CIE(), geom = "rect", alpha = 0.5,
                unit.in = "energy", time.unit = "second") +
  stat_wb_sirrad(w.band = CIE(), geom = "text",
                unit.in = "energy", time.unit = "second") +
  geom_line() + 
  scale_color_identity() + 
  scale_fill_identity() +  
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.daily.spct) +
  stat_wb_sirrad(w.band = CIE(), geom = "rect", alpha = 0.5,
                unit.in = "energy", time.unit = "day") +
  stat_wb_sirrad(w.band = CIE(), geom = "text",
                unit.in = "energy", time.unit = "day", label.mult = 1e-3) +
  geom_line() + 
  scale_color_identity() + 
  scale_fill_identity() +  
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.spct) +
  stat_wb_sirrad(w.band = VIS_bands(), geom = "rect", alpha = 0.5,
                unit.in = "energy", time.unit = "second") +
  stat_wb_sirrad(w.band = VIS_bands(), angle = 90, geom = "text",
                unit.in = "energy", time.unit = "second", ypos.fixed = 0.05, hjust = 0,
                 aes(label = paste(..wb.name.., ..y.label.., sep = " = "),
                     color = "black")) +
  geom_line() + 
  scale_color_identity() + 
  scale_fill_identity() +  
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.spct, unit.out = "photon") +
  stat_wb_sirrad(w.band = VIS_bands(), geom = "rect", alpha = 0.5,
                unit.in = "photon", time.unit = "second") +
  stat_wb_irrad(w.band = VIS_bands(), angle = 90, geom = "text", label.mult = 1e6,
                unit.in = "photon", time.unit = "second", ypos.fixed = 1e-7, hjust = 0,
                 aes(label = paste(..wb.name.., ..y.label.., sep = " = "),
                     color = "black")) +
  geom_line() + 
  scale_color_identity() + 
  scale_fill_identity() +  
  theme_bw()

## ------------------------------------------------------------------------
my.bands <- split_bands(c(300,800), length.out = 10)
ggplot(sun.spct) +
  stat_wb_sirrad(w.band = my.bands, geom = "rect", alpha = 0.5,
                unit.in = "energy", time.unit = "second") +
  stat_wb_irrad(w.band = my.bands, angle = 90, geom = "text",
                unit.in = "energy", time.unit = "second", ypos.fixed = 0.05, hjust = 0,
                color = "black") +
  geom_line() + 
  scale_color_identity() + 
  scale_fill_identity() +  
  theme_bw()

## ------------------------------------------------------------------------
my.bands <- split_bands(c(300,800), length.out = 10)
ggplot(sun.spct, unit.out = "photon") +
  stat_wb_sirrad(w.band = my.bands, geom = "rect", alpha = 0.5,
                unit.in = "photon", time.unit = "second") +
  stat_wb_sirrad(w.band = my.bands, geom = "text", angle = 90,
                unit.in = "photon", time.unit = "second", label.mult = 1e6,
                color = "black") +
  geom_line() + 
  scale_color_identity() + 
  scale_fill_identity() +  
  theme_bw()

## ------------------------------------------------------------------------
my.data <- data.frame(x = 300:800)
ggplot(my.data, aes(x)) + stat_wl_strip(ymin = -1, ymax = 1) +
    scale_fill_identity()

## ------------------------------------------------------------------------
ggplot(sun.spct) +
  geom_line() +
  stat_wl_strip(ymin = -Inf, ymax = -0.025) + 
  scale_fill_identity() +
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.spct) +
  geom_line() + 
  stat_wl_strip(w.band = VIS_bands(), ymin = -Inf, ymax = -0.025) + 
  scale_fill_identity() +
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.spct) +
  stat_wl_strip(w.band = VIS_bands(), ymin = -Inf, ymax = Inf, alpha = 0.4) + 
  scale_fill_identity() +
  geom_line() +
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.spct) +
  stat_wl_strip(alpha = 0.4, ymin = -Inf, ymax = Inf) + 
  scale_fill_identity() +
  geom_line() +
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.spct, unit.out = "photon") +
  stat_wl_strip(alpha = 0.4, ymin = -Inf, ymax = Inf) + 
  stat_wb_total(w.band = PAR(), color = "black", geom = "rect") +
  stat_wb_total(w.band = PAR(), label.mult = 1e6,
                geom = "text", color = "black",
                aes(label = paste(..wb.name.., " = ", ..y.label.., sep = ""))) +  
  geom_line() + 
  scale_fill_identity() +
  theme_bw()

## ------------------------------------------------------------------------
ggplot(sun.spct) + geom_spct()

## ------------------------------------------------------------------------
ggplot(sun.spct) + 
  geom_spct(fill = color(sun.spct))

## ------------------------------------------------------------------------
ggplot(sun.spct * yellow_gel.spct) + 
  geom_spct(fill = color(sun.spct * yellow_gel.spct))

## ------------------------------------------------------------------------
ggplot(sun.spct) + 
  wl_guide(alpha = 0.4) +
  geom_line()

## ------------------------------------------------------------------------
ggplot(sun.spct) + 
  wl_guide(ymax = -0.025) +
  geom_line() 

## ------------------------------------------------------------------------
ggplot(sun.spct) + 
  wl_guide() +
  geom_spct(alpha = 0.75, colour = "white", size = 1) 

## ------------------------------------------------------------------------
ggplot(sun.spct) + 
  wl_guide() +
  geom_line(size = 2, colour = "white") +
  geom_line(size = 1, colour = "black") +
  geom_hline(yintercept = 0, colour = "grey92") +
  theme_bw()

## ------------------------------------------------------------------------
plot(sun.spct)

## ------------------------------------------------------------------------
plot(sun.spct, label.qty = "mean")

## ------------------------------------------------------------------------
plot(sun.spct, label.qty = "contribution")

## ------------------------------------------------------------------------
plot(sun.spct, label.qty = "relative")

## ------------------------------------------------------------------------
plot(sun.spct, label.qty = "mean", unit.out = "photon")

## ------------------------------------------------------------------------
plot(sun.spct, label.qty = "relative.pc")

## ------------------------------------------------------------------------
plot(sun.spct, unit.out = "photon")

## ------------------------------------------------------------------------
plot(sun.spct, annotations = c("segments", "labels", "color.guide"))

## ------------------------------------------------------------------------
plot(sun.spct, 
     annotations = c("segments", "labels", "summaries", "color.guide"))

## ------------------------------------------------------------------------
plot(sun.spct, range = VIS())

## ------------------------------------------------------------------------
plot(sun.spct, w.band = PAR())

## ------------------------------------------------------------------------
plot(sun.spct, w.band = CIE())

## ------------------------------------------------------------------------
plot(sun.daily.spct)

## ------------------------------------------------------------------------
plot(two_suns.spct, label.qty = "mean") + facet_wrap(~spct.idx)

## ------------------------------------------------------------------------
plot(two_suns.spct) + aes(linetype = spct.idx)

## ------------------------------------------------------------------------
plot(yellow_gel.spct, pc.out = TRUE)

## ------------------------------------------------------------------------
plot(sun.spct) + geom_spct(fill = color(sun.spct)) + 
  geom_spct(data = yellow_gel.spct * sun.spct, color = "black", 
            fill = color(yellow_gel.spct * sun.spct))

## ------------------------------------------------------------------------
plot(sun.spct) + geom_spct(fill = color(sun.spct)) + 
  geom_spct(data = yellow_gel.spct * sun.spct, color = "black", 
            fill = color(yellow_gel.spct * sun.spct)) +
  stat_peaks(data = yellow_gel.spct * sun.spct, color = "yellow", 
             ignore_threshold = 0.1, span = 21)

## ------------------------------------------------------------------------
theme_set(theme_bw(8))
plot(yellow_gel.spct, annotations = "colour.guide")

## ------------------------------------------------------------------------
two_suns.mspct <- source_mspct(list(sun1 = sun.spct, sun2 = sun.spct * 2))

## ---- fig.width = 7, fig.height = 8--------------------------------------
multiplot(plotlist = mslply(two_suns.mspct, plot))

## ---- fig.width = 7, fig.height = 8--------------------------------------
plot(rbindspct(two_suns.mspct)) + facet_wrap(~spct.idx, ncol = 1)

## ------------------------------------------------------------------------
plot(rbindspct(two_suns.mspct), 
     annotations = c("color.guide", "labels", "boxes", "peaks")) + 
  aes(linetype = spct.idx)

## ------------------------------------------------------------------------
ggplot(rbindspct(two_suns.mspct)) + 
  aes(linetype = spct.idx) +
  wl_guide(ymax = -0.05) +
  geom_line()

## ------------------------------------------------------------------------
plot(VIS())

## ------------------------------------------------------------------------
plot(CIE(), range = CIE())

## ------------------------------------------------------------------------
plot(DNA_N(), range = c(270, 420))

