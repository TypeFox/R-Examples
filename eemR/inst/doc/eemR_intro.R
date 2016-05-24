## ----setup, warning=FALSE, message=FALSE, echo=FALSE---------------------
library(extrafont)
loadfonts()

library(ggplot2)
library(dplyr)
library(tidyr)
library(eemR)
library(plot3D)

theme_set(theme_bw(base_size = 8, base_family = "Open Sans"))

## ----fig1, echo = FALSE, fig.width=5, fig.height=5, dev='png', fig.align='center', fig.cap = "Example of an excitation-emission fluorescence matrix (EEM). The diagonal structure with high fluorescence corresponds to the first order of Rayleigh scattering."----
file <- system.file("extdata/cary/scans_day_1", "sample3.csv", package = "eemR")

eem <- eem_read(file)

persp3D(x = eem[[1]]$ex,
        y = eem[[1]]$em,
        z = t(eem[[1]]$x),
        theta = -37.5,
        phi = 20,
        facets = TRUE,
        xlab = "Excitation (nm)",
        ylab = "Emisison (nm)",
        zlab = "Fluorescence (A.U.)",
        ticktype = "detailed",
        box = TRUE,
        expand = 0.5)

## ----fig2, echo=FALSE, fig.height=3, fig.width=5, warning=FALSE, dev='png', fig.align='center', fig.cap = "Emission fluorescence emitted at excitation $ex = 350$. First order of Rayleigh and Raman scattering regions are identified in blue and red."----

file <- system.file("extdata/cary/scans_day_1", "sample3.csv", package = "eemR")

x_raw <- eem_read(file, recursive = TRUE)

x_cor <- eem_remove_scattering(x_raw, "rayleigh", 1, 12)
x_cor <- eem_remove_scattering(x_cor, "raman", 1, 5)

ex <- 350
em <- x_raw[[1]]$em

em_raw <- x_raw[[1]]$x[, which(x_raw[[1]]$ex == 350)]
em_cor <- x_cor[[1]]$x[, which(x_cor[[1]]$ex == 350)]

df <- data.frame(em, em_raw, em_cor)
df$em_raman <- df$em_raw
df$em_raman[df$em <= 375] <- NA
df$em_rayleigh <- df$em_raw
df$em_rayleigh[df$em > 375] <- NA

ggplot(df, aes(x = em)) +
  geom_line(aes(y = em_rayleigh, color = "Rayleigh scattering"), size = 0.75) +
  geom_line(aes(y = em_raman, color = "Raman scattering"), size = 0.75) +
  geom_line(aes(y = em_cor, color = "Fluorescence signal"), size = 0.75) +
  labs(color = "") +
  xlab("Emission (nm)") +
  ylab("Fluorescence (A.U)") +
  scale_color_manual(values = c("black", "#D55E00", "#0072B2")) +
  theme(legend.key = element_blank()) +
  theme(legend.justification = c(1, 1), legend.position = c(0.9, 0.9))

## ----fig3, echo = FALSE, message=FALSE, dev='png', fig.align='center', fig.width=7, fig.height=5, fig.cap = 'Surface plot of an EEM with first order of Raman and Rayleigh scattering removed. Missing values (`NA`) have been placed in both diagonals using a bandwidth of 10 nm.', results="hide"----

file <- system.file("extdata/cary/scans_day_1/", "sample3.csv", package = "eemR")

eem <- eem_read(file, recursive = FALSE)

eem <- eem_remove_scattering(eem, "rayleigh", 1, 10)
eem <- eem_remove_scattering(eem, "raman", 1, 10)

data("absorbance")

eem_scatter_removed <- eem_inner_filter_effect(eem, absorbance, 1)

file <- "/media/persican/Philippe Massicotte/Phil/articles/eemR/jss/figures/ife.rds"
ife <- readRDS(file)


jet.colors <- colorRampPalette(c("#00007F",
                                 "blue",
                                 "#007FFF",
                                 "cyan",
                                 "#7FFF7F",
                                 "yellow",
                                 "#FF7F00",
                                 "red",
                                 "#7F0000"))

par(mfrow = c(2,2), mar = c(4, 4, 1, 2) + 0.1, oma = c(0,0,0,2))

plot3D::image2D(y = eem[[1]]$em,
                x = eem[[1]]$ex,
                z = t(eem[[1]]$x),
                xlab = "Excitation (nm.)",
                ylab = "Emission (nm.)",
                legend.lab = "IFE correction factors",
                col = jet.colors(255))

legend("topleft",
       expression(bold("A")),
       text.font = 1,
       cex = 1.5,
       bty = "n",
       adj = c(1.5, 0.25))


plot(absorbance$wavelength,
     absorbance$sample3,
     type = "l",
     lwd = 2,
     col = "red",
     xlab = "Wavelength (nm)",
     ylab = "Absorbance")

legend("topleft",
       expression(bold("B")),
       text.font = 1,
       cex = 1.5,
       bty = "n",
       adj = c(1.5, 0.25))

plot3D::image2D(y = eem[[1]]$em,
                x = eem[[1]]$ex,
                z = t(ife),
                xlab = "Excitation (nm.)",
                ylab = "Emission (nm.)",
                legend.lab = "IFE correction factors",
                col = rev(jet.colors(255)))

legend("topleft",
       expression(bold("C")),
       text.font = 1,
       cex = 1.5,
       bty = "n",
       adj = c(1.5, 0.25))

plot3D::image2D(y = eem_scatter_removed[[1]]$em,
                x = eem_scatter_removed[[1]]$ex,
                z = t(eem_scatter_removed[[1]]$x),
                xlab = "Excitation (nm.)",
                ylab = "Emission (nm.)",
                legend.lab = "IFE correction factors",
                col = jet.colors(255))

legend("topleft",
       expression(bold("D")),
       text.font = 1,
       cex = 1.5,
       bty = "n",
       adj = c(1.5, 0.25))

## ----fig4, echo = FALSE, warning=FALSE, dev="png", fig.width = 7, fig.align="center", fig.cap = "IFE correction process. Panel (A) shows an uncorrected EEM (the color bar is the florescence intensity in A.U.). Panel (B) is the corresponding absorbance spectra measured on the same sample. Panel (C) shows the IFE correction factors corresponding to the values of the denominator in Equation (@ife) with values close to 1 indicating less pronounced correction. Panel (D) shows the corrected sample (the color bar is the fluorescence intensity in A.U.)."----

folder <- "/media/persican/Philippe Massicotte/Phil/Doctorat/PARAFAC/PARAFAC Files/Raw Data/Daniel/CDOM/3D/BlancNano5.csv"

eem <- eem_read(folder)

raman <- eem[[1]]$x[, which(eem[[1]]$ex == 350)]

df <- data.frame(em = eem[[1]]$em, raman)

shade <- rbind(c(371,0), subset(df, em >= 371 & em <= 428), c(df[nrow(df), "X"], 0))

ggplot(df, aes(x = em, y = raman)) +
  geom_line(size = 0.2) +
  xlab("Emission (nm)") +
  ylab("Fluorescence (A.U)") +
  geom_polygon(data = shade, aes(em, raman)) +
  xlim(300, 450) +
  annotate("text", x = 400, y = 1, label = "A[rp]", size = 2, parse = TRUE) +
  annotate("text",
           x = 390,
           y = 7.5,
           label = "First order Rayleigh scattering",
           size = 2)

## ------------------------------------------------------------------------
library(eemR)
ls("package:eemR")

## ------------------------------------------------------------------------
file <- system.file("extdata/cary/scans_day_1", package = "eemR")
eems <- eem_read(file)

## ------------------------------------------------------------------------
summary(eems)

## ---- eval = FALSE-------------------------------------------------------
#  plot(eems, which = 3)

## ------------------------------------------------------------------------
file <- system.file("extdata/cary/scans_day_1", "nano.csv", package = "eemR")
blank <- eem_read(file)

eems <- eem_remove_blank(eems, blank)

## ----fig5, tidy=FALSE, fig.height=4, fig.width=5, warning=FALSE, dev='png', fig.align='center', fig.cap = 'Fluorescence profile of a pure water sample at excitation 350 nm between 300 and 450 nm emission. The area of the Raman peak is identified by the shaded polygon and is calculated using Equation (@arp).'----

eems <- eem_remove_scattering(eems, "rayleigh", 1, 10) %>%
  eem_remove_scattering("raman", 1, 10)

plot(eems, which = 3)

## ------------------------------------------------------------------------
data("absorbance")
head(absorbance)

## ------------------------------------------------------------------------
eem_names(eems)

## ---- fig.keep="none"----------------------------------------------------
eems <- eem_inner_filter_effect(eem = eems,
                                absorbance = absorbance,
                                pathlength = 1)

plot(eems, which = 3)

## ---- fig.keep="none"----------------------------------------------------
eems <- eem_raman_normalisation(eems, blank)

plot(eems, which = 3)

## ------------------------------------------------------------------------
summary(eems)

## ---- eval=FALSE---------------------------------------------------------
#  eem_export_matlab("myfile.mat", eems)

## ------------------------------------------------------------------------
file <- system.file("extdata/cary/scans_day_1", package = "eemR")
eems <- eem_read(file)

eem_coble_peaks(eems, verbose = FALSE)

## ------------------------------------------------------------------------
eem_fluorescence_index(eems, verbose = FALSE)

eem_humification_index(eems, verbose = FALSE)

eem_biological_index(eems, verbose = FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  library(magrittr)
#  
#  file <- system.file("extdata/cary/scans_day_1/", package = "eemR")
#  file %>%
#    eem_read(recursive = TRUE) %>%
#    eem_raman_normalisation() %>%
#    eem_remove_scattering(type = "raman", order = 1, width = 10) %>%
#    eem_remove_scattering(type = "rayleigh", order = 1, width = 10) %>%
#    plot(2)

