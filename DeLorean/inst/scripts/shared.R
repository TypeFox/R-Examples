#
# Code that is shared amongst the vignettes.
#
library(ggplot2)
library(dplyr)

#
# Install extrafont and import fonts from system.
# We just need to do this once.
#
# install.packages("extrafont")
# library(extrafont)
# font_import()
# fonts()

# Load fonts for plotting
library(extrafont)
loadfonts()
loadfonts(device="postscript")

# font_import(pattern="[A/a]rial")
# loadfonts(device="win")

#
# Widths for saving figures
#
text.width <- 4.7  # LaTeX width in inches
golden.ratio <- 1.618  # Pleasing ratio
fig.width <- text.width
fig.height <- text.width / golden.ratio
plos.width <- 7.2 # LaTeX width for PLoS paper in inches
plos.height <- plos.width / golden.ratio
html5 <- list(width=1300, height=700)  # in pixels
html5$ratio <- with(html5, width / height)
slide.fig.width <- 7
bioinf.single.w <- 3.38  # Width of single bioinformatics column in inches
bioinf.double.w <- 7.00  # Width of double bioinformatics column in inches
bioinf.single.h <- bioinf.single.w / golden.ratio
bioinf.double.h <- bioinf.double.w / golden.ratio

#
# Theme for Bioinformatics
#
base.family <- "Helvetica"
bioinf.theme <- ggplot2::theme_grey(base_size=8, base_family=base.family)
bioinf.config <- list(
    width = bioinf.single.w,
    height = bioinf.single.h,
    theme = bioinf.theme)

#
# Theme for PLoS
#
base.family <- "Times"
# base.family <- "Arial"  # Arial may not be on all systems
# base.family <- "Liberation Sans"
plos.theme <- ggplot2::theme_classic(base_size=8, base_family=base.family)

#
# Fit model if not defined
#
if (! exists('fit.model')) fit.model <- TRUE

#
# Create data directory if not there
#
if( ! file.exists('Data')) dir.create('Data')
