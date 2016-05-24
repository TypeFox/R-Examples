## ----pphpc1, eval = FALSE------------------------------------------------
#  # Load library
#  library(micompr)
#  
#  # Output names
#  outputs <- c("$P^s$", "$P^w$", "$P^c$", "$\\mean{E}^s$",
#               "$\\overline{E}^w$", "$\\overline{C}$",
#               "$\\widetilde{A}$")
#  
#  # Outputs from the NetLogo implementation
#  dir_nl_ok <- paste0(dir_data, "nl_ok")
#  # Outputs from the Java implementation, first configuration
#  dir_jex_ok <- paste0(dir_data, "j_ex_ok")
#  # Outputs from the Java implementation, second configuration
#  dir_jex_noshuff <- paste0(dir_data, "j_ex_noshuff")
#  # Outputs from the Java implementation, third configuration
#  dir_jex_diff <- paste0(dir_data, "j_ex_diff")
#  
#  # Files for model size 400, parameter set 1
#  filez <- glob2rx("stats400v1*.txt")
#  
#  # Perform the three comparison cases
#  mic <- micomp(outputs,
#                ve_npcs = 0.75,
#                list(list(name = "I",
#                          folders = c(dir_nl_ok, dir_jex_ok),
#                          files = c(filez, filez),
#                          lvls = c("NLOK", "JEXOK")),
#                     list(name = "II",
#                          folders = c(dir_nl_ok, dir_jex_noshuff),
#                          files = c(filez, filez),
#                          lvls = c("NLOK", "JEXNS")),
#                     list(name = "III",
#                          folders = c(dir_nl_ok, dir_jex_diff),
#                          files = c(filez, filez),
#                          lvls = c("NLOK","JEXDIF"))),
#                concat = T)

## ----pphpc2, eval = FALSE------------------------------------------------
#  toLatex(mic,
#          booktabs = T,
#          data_show = c("npcs-1", "mnvp-1", "parp-1", "scoreplot"),
#          data_labels = c("$\\#$PCs", "MNV", "$t$-test", "PCS"),
#          col_width = T,
#          pvalf_params = list(minval = 1e-8, na_str = "*"),
#          label = "tab:pphpc",
#          caption = paste("Comparison of a NetLogo implementation of",
#                          "the PPHPC model against three configurations",
#                          "of a parallel Java implementation."))

## ----sunspot1, results = 'hide', warning = FALSE-------------------------
# Load library
library(micompr)

# Months in the 1749-1859 interval (110 years)
# Months in the 1902-2012 interval (110 years)
m <- sunspot.month[c(1:1320, 1837:3156)]
m <- matrix(m, nrow = 20)

# Factor vector, two levels:
# a) ten 11-year cycles from 1749 to 1859
# b) ten 11-year cycles from 1902 to 2012
groups <- factor(c(rep("A", 10), rep("B", 10)))

# Compare the two groups, use 9 PCs for MANOVA
cmp <- cmpoutput("SunSpots", 9, m, groups)

## ----sunspot2, results = 'markup', warning = FALSE-----------------------
cmp

## ----sunspot3, results = 'markup', warning = FALSE-----------------------
assumptions(cmp)

## ----sunspot4, fig.show = 'asis', fig.env = 'figure', fig.cap = 'Plots produced by sunspots example.'----
plot(cmp)

## ----saugeen1, results = 'markup', warning = FALSE-----------------------
# Load libraries
library(micompr)
library(deseasonalize)

# Unique years
years <- unique(sapply(rownames(SaugeenDay), substr, 1, 4))

# Number of days in each year
ndays <- sapply(years, function(x)
                sum(substr(rownames(SaugeenDay), 1, 4) == x))

# Indexes of last day in each year
lastdays <- cumsum(ndays)

# Prepare data for PCA
saugdata <- t(mapply(
  function(nd, ld) {
    rflows <- rep(NA, 366)
    rflows[1:nd] <- SaugeenDay[(ld - nd + 1):ld]
    # Discard last day in leap years
    rflows[1:365]
  },
  ndays, lastdays))

# Consider first 30 years and last 30 years (discard 5 years in between)
saugdata <- saugdata[c(1:30, 36:65), ]

# Factor vector, two levels: first 30 years and last 30 years
groups <- factor(c(rep("A", 30), rep("B", 30)))

# Compare
cmp <- cmpoutput("SaugeenFlow", 0.9, saugdata, groups)

## ----saugeen2, results = 'markup', warning = FALSE-----------------------
cmp

## ----saugeen3, results = 'markup', warning = FALSE-----------------------
assumptions(cmp)

## ----saugeen4, fig.show = 'asis', fig.env = 'figure', fig.cap = 'Plots produced by the Saugeen river flow example.'----
plot(cmp)

## ----derma1, eval = FALSE------------------------------------------------
#  # Load libraries
#  library(bmp)
#  library(micompr)
#  
#  # Image definitions
#  imgs <- dir(imgfolder)
#  nimgs <- length(imgs)
#  npixels <- 760 * 570
#  
#  # Specify image groups (Common nevi, atypical nevi,
#  # melanomas).
#  f <- read.table(grpsfile, row.names = 1)
#  grps <- f[order(row.names(f)), ]
#  
#  # Read images from disk
#  # Use different color channels as outputs, and also
#  # use a concatenated output
#  rimgs <- matrix(nrow = nimgs, ncol = npixels)
#  gimgs <- matrix(nrow = nimgs, ncol = npixels)
#  bimgs <- matrix(nrow = nimgs, ncol = npixels)
#  rgbimgs <- matrix(nrow = nimgs, ncol = npixels * 3)
#  
#  for (i in 1:nimgs) {
#  
#    cimg <- read.bmp(paste0(imgfolder, imgs[i]))
#    rimgs[i, ] <- c(cimg[ , , 1])
#    gimgs[i, ] <- c(cimg[ , , 2])
#    bimgs[i, ] <- c(cimg[ , , 3])
#    rgbimgs[i, ] <- c(cimg[ , , 1], cimg[ , , 2], cimg[ , , 3])
#  
#  }
#  
#  # Perform multivariate independent comparison of images
#  mic <-
#    micomp(outputs = c("R", "G", "B", "RGB"),
#           ve_npcs = 0.9,
#           comps = list(
#             list(name = "1v2",
#                  grpout = list(
#                    data = list(R = rimgs[grps != 3, ],
#                                G = gimgs[grps != 3, ],
#                                B = bimgs[grps != 3, ],
#                                RGB = rgbimgs[grps != 3, ]),
#                    obs_lvls = factor(grps[grps != 3]))),
#             list(name = "1v3",
#                  grpout = list(
#                    data = list(R = rimgs[grps != 2, ],
#                                G = gimgs[grps != 2, ],
#                                B = bimgs[grps != 2, ],
#                                RGB = rgbimgs[grps != 2, ]),
#                    obs_lvls = factor(grps[grps != 2]))),
#             list(name = "2v3",
#                  grpout = list(
#                    data = list(R = rimgs[grps != 1, ],
#                                G = gimgs[grps != 1, ],
#                                B = bimgs[grps != 1, ],
#                                RGB = rgbimgs[grps != 1, ]),
#                    obs_lvls = factor(grps[grps != 1])))))

## ----derma2, eval = FALSE------------------------------------------------
#  toLatex(mic,
#          booktabs = T,
#          data_show = c("parp-1", "nparp-1", "scoreplot"),
#          data_labels = c("$t$-test", "$U$ test", "PCS"),
#          pvalf_params = list(minval = 1e-8, na_str = "*"),
#          label = "tab:ph2",
#          caption = paste("Comparison of PH$^2$ dataset images",
#                          "grouped by lesion type."))

