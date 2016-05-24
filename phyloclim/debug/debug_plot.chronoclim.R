require(Roxalis)

setwd("/Users/Stoffi/R_files/hordeum/climatic_niche_Hordeum")
load("climtol_grass.RData")
pdfname <- "Fig3_climTol.pdf"

# extract layers
layers <- colnames(x$data)
ylab <- c(
"Mean annual temperature (°C)", 
"Precipitation of warmest quarter (cm)",
"Precipitation seasonality",
"Minimum temperature of coldest month (°C)")

# set tip colors according to distribution:
# ----------------------------------------
clades <- list(L1 = c("euclaston", "intercedens", "stenostachys", "pusillum"), L2 = c("comosum", "pubiflorum", "patagonicum"), L3 = c("cordobense", "muticum"), L4 = c("flexuosum", "chilense"), L5 = "californicum")

col <- c("orange", "red", "blue", "green", "purple")

ii <- 1

layer <- layers[6]
density <- TRUE

lwd = 2
xspace = c(-0.5, 0.1)

ylab <- ylab[ii]
