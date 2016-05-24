require(Roxalis)

ENM <- "/Users/Stoffi/R_files/maxent/output_diversity_TG_RST"
ENM <- list.files(ENM, pattern = ".asc", full.names = TRUE)
ENM <- ENM[-grep("clamping", ENM)]
ENM <- ENM[grep("virgosa|gigantea", ENM)]

gig <- import.asc(ENM[1], type = "numeric")
vir <- import.asc(ENM[2], type = "numeric")

# x: a pno profile
x <- read.table("/Users/Stoffi/R_files/hordeum/climatic_niche_Hordeum/pno_grass/01Annual_Mean_Temperature.txt")

# y: a vector of filenames
setwd("/Users/Stoffi/R_files/oxalis/climatic_niche_Alpinae")
load("alpinaetree.RData")
phy <- target
ENM <- "/Users/Stoffi/R_files/maxent/output_diversity_TG_RST"
ENM <- list.files(ENM, pattern = ".asc", full.names = TRUE)
ENM <- ENM[-grep("clamping", ENM)]
tips <- paste(phy$tip.label, collapse = "|")
y <- ENM[grep(tips, ENM)]

# z: a list of asc objects
z <- lapply(y, import.asc, type = "numeric")
names(z) <- phy$tip.label

