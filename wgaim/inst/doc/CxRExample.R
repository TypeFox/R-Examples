## 2012 wgaim 1.3+ script

############ Cascades x RAC875-2 Zinc Experiment
############ Example 1:

## exploratory analysis

data("phenoCxR", package = "wgaim")
xyplot(znconc ~ shoot, data = phenoCxR, type = "p", xlab = "Shoot Length", ylab = "Zinc Concentration")
zincm <- with(phenoCxR, aggregate(cbind(znconc, shoot), list(Block = Block, Type= Type), mean))
zincd <- cbind(val = c(zincm$znconc, zincm$shoot), trait = rep(c("Zn Conc", "Shoot"), each = 6), rbind(zincm[,1:2]), row.names = NULL)
barchart(Type~ val | trait, groups = Block, data = zincd, scales = list(relation = "free"), col = 4:5, xlab = "")

## base model

sh.fm <- asreml(shoot ~ Type, random = ~ Block + id, data = phenoCxR)
summary(sh.fm)$varcomp
wald(sh.fm)

## genetic data

# read in genetic data to an R/qtl "cross" object

data("genoCxR", package = "wgaim")
wgpath <- system.file("extdata", package = "wgaim")
read.csv(paste(wgpath, "/genoCxR.csv", sep = ""), header = FALSE)[1:10,1:10]
genoCxR <- read.cross("csvr", file="genoCxR.csv", genotypes=c("AA","BB"),
               dir = wgpath, na.strings = c("-", "NA"))

summary(genoCxR)
class(genoCxR)
names(genoCxR$geno)
names(genoCxR$geno$"1B")
genoCxR$geno$"1B"$data[1:8,1:8]
genoCxR$geno$"1B"$map

# convert "cross" object to a "interval" object for use in wgaim

genoCxR <- cross2int(genoCxR, missgeno="Mart", id = "id", rem.mark = FALSE)
class(genoCxR)
names(genoCxR$geno$"1B")
genoCxR$geno$"1B"$imputed.data[1:8,1:8]
genoCxR$geno$"1B"$intval[1:6,1:6]

## wgaim QTL analyses

# QTL interval analysis

sh.qtlI <- wgaim(sh.fm, phenoData = phenoCxR, intervalObj = genoCxR, merge.by = "id", na.method.X = "include", gen.type = "interval", method = "fixed", selection = "interval")

summary(sh.qtlI, genoCxR, LOD = FALSE)

# QTL marker analysis

sh.qtlM <- wgaim(sh.fm, phenoData = phenoCxR, intervalObj = genoCxR, merge.by = "id", na.method.X = "include", gen.type = "marker", method = "fixed", selection = "interval")

summary(sh.qtlM, genoCxR)

##########
##########
