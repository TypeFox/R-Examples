
########### RAC875-Kukri Example
########### Example 1:

####### base model and diagnostics

wgpath <- system.file("extdata", package = "wgaim")

data(phenoRxK, package = "wgaim")

rkyld.asi <- asreml(yld ~ Type, random = ~ Genotype + Rep,
                   rcov = ~ ar1(Range):ar1(Row), data = phenoRxK)

summary(rkyld.asi)$varcomp
plot(rkyld.asi)
plot(variogram(rkyld.asi))

row.ind <- c(1,seq(4, 20, by = 3))
xyplot(resid(rkyld.asi) ~ Range | Row, data = phenoRxK, type = "b",
     panel = function(x, y, ...){
           panel.abline(h = 0, lty = 2)
           panel.xyplot(x, y, ...)},
     aspect = 4/5, layout = c(9,3), ylab = "Residuals",
     scales = list(x = list(at = row.ind,labels = phenoRxK$Row[row.ind])))

rkyld.asf <- asreml(yld ~ Type + lrange, random = ~ Genotype + Rep + Range,
                   rcov = ~ ar1(Range):ar1(Row), data = phenoRxK)

summary(rkyld.asf)$varcomp
plot(rkyld.asf)
plot(variogram(rkyld.asf))

xyplot(resid(rkyld.asf) ~ Row | Range, data = phenoRxK, type = "b",
     panel = function(x, y, ...){
           panel.abline(h = 0, lty = 2)
           panel.xyplot(x, y, ...)},
     aspect = 4/5, layout = c(9,3), ylab = "Residuals",
     scales = list(x = list(at = row.ind, labels = phenoRxK$Row[row.ind])))

######## genetic data

data(genoRxK, package = "wgaim")
genoRxK <- read.cross("csvr", file="genoRxK.csv", genotypes=c("AA","BB"),
        dir = wgpath, na.strings = c("-", "NA"))

names(genoRxK)
genoRxK$pheno[["Genotype"]][1:18]
summary(genoRxK)

names(genoRxK$geno$"3D")
genoRxK$geno$"3D"$data[200:208,1:8]
genoRxK$geno$"3D"$map

genoRxK <- cross2int(genoRxK, missgeno = "Mart", id = "Genotype", rem.mark = TRUE)
class(genoRxK)
names(genoRxK$geno$"3D")
genoRxK$geno$"3D"$map
genoRxK$geno$"3D"$imputed.data[200:208,1:8]
genoRxK$geno$"3D"$intval[200:208,1:6]

####### QTL analysis and breakout

rkyld.qtl0 <- wgaim(rkyld.asf, phenoData = phenoRxK, intervalObj = genoRxK,
merge.by = "Genotype", trace = TRUE, na.method.X = "include", breakout = 1,
exclusion.window = 0)

out.stat(rkyld.qtl0, genoRxK, iter = 1, stat= "blups")
out.stat(rkyld.qtl0, genoRxK, iter = 1, stat= "os")

asreml:::summary.asreml(rkyld.qtl0)$varcomp

####### QTL analysis and summary

rkyld.qtl1 <- wgaim(rkyld.asf, phenoData = phenoRxK, intervalObj = genoRxK,
merge.by = "Genotype", trace = TRUE, na.method.X = "include", exclusion.window = 0)

summary(rkyld.qtl1, genoRxK, LOD = FALSE)
out.stat(rkyld.qtl1, genoRxK, iter = 1:5, stat= "blups")

############ end script
