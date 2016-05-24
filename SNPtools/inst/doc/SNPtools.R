### R code from vignette source 'SNPtools.Rnw'

###################################################
### code chunk number 1: SNPtools.Rnw:51-52
###################################################
library(SNPtools)


###################################################
### code chunk number 2: SNPtools.Rnw:100-102
###################################################
snp.file = "http://cgd.jax.org/tools/SNPtools/Build38/sanger.snps.NCBI38.txt.gz" 
mgi.file = "http://cgd.jax.org/tools/SNPtools/MGI/MGI.20130305.sorted.txt.gz"


###################################################
### code chunk number 3: SNPtools.Rnw:111-114
###################################################
available.strains = get.strains(snp.file)
strains = available.strains[c(4, 2, 8,  15, 16, 10, 17, 18)]
strains


###################################################
### code chunk number 4: SNPtools.Rnw:120-122
###################################################
data(qtl)
head(qtl)


###################################################
### code chunk number 5: SNPtools.Rnw:131-134
###################################################
cat.snps = variant.plot(var.file = snp.file, mgi.file = mgi.file,
           chr = 7, start = 103.3, end = 104.3, strains = strains,
           pattern = strains[c(3:5)], qtl = qtl)


###################################################
### code chunk number 6: SNPtools.Rnw:142-143
###################################################
head(cat.snps)


###################################################
### code chunk number 7: SNPtools.Rnw:161-166
###################################################
available.strains = get.strains(snp.file)
strains = available.strains[c(4, 2, 8,  15, 16, 10, 17, 18)]
snps = get.variants(chr = 7, start = 103.3, end = 104.3,
       strains = strains, type = "snp")
nrow(snps)


###################################################
### code chunk number 8: SNPtools.Rnw:171-173
###################################################
nsnps = convert.variants.to.numeric(snps)
snp.plot(nsnps, ref = "C57BL/6J", cluster = T)


###################################################
### code chunk number 9: SNPtools.Rnw:178-181
###################################################
strain.subset = c("C57BL/6J", "NOD/ShiLtJ", "NZO/HlLtJ")
snp.subset = get.pattern.variants(snps, strain.subset)
nrow(snp.subset)


###################################################
### code chunk number 10: SNPtools.Rnw:187-189
###################################################
snp.plot(nsnps, ref = "C57BL/6J", cluster = T, highlight = strain.subset,
         pattern = snp.subset)


###################################################
### code chunk number 11: SNPtools.Rnw:195-197
###################################################
data(qtl)
head(qtl)


###################################################
### code chunk number 12: SNPtools.Rnw:201-203
###################################################
snp.plot(nsnps, ref = "C57BL/6J", cluster = T, highlight = strain.subset,
         pattern = snp.subset, qtl = qtl)


###################################################
### code chunk number 13: SNPtools.Rnw:214-220
###################################################
mgi = get.mgi.features(chr = 7, start = 103, end = 105,
      source = "MGI", type = "gene")
mgi = mgi[-grep("^Gm", mgi$Name),]
snp.plot(nsnps, ref = "C57BL/6J", cluster = T,highlight =
         strain.subset, pattern = snp.subset, qtl = qtl, 
         mgi = mgi)


###################################################
### code chunk number 14: SNPtools.Rnw:230-234
###################################################
snps = get.variants(chr = 7, start = 103.3, end = 104.3,
       strains = strains, type = "snp")
snp.type = categorize.variants(snps)
head(snp.type[grep("Hbb", snp.type$symbol),])


###################################################
### code chunk number 15: SNPtools.Rnw:242-245
###################################################
indels = get.variants(file = "http://cgd.jax.org/tools/SNPtools/Build38/cc.indels.NCBI38.txt.gz",
       chr = 7, start = 103.3, end = 104.3, type = "indel")
head(indels)


###################################################
### code chunk number 16: SNPtools.Rnw:253-258
###################################################
mgi = get.mgi.features(chr = 7, start = 103.3, end = 104.3,
      source = "MGI", type = "gene")
col = rep(1, nrow(mgi))
col[grep("Hbb", mgi$Name)] = 2
gp = gene.plot(mgi, col = col)


