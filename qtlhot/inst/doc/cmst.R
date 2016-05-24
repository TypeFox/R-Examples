### R code from vignette source 'cmst.Rnw'

###################################################
### code chunk number 1: cmst.Rnw:112-113
###################################################
library(qtlhot)


###################################################
### code chunk number 2: cmst.Rnw:118-130
###################################################
set.seed(987654321)
CMSTCross <- SimCrossCausal(n.ind = 100,
                        len = rep(100, 3),
                        n.mar = 101,
                        beta = rep(0.5, 2),
                        add.eff = 1,
                        dom.eff = 0,
                        sig2.1 = 0.4,
                        sig2.2 = 0.1,
                        eq.spacing = FALSE,
                        cross.type = "bc",
                        normalize = TRUE)


###################################################
### code chunk number 3: cmst.Rnw:135-136
###################################################
CMSTCross <- calc.genoprob(CMSTCross, step = 1)


###################################################
### code chunk number 4: cmst.Rnw:141-145
###################################################
Scan <- scanone(CMSTCross, pheno.col = 1 : 3, method = "hk")
summary(Scan[, c(1, 2, 3)], thr = 3)
summary(Scan[, c(1, 2, 4)], thr = 3)
summary(Scan[, c(1, 2, 5)], thr = 3)


###################################################
### code chunk number 5: lodprofiles
###################################################
plot(Scan, lodcolumn = 1 : 3, ylab = "LOD")


###################################################
### code chunk number 6: cmst.Rnw:160-170
###################################################
commqtls <- GetCommonQtls(CMSTCross,
                          pheno1 = "y1",
                          pheno2 = "y3",
                          thr = 3,
                          peak.dist = 5,
                          addcov1 = NULL,
                          addcov2 = NULL,
                          intcov1 = NULL,
                          intcov2 = NULL)
commqtls


###################################################
### code chunk number 7: cmst.Rnw:175-187
###################################################
nms <- names(CMSTCross$pheno)
out1 <- CMSTtests(CMSTCross,
                  pheno1 = nms[1],
                  pheno2 = nms[2],
                  Q.chr = 1,
                  Q.pos = 55,
                  addcov1 = NULL,
                  addcov2 = NULL,
                  intcov1 = NULL,
                  intcov2 = NULL,
                  method = "all",
                  penalty = "both")


###################################################
### code chunk number 8: cmst.Rnw:192-193
###################################################
out1[1:3]


###################################################
### code chunk number 9: cmst.Rnw:198-199
###################################################
out1[4]


###################################################
### code chunk number 10: cmst.Rnw:204-205
###################################################
out1[5]


###################################################
### code chunk number 11: cmst.Rnw:210-211
###################################################
out1[6]


###################################################
### code chunk number 12: cmst.Rnw:216-217
###################################################
out1[7]


###################################################
### code chunk number 13: cmst.Rnw:222-223
###################################################
out1[8]


###################################################
### code chunk number 14: cmst.Rnw:228-229
###################################################
out1[9]


###################################################
### code chunk number 15: cmst.Rnw:236-237
###################################################
out1[10:12]


###################################################
### code chunk number 16: cmst.Rnw:242-243
###################################################
out1[13:17]


###################################################
### code chunk number 17: cmst.Rnw:248-260
###################################################
out2 <- CMSTtests(CMSTCross,
                  pheno1 = nms[1],
                  pheno2 = nms[-1],
                  Q.chr = 1,
                  Q.pos = 55.5,
                  addcov1 = NULL,
                  addcov2 = NULL,
                  intcov1 = NULL,
                  intcov2 = NULL,
                  method = "all",
                  penalty = "both")
out2


###################################################
### code chunk number 18: cmst.Rnw:269-271
###################################################
CMSTscan <- scanone(CMSTCross, pheno.col = 1:3, method = "hk")
CMSThigh <- highlod(CMSTscan)


###################################################
### code chunk number 19: cmst.Rnw:276-282
###################################################
traits <- names(CMSTCross$pheno)
annot <- data.frame(name = traits, traits = traits, chr = rep(1, 3),
 Mb.pos = c(55,10,100))
annot$cM.pos <- annot$Mb.pos
annot
targets <- list(y1 = c("y2","y3"))


###################################################
### code chunk number 20: cmst.Rnw:287-293
###################################################
cand.reg <- GetCandReg(CMSThigh, annot, traits)
cand.reg
cis.cand.reg <- GetCisCandReg(CMSThigh, cand.reg)
cis.cand.reg
comap.targets <- GetCoMappingTraits(CMSThigh, cand.reg)
comap.targets


###################################################
### code chunk number 21: cmst.Rnw:298-308
###################################################
tests <- list()
for(k in seq(names(comap.targets))) {
  tests[[k]] <- FitAllTests(CMSTCross, pheno1 = names(comap.targets)[k],
                      pheno2 = comap.targets[[k]],
                      Q.chr = cand.reg[k, 4],
                      Q.pos = cand.reg[k, 5])
}
names(tests) <- names(comap.targets)
tests <- JoinTestOutputs(comap.targets, tests)
tests


###################################################
### code chunk number 22: cmst.Rnw:313-316
###################################################
PrecTpFpMatrix(alpha = seq(0.01, 0.10, by = 0.01),
  val.targets = targets, all.orfs = CMSThigh$names, tests = tests,
  cand.reg = cand.reg, cis.cand.reg = cis.cand.reg)


