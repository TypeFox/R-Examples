# ######### Simulate mQTL and vQTL modeled with scanonevar
# library(qtl)
# library(vqtl)
#
# B6.C58.cross <- readRDS('~/Dropbox (ValdarLab)/Aim1_for_local/vQTL_reanalysis/rds/B6_C58_cross_rem2ol_transformed_and_fa.RDS')
#
#
#
# a <- scanonevar(cross = B6.C58.cross,
#            mean.formula = formula('TOTDIST ~ sex + mean.QTL.add + mean.QTL.dom'),
#            var.formula = formula('~ sex + var.QTL.add + var.QTL.dom'),
#            chrs = c(19, 'X'))
#
# b <- scanone(cross = B6.C58.cross, chr = c(19, 'X'), pheno.col = 10)
# plot(a, b)
#
# set.seed(27599)
#
# data(fake.f2)
# fake.f2 <- calc.genoprob(fake.f2, step = 2)
#
# N = nind(fake.f2)
# fake.f2$pheno$sex <- rbinom(n = N, size = 1, prob = 0.5)
# fake.f2$pheno$age <- rnorm(n = N, mean = 10, sd = 1)
#
# col.sex <- rep('', nind(fake.f2))
# col.sex[fake.f2$pheno$sex == 0] <- 'red'
# col.sex[fake.f2$pheno$sex == 1] <- 'blue'
#
#
# # NULL PHENOTYPE
# fake.f2$pheno$phen1 <- rnorm(n = N, mean = 20, sd = 2)
#
# margin.plot(cross = fake.f2,
#             focal.phenotype.name = 'phenotype',
#             marginal.phen.names = list('sex', 'age'),
#             marginal.marker.names = 'DXM66',
#             pch = 19, col = col.sex,
#             subset = (col.sex == 'blue'))
#
#
# predictive.plot(cross = fake.f2,
#                 mean.formula = formula('phen1 ~ sex + age + mean.QTL.add + mean.QTL.dom'),
#                 var.formula = formula('~sex + age + var.QTL.add + var.QTL.dom'),
#                 marker.name = 'D3M17', phen.name = 'sex')
#
# # scanone
# scan1 <- scanone(cross = fake.f2, chr = c(17:19, 'X'), pheno.col = 'phen1')
# plot(x = scan1, bandcol = 'gray')
#
# # scanonevar
# varscan1 <- scanonevar(cross = fake.f2,
#                        mean.formula = formula('phen1 ~ sex + age  + D17M66 + mean.QTL.add + mean.QTL.dom'),
#                        var.formula = formula('~sex + age + var.QTL.add + var.QTL.dom'),
#                        chrs = c(17:19, 'X'))
# plot(x = varscan1, y = scan1, main = 'test')
#
#
# # # do permutations and convert to empirical p values
# # varscan1.perms <- scanonevar.perm(cross = fake.f2,
# #                                   mean.formula = formula('phen1 ~ sex + age + D17M66 + mean.QTL.add + mean.QTL.dom'),
# #                                   var.formula = formula('~sex + age + var.QTL.add + var.QTL.dom'),
# #                                   chrs = c(17:19, 'X'),
# #                                   n.perms = 20)
#
# # varscan1b <- scanonevar.to.p.values(scanonevar = varscan1, perm.scan.maxes = varscan1.perms)
#
# # plot(x = varscan1b, ylim = c(0, 3))
#
# # attr(varscan1b, 'units') <- 'lods'
# # plot(x = varscan1b, y = scan1)
#
# # MQTL ON CHROMOSOME 19
# marker2.name <- colnames(fake.f2$geno$`19`$data)[2]
# marker2.vals <- get.genotypes.by.marker.name(cross = fake.f2, marker.name = marker2.name, as.matrix = FALSE) - 2
# marker2.vals[is.na(marker2.vals)] <- 0
#
# fake.f2$pheno$phen2 <- rnorm(n = N, 25 + marker2.vals, 3)
#
# scan2 <- scanone(cross = fake.f2, chr = c(15:19, 'X'), pheno.col = 'phen2')
# varscan2 <- scanonevar(cross = fake.f2,
#                        mean.formula = formula('phen2 ~ sex + age + mean.QTL.add + mean.QTL.dom'),
#                        var.formula = formula('~sex + age + var.QTL.add + var.QTL.dom'),
#                        chrs = c(15:19, 'X'))
# summary(varscan2)
# plot(x = varscan2, y = scan2, chrs = c(15, 17, 'X'))
# plot(x = varscan2, y = scan2, chrs = 15)
#
#
# predictive.plot(cross = fake.f2,
#                 mean.formula = formula('phen2 ~ sex + age + mean.QTL.add + mean.QTL.dom'),
#                 var.formula = formula('~sex + age + var.QTL.add + var.QTL.dom'),
#                 marker.name = marker2.name,
#                 phen.name = 'sex')
#
#
#
#
#
#
# # an additive vQTL on chr 15
# marker3.name <- colnames(fake.f2$geno$`19`$data)[2]
# marker3.vals <- get.genotypes.by.marker.name(cross = fake.f2, marker.name = marker3.name, as.matrix = FALSE) - 2
# marker3.vals[is.na(marker3.vals)] <- 0
#
# fake.f2$pheno$phen3 <- rnorm(n = N, 25, exp(0.5*marker3.vals))
#
# varscan3 <- scanonevar(cross = fake.f2,
#                        mean.formula = formula('phen3 ~ sex + age + mean.QTL.add + mean.QTL.dom'),
#                        var.formula = formula('~sex + age + var.QTL.add + var.QTL.dom'),
#                        chrs = c(15:19, 'X'))
#
# summary(varscan3)
# plot(varscan3)
#
# predictive.plot(cross = fake.f2,
#                 mean.formula = formula('phen3 ~ age + sex*mean.QTL.add + sex*mean.QTL.dom'),
#                 var.formula = formula('~sex + age + sex*var.QTL.add + sex*var.QTL.dom'),
#                 marker.name = marker3.name,
#                 phen.name = 'sex')
#
# #
# # perms <- scanonevar.perm(cross = fake.f2,
# #                          mean.formula = formula('phen3 ~ sex + age + mean.QTL.add + mean.QTL.dom'),
# #                          var.formula = formula('~sex + age + var.QTL.add + var.QTL.dom'),
# #                          n.perms = 50,
# #                          chrs = c(15:19, 'X'))
#
# #
# # varscan3b <- convert.scanonevar.to.empirical.ps(scan = varscan3,
# #                                                 null.scan.maxes = perms)
#
# # if the effects are really strong, the empricical p value will underflow R's float
# # and it's log will be -Inf, so we can't plot....maybe should replace 0's with .Machine$double.eps?
# # not likely to come up in real work, so I'll leave it as an error-throwing case for now
# # plot(varscan3b)
#
#
#
# fake.f2$pheno$numSomething <- rpois(n = nind(fake.f2), lambda = 8)
#
# # scanonevar with poisson regression
# varscan1b <- scanonevar(cross = fake.f2,
#                         mean.formula = formula('numSomething ~ sex + age  + D17M66 + mean.QTL.add + mean.QTL.dom'),
#                         var.formula = formula('~sex + age + var.QTL.add + var.QTL.dom'),
#                         chrs = c(17:19, 'X'),
#                         family = 'poisson')
# plot(x = varscan1b, y = scan1, main = 'test')
