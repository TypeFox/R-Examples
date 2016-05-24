#3-6-2011 MRC-Epid JHZ

library(gap)
library(genetics)
apoe <- genotype(PD$apoe)
rs10506151 <- genotype(PD$rs10506151)
rs10784486 <- genotype(PD$rs10784486)
rs1365763 <- genotype(PD$rs1365763)
rs1388598 <- genotype(PD$rs1388598)
rs1491938 <- genotype(PD$rs1491938)
rs1491941 <- genotype(PD$rs1491941)
m770 <- genotype(PD$m770)
int4 <- genotype(PD$int4)
snca <- genotype(PD$snca)
apoe.a1 <- allele(apoe)[,1]
apoe.a2 <- allele(apoe)[,2]
rs10506151.a1 <- allele(rs10506151)[,1]
rs10506151.a2 <- allele(rs10506151)[,2]
rs10784486.a1 <- allele(rs10784486)[,1]
rs10784486.a2 <- allele(rs10784486)[,2]
rs1365763.a1 <- allele(rs1365763)[,1]
rs1365763.a2 <- allele(rs1365763)[,2]
rs1388598.a1 <- allele(rs1388598)[,1]
rs1388598.a2 <- allele(rs1388598)[,2]
rs1491938.a1 <- allele(rs1491938)[,1]
rs1491938.a2 <- allele(rs1491938)[,2]
rs1491941.a1 <- allele(rs1491941)[,1]
rs1491941.a2 <- allele(rs1491941)[,2]
m770.a1 <- allele(m770)[,1]
m770.a2 <- allele(m770)[,2]
int4.a1 <- allele(int4)[,1]
int4.a2 <- allele(int4)[,2]
snca.a1 <- allele(snca)[,1]
snca.a2 <- allele(snca)[,2]

library(haplo.stats)
geno.lrrk2 <- setupGeno(data.frame(
rs10506151.a1,
rs10506151.a2,
rs10784486.a1,
rs10784486.a2,
rs1365763.a1,
rs1365763.a2,
rs1388598.a1,
rs1388598.a2,
rs1491938.a1,
rs1491938.a2,
rs1491941.a1,
rs1491941.a2))

geno.snca <- setupGeno(data.frame(
m770.a1,
m770.a2,
int4.a1,
int4.a2,
snca.a1,
snca.a2))

label.snca <- c("m770", "int4", "snca")
label.lrrk2 <- c("rs10506151", "rs10784486", "rs1365763", "rs1388598", "rs1491938", "rs1491941")
gender <- as.factor(PD$sex)
data.lrrk2 <- with(PD,data.frame(geno=geno.lrrk2,gender,abc,pd,diag,aon,race,apoe2,apoe3,apoe4,apoe234))
data.snca <- with(PD,data.frame(geno=geno.snca,gender,abc,pd,diag,aon,race,apoe2,apoe3,apoe4,apoe234))
# lrrk2
apoe_lrrk2 <- haplo.glm(formula=pd~gender+apoe2*geno.lrrk2,family="binomial",data=data.lrrk2,locus.label=label.lrrk2)
apoe_lrrk2
apoe_lrrk2 <- haplo.glm(formula=pd~gender+apoe4*geno.lrrk2,family="binomial",data=data.lrrk2,locus.label=label.lrrk2)
apoe_lrrk2
apoe_lrrk2 <- haplo.glm(formula=pd~gender+apoe234*geno.lrrk2,family="binomial",data=data.lrrk2,locus.label=label.lrrk2)
apoe_lrrk2
# snca
apoe_snca <- haplo.glm(formula=pd~gender+apoe2*geno.snca,family="binomial",data=data.snca,locus.label=label.snca)
apoe_snca
apoe_snca <- haplo.glm(formula=pd~gender+apoe4*geno.snca,family="binomial",data=data.snca,locus.label=label.snca)
apoe_snca
apoe_snca <- haplo.glm(formula=pd~gender+apoe234*geno.snca,family="binomial",data=data.snca,locus.label=label.snca)
apoe_snca
