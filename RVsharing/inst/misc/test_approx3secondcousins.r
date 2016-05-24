source("approx3secondcousins.R")
load("extdata/ped.list.Rdata")
library(kinship2)

# Test of the approximation
sharing3secondcousins.obj = sharing3secondcousins()

# Get the ids of the three second cousins family (13th family in ped.list)
fam.id = ped.list[[13]]$id
fam.dadid = fam.momid = numeric(length(fam.id))
fam.dadid[ped.list[[13]]$findex>0] = ped.list[[13]]$id[ped.list[[13]]$findex]
fam.momid[ped.list[[13]]$findex>0] = ped.list[[13]]$id[ped.list[[13]]$mindex]

nf=8
phihat = 0.01

# Record the seed
seed.approx = .Random.seed
#.Random.seed = seed.approx

# Approximation
pop2.approx.sampest = approxsharing(nf,ord=5,phihat=phihat,sharing3secondcousins.obj,fam.id,fam.dadid,fam.momid,nsim=100000)