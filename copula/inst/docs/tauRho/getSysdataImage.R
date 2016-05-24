##########################################################################
## This script generates the sysdata.rda file used in the copula package
##########################################################################

## this part was done while working on KojYan2010 (IME)
## load images generated from trpsrho.R
load("claytonRho.rda") 
load("gumbelRho.rda")
## load image generated from trpstau.R
load("plackettTau.rda")

## load images generated from evtrps.R
## this part was done while working on GKNY 2011 (Bernoulli)
load("galambos.rda")
load("huslerReiss.rda")
load("tev.rda")

save(.claytonRhoNeg, .claytonRhoPos, ## claytonRhoFun, claytondRho,
     .gumbelRho, ## gumbelRhoFun, gumbeldRho,
     .plackettTau, ## plackettTauFun, plackettdTau,
     .galambosTau, .galambosRho,
     .huslerReissTau, .huslerReissRho,
     .tevTau, .tevRho,
     file = "sysdata.rda", compress=TRUE)

#################################################################################
## NOTE:
## gumbelRhoFun should not be defined in sysdata.rda.
## Otherwise .gumbelRho would not be found in R CMD check
#################################################################################
