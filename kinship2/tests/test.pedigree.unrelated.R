#######################################
## Name: test.pedigree.unrelated.r
## Purpose: Test Suite for pedigree.unrelated
## Created: 3/29/2011
## Last Updated: 7/13/2011
## Author: Jason Sinnwell, MS
########################################

## examples from help file, available with
## > example(pedigree.unrelated)

require(kinship2)
#library(kinship2, lib.loc="~/Rdir/library")
#citation("kinship2")

data(sample.ped)


pedAll <- pedigree(sample.ped$id, sample.ped$father, sample.ped$mother,
                   sample.ped$sex, famid=sample.ped$ped,
       affected=cbind(sample.ped$affected, sample.ped$avail))
                   
ped1 <- pedAll['1']
ped2 <- pedAll['2']

## to see plot:
## plot.pedigree(ped1, align=FALSE)
set.seed(10)
id1 <- pedigree.unrelated(ped1, avail=ped1$affected[,2])

## some possible vectors
id1
# "109" "113" "133"
# "109" "110" "130"
# "109" "118" "141"

set.seed(10)
id2 <- pedigree.unrelated(ped2, avail=ped2$affected[,2])

## some possible vectors
id2
##[1] "203" "206"
##[1] "203" "213"
##[1] "203" "204"

