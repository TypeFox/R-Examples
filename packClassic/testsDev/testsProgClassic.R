#################################################################
###                    testProgClassique.R                    ###
#################################################################
load("../data/dataAges.rda")
source("../R/progClassic.R")
publicA(1)
privateA(2)
.publicB(3)
.privateB(4)
publicC(5)
publicC(dataAges$age)
privateC(6)
