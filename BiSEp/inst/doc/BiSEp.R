### R code from vignette source 'BiSEp.Snw'

###################################################
### code chunk number 1: BiSEp.Snw:18-23
###################################################

require(BiSEp)
data(INPUT_data)

INPUT_data[1:2,1:6]


###################################################
### code chunk number 2: packages
###################################################
BISEP_data <- BISEP(INPUT_data)
biIndex <- BISEP_data$BI
bisepIndex <- BISEP_data$BISEP


###################################################
### code chunk number 3: BiSEp.Snw:40-43
###################################################

biIndex[1:10,]
bisepIndex[1:10,]


###################################################
### code chunk number 4: fig1
###################################################
plot(density(INPUT_data["TUSC3",]), main="TUSC3 Density Distribution")


###################################################
### code chunk number 5: fig2
###################################################
plot(density(INPUT_data["MLH1",]), main="MLH1 Density Distribution")


###################################################
### code chunk number 6: packages
###################################################
plot(density(INPUT_data["MLH1",]), main="MLH1 Density Distribution")


###################################################
### code chunk number 7: BiSEp.Snw:72-73
###################################################
BIGEE_out <- BIGEE(BISEP_data, sampleType="cell_line")


###################################################
### code chunk number 8: BiSEp.Snw:78-79
###################################################
BIGEE_out[1:4,]


###################################################
### code chunk number 9: fig3
###################################################
expressionPlot(BISEP_data, gene1="SMARCA4", gene2="SMARCA1")


###################################################
### code chunk number 10: fig4
###################################################
expressionPlot(BISEP_data, gene1="MTAP", gene2="MLH1")


###################################################
### code chunk number 11: BiSEp.Snw:101-103
###################################################
data(MUT_data)
MUT_data[1:4,1:10]


###################################################
### code chunk number 12: BiSEp.Snw:106-107
###################################################
BEEMout <- BEEM(BISEP_data, mutData=MUT_data, sampleType="cell_line",  minMut=40)


###################################################
### code chunk number 13: BiSEp.Snw:112-113
###################################################
BEEMout


###################################################
### code chunk number 14: fig5
###################################################
waterfallPlot(BISEP_data, MUT_data, expressionGene="MICB", 
mutationGene="PBRM1")


###################################################
### code chunk number 15: fig6
###################################################
waterfallPlot(BISEP_data, MUT_data, expressionGene="BOK", 
mutationGene="BRCA2")


###################################################
### code chunk number 16: packages
###################################################
fOut <- FURE(BIGEE_out[1,], inputType="BIGEE")
frPairs <- fOut$funcRedundantPairs
allPairs <- fOut$allPairs


###################################################
### code chunk number 17: BiSEp.Snw:145-146
###################################################
allPairs[1,]


