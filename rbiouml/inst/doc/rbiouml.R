### R code from vignette source 'rbiouml.Rnw'

###################################################
### code chunk number 1: rbiouml.Rnw:35-37
###################################################
library(rbiouml)
biouml.login("ie.biouml.org")


###################################################
### code chunk number 2: rbiouml.Rnw:58-59
###################################################
biouml.ls("databases")


###################################################
### code chunk number 3: rbiouml.Rnw:62-63
###################################################
biouml.ls("data/Examples/Optimization/Data/Experiments")


###################################################
### code chunk number 4: rbiouml.Rnw:67-70
###################################################
x <- biouml.get("data/Examples/Optimization/Data/Experiments/exp_data_1")
class(x)
head(x)


###################################################
### code chunk number 5: rbiouml.Rnw:77-80
###################################################
x[,5] <- x[,3] + x[,4]
biouml.put("data/Collaboration/Demo/Data/Rtest/exp_data_1_sum", x)
biouml.ls("data/Collaboration/Demo/Data/Rtest")


###################################################
### code chunk number 6: rbiouml.Rnw:86-87
###################################################
summary( biouml.analysis.list() )


###################################################
### code chunk number 7: rbiouml.Rnw:91-92
###################################################
biouml.analysis.parameters("Filter table")


###################################################
### code chunk number 8: rbiouml.Rnw:95-100
###################################################
biouml.analysis("Filter table", list(
  inputPath="data/Examples/Optimization/Data/Experiments/exp_data_1",
  filterExpression="time < 40",
  outputPath="data/Collaboration/Demo/Data/Rtest/exp_data_1 filtered"
))


###################################################
### code chunk number 9: rbiouml.Rnw:107-108
###################################################
head( biouml.importers() )


###################################################
### code chunk number 10: rbiouml.Rnw:111-115
###################################################
hiv.genome <- system.file("extdata","HIV-1.fa", package="rbiouml")
output.folder <- "data/Collaboration/Demo/Data/Rtest"
biouml.import(hiv.genome, output.folder,  importer="Fasta format (*.fasta)")
biouml.ls(output.folder)


###################################################
### code chunk number 11: rbiouml.Rnw:119-123
###################################################
head( biouml.exporters() )
biouml.export("data/Collaboration/Demo/Data/Rtest/HIV-1",
  exporter="Fasta format (*.fasta)", target.file="HIV-1.fa")
file.exists("HIV-1.fa")


###################################################
### code chunk number 12: rbiouml.Rnw:128-129
###################################################
biouml.logout()


