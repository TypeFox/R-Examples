### R code from vignette source 'msarc.Rnw'

###################################################
### code chunk number 1: msarc.Rnw:69-88
###################################################
library(msarc)
# load data
data(sample_df,package="msarc")
data(control_df,package="msarc")
# create the msarc objects
sample <- msarc(sample_df)
control <- msarc(control_df)
# subtract the control (i.e. remove UniProt IDs from sample that are also in
# control)
sample <- msarc.subtract(sample,control)
# generate the list of GO terms
sample <- msarc.findGOterms(sample,minCount=5)
# get the list as a data frame
term_df <- msarc.getTerms(sample)
term_df <- term_df[c("GO:0008092","GO:0017076","GO:0097159","GO:1901265",
                     "GO:0016787","GO:0017111"),]
sample <- msarc.filterTerms(sample,term_df)
# generate the plot
msarc.plotSVG(sample,file="thing.svg")


###################################################
### code chunk number 2: setup
###################################################
sessionInfo()


