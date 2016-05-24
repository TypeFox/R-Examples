### R code from vignette source 'sc-intro.Rnw'

###################################################
### code chunk number 1: sc-intro.Rnw:55-70
###################################################
dir.create('Figures', showWarnings=FALSE)
capture.output(source('Code/ex-compare-samplers.Rfrag'),
  file='Figures/ex-compare-samplers.txt')
capture.output(source('Code/ex-compare-results.Rfrag'),
  file='Figures/ex-compare-results.txt')

pdf('Figures/ex-comparison-plot.pdf', height=4, width=4)
print(comparison.plot(sampler.comparison, base_size=8))
dev.off()

source('Code/ex-metropolis.Rfrag')
source('Code/ex-beta.Rfrag')
source('Code/ex-final.Rfrag')
capture.output(source('Code/ex-final-compare.Rfrag'),
  file='Figures/ex-final-compare.txt')


