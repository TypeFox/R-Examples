## ----setup, echo=FALSE--------------------------------------------------------
options(tidy=TRUE, width=80)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
library(metaMix)
###Location of input files.
datapath <- system.file("extdata", package="metaMix")
blastOut.default<-file.path(datapath, "blastOut_default.tab")
read.table(blastOut.default, nrows=2, sep="\t")

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
read.lengths<-file.path(datapath, "read_lengths.tab")
read.weights<-file.path(datapath, "read_weights.tab")
taxon.file<-file.path(datapath, "gi_taxid_prot_example.dmp")

read.table(read.lengths, nrows=2, sep="\t")
read.table(read.weights, nrows=2, sep="\t")
read.table(taxon.file, nrows=2, sep="\t")

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
blastOut.custom<-file.path(datapath, "blastOut_custom.tab")
read.table(blastOut.custom, nrows=2, sep="\t")

## ----echo=TRUE----------------------------------------------------------------
  step1 <-generative.prob(blast.output.file = blastOut.custom,
                          contig.weight.file=read.weights,
                          blast.default=FALSE,
                          outDir=NULL)

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  step1 <-generative.prob(blast.output.file = blastOut.default,
#                          read.length.file=read.lengths,
#                          contig.weight.file=read.weights,
#                          gi.taxon.file = taxon.file,
#                          blast.default=TRUE,
#                          outDir=NULL)

## ----echo=TRUE----------------------------------------------------------------
###The resulting list consists of five elements
names(step1)

### The sparse matrix of generative probs. 
step1$pij.sparse.mat[1:5,c("374840",  "258", "unknown")]

### There are that many potential species in the sample:
nrow(step1$ordered.species)

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  step2 <- reduce.space(step1=step1)

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  step2 <- reduce.space(step1="/pathtoFile/step1.RData")

## ----echo=FALSE, eval=TRUE----------------------------------------------------
data(step2)

## ----echo=TRUE----------------------------------------------------------------
##These are the elements of the step2 list.
names(step2)

## After this approximating step, there are now that many potential species in
##the sample:
nrow(step2$ordered.species)

## And these are:
step2$ordered.species

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  step3<-parallel.temper(step2=step2)

## ----echo=FALSE, eval= TRUE---------------------------------------------------
data(step3)

## ----echo=TRUE----------------------------------------------------------------
##These are the elements of the step3 list.
names(step3)

## ----echo=TRUE----------------------------------------------------------------
## Steps MCMC took during some iterations.
step3$result$slave1$record[10:15,]

## ----echo=TRUE----------------------------------------------------------------
## Location of the taxonomy names file.
taxon.file<-file.path(datapath, "names_example.dmp")

step4<-bayes.model.aver(step2=step2,
                        step3=step3,
                        taxon.name.map=taxon.file)

## ----echo=TRUE----------------------------------------------------------------
##These are the elements of the step4 list.
names(step4)

##This is the species summary
print(step4$presentSpecies.allInfo)

## ----echo=FALSE, eval=TRUE, out.width='.6\\linewidth', out.height='.6\\linewidth', fig.align='center'----
PTastro<-file.path(datapath, "PT_plots.RData")
load(PTastro)
nIter<- length(PTresult$result$slave1$record[,'logL'])
plot(PTresult$result$slave1$record[(nIter/5):nIter,'logL'], type='l', col='dodgerblue', xlab='Last 80% of iterations', ylab='Log-likelihood', main='Parallel Tempering - Coldest Chain', lwd=1.5)


## ----echo=TRUE----------------------------------------------------------------
sessionInfo()

