## ----loadseq-------------------------------------------------------------
library("sequences")
fastafilename <- dir(system.file(package="sequences", dir="extdata"),
                     full.name=TRUE,
                     pattern="fasta$")
fastafilename
myseq <- readFasta(fastafilename[1])
myseq

## ----printseq------------------------------------------------------------
print(myseq)

## ----transcribe----------------------------------------------------------
transcribe(myseq)

## ----gccount, dev='pdf', echo=TRUE,fig.width=6, fig.height=4-------------
barplot(gccount(seq(myseq)))

## ----sessioninfo, results='asis', echo=FALSE, cache=FALSE----------------
toLatex(sessionInfo())

