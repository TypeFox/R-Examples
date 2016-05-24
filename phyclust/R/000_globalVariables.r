### Do not delete this file and file name!!
### This file should be loaded before all other *.r files.

### This is to avoid the fale positive messages from R CMD check.
###   "no visible binding for global variable"
### Suggested by Prof Brian Ripley
### ?globalVariables

if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(
    ".Color",
    ".EMC",
    ".amino.acid",
    ".boundary.method",
    ".code.type",
    ".codon",
    ".edist.model",
    ".em.method",
    ".genetic.code",
    ".identifier",
    ".init.method",
    ".init.procedure",
    ".label.method",
    ".missing.code",
    ".nucleotide",
    ".se.model",
    ".snp",
    ".substitution.model"
    ))
}
