
####################################################
####  S3 CLASSES DEFINED AS S4 AND CLASS UNIONS ####
####################################################

## DNA SEQUENCES
setOldClass("DNAbin")
setOldClass("phylo")
setOldClass("multiPhylo")

## DATES
setOldClass("POSIXct")
setOldClass("Date")

## CONTACT NETWORKS
setOldClass("networkDynamic")
setOldClass("network")
setClassUnion("networkDynamicOrNetwork",c("network","networkDynamic"))

## ALLOW FOR SLOTS TO HAVE A TYPE, OR NULL
setClassUnion("characterOrNULL", c("character","NULL"))
setClassUnion("integerOrNULL", c("integer","NULL"))
setClassUnion("factorNULL", c("factor","NULL"))
setClassUnion("numericOrNULL", c("numeric","NULL"))
setClassUnion("matrixOrNULL", c("matrix","NULL"))
setClassUnion("listOrNULL", c("list","NULL")) # defined in adegenet
setClassUnion("DNAbinOrNULL", c("DNAbin","NULL"))
setClassUnion("POSIXctOrNULL", c("POSIXct","NULL"))
setClassUnion("networkDynamicOrNetworkOrNULL",c("networkDynamicOrNetwork","NULL"))
setClassUnion("DateOrNULL", c("Date", "NULL"))
setClassUnion("data.frameOrNULL", c("data.frame", "NULL"))
setClassUnion("multiPhyloOrNULL", c("multiPhylo", "NULL"))

