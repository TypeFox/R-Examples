ready <-
function(){setClass("Abstracts", representation(Journal = "character", Abstract = "character",PMID = "numeric"), where = ".GlobalEnv"); 
setClass("HGNC", representation(HGNCID = "character", ApprovedSymbol = "character", ApprovedName = "character", Status = "character", PreviousSymbols = "character", Aliases = "character", Chromosome = "character", AccessionNumbers = "character", RefSeqIDs = "character"), where = ".GlobalEnv")}

