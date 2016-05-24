ListTaxa <-
function(name){

options (warn = -1)

inFileName <- paste0(name, ".xml")

inFile <- readLines(inFileName)

endTaxaLine <- grep(pattern = "</taxa>", x = inFile, value = F)

taxonLine <- grep(pattern = "<taxon id=", x = inFile, value = T)

taxonLinePosition <- grep(pattern = "<taxon id=", x = inFile, value = F)

taxonLine <- unlist(strsplit(taxonLine, "\""))

taxon <- taxonLine[c(F, T, F)]

numberTaxa <- length(taxon)

if (length(taxon) == 0) {stop("No date info found, check input files.")}  
    
matchFileName <- grep(pattern = "fileName", x = inFile, value = T)

matchFileNamePosition <- grep(pattern = "fileName", x = inFile, value = F)

# Loop

for (i in 1 : numberTaxa) { 

cat ("Taxon", i, "is", taxon[i], "\n")
}
}