TaxaOut <-
function(name) {

options (warn = -1)  

inFileName <- paste0(name, ".xml")

inFile <- readLines(inFileName)

endTaxaLine <- grep(pattern = "</taxa>", x = inFile, value = F)

taxonLine <- grep(pattern = "<taxon id=", x = inFile, value = T)

taxonLinePosition <- grep(pattern = "<taxon id=", x = inFile, value = F)

taxonLine <- unlist(strsplit(taxonLine, "\""))

taxon <- taxonLine[c (F, T, F)]

numberTaxa <- length(taxon)

if (length(taxon) == 0) {stop(
  "No date info found, check the integrity of files and Beast version.")}  

matchFileName <- grep(pattern = "fileName", x = inFile, value = T)

matchFileNamePosition <- grep(pattern = "fileName", x = inFile, value = F)

# Loop

for (i in 1 : numberTaxa) { 
  
newFile <- inFile

taxa <- taxon[i]

add1 <- grep(pattern = "</taxa>", x = inFile, value = F)

newFile [add1] <- paste0("\t</taxa>\n\n","\t<taxa id=\"leave_out\">\n",
  "\t\t<taxon idref=\"",taxa,"\"/>\n\t</taxa>")

add2 <- grep(pattern = "</treeModel>", x = inFile, value = F)

newFile [add2] <- paste0(
  "\n\t\t<!-- START Tip date sampling\t\t\t\t\t\t\t\t\t\t\t-->
  \n\t\t<leafHeight taxon=\"",taxa,"\">\n\t\t\t<parameter id=\"age(",
  taxa,"\"/>\n\t\t</leafHeight>
  \n\t\t<!-- END Tip date sampling\t\t\t\t\t\t\t\t\t\t\t\t-->
  \n\t</treeModel>","\n\n\t<!-- Taxon Sets\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t-->
  \t<tmrcaStatistic id=\"tmrca(leave_out)\" includeStem=\"false\">
  \t\t<mrca>\n\t\t\t<taxa idref=\"leave_out\"/>\n\t\t</mrca>
  \t\t<treeModel idref=\"treeModel\"/>\n\t</tmrcaStatistic>")

add34 <- grep(pattern = "</operators>", x = inFile, value = F)

newFile [add34] <- paste0(
  "\t\t<randomWalkOperator windowSize=\"1.0\" weight=\"1\">
  \t\t\t<parameter idref=\"age(",taxa,")\"/>\n\t\t</randomWalkOperator>
  \t</operators>")

add5 <- grep(pattern = "<prior id=\"prior\">", x = inFile, value = F)

newFile [add5] <- paste0(
  "\t\t\t<prior id=\"prior\">
  \t\t\t\t<uniformPrior lower=\"0.0\" upper=\"1.0E100\">
  \t\t\t\t\t<parameter idref=\"age(",taxa,")\"/>\n\t\t\t\t</uniformPrior>")

add6 <- grep(pattern = "<log id=\"fileLog\"", x = inFile, value = F)

keep6 <- inFile[add6 + 1]

newFile [add6 + 1] <- paste0(
  "\n\t\t\t<!-- START Tip date sampling\t\t\t\t\t\t\t-->
  \t\t\t<parameter idref=\"age(",taxa,")\"/>
  \t\t\t<!-- END Tip date sampling\t\t\t\t\t\t\t\t-->\n\n",keep6)

fileRep <- paste0(" fileName=\"Taxa", i, ".") # set the base name for Taxa

newFile [matchFileNamePosition] <- gsub(
  " fileName=\"", fileRep, matchFileName) 

outName <- paste0(name, "_Taxa", i) 

out <- paste0(name, ".Taxa", i, ".xml")

cat (newFile, file = out, sep = "\n")

cat ("Taxon", i, "processed \n")
}
}