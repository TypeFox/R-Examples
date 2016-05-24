RandomDates <-
function(name, reps = 20) {

options (warn = -1)

inFileName <- paste0(name, ".xml")

inFile <- readLines(inFileName)

matchLines <- grep(pattern = "date value=", x = inFile, value = T)

if (length(matchLines) == 0) {stop(
"No dates found, check the integrity of input files and Beast version.")}  

numberDates <- length(grep("<sequence id", inFile))

matchLinesPosition <- grep(pattern = "date value=", x = inFile, value = F)

matchFileName <- grep(pattern = "fileName", x = inFile, value = T)
matchFileNamePosition <- grep(pattern = "fileName", x = inFile, value = F)

# Loop

for (i in 1 : reps){
  
newFile <- inFile

fileRep <- paste0(" fileName=\"Rep", i, ".")

newFile [matchFileNamePosition] <- gsub(" fileName=\"", fileRep, matchFileName) 

randLines <- sample(matchLines)
newFile[matchLinesPosition] <- randLines

outName <- paste0(name, "_Rep", i) 

out <- paste0(name, ".Rep", i, ".xml")

cat (newFile, file = out, sep = "\n")

}
cat ("Replicates done:", i,"\n")
}

