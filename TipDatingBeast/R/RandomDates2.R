RandomDates2 <-
function(name, reps = 20) {

options (warn = -1)

inFileName <- paste0(name, ".xml")

inFile <- readLines(inFileName)

numberDates <- length(grep("<sequence id", inFile))

matchTraitDateLine <- grep(pattern = "dateTrait", x = inFile, value = F)

if (length(matchTraitDateLine) == 0) {stop(
 "No date info found, check the integrity of files and Beast version.")}  

datePositions <- seq(matchTraitDateLine + 1, matchTraitDateLine + numberDates)

dateLines <- inFile[datePositions]

date <- unlist(strsplit(dateLines, "="))

dateHap <- date[c(T, F)]

dateHap <- dateHap[1: numberDates]

dateValues <- date[c(F, T)]

lastLine <- length(grep("<taxa", dateValues))

dateValues <- na.omit(as.numeric (gsub("[^\\d]+", "", dateValues, perl = T)))

# Loop

for(i in 1 :reps) {

newFile <- inFile 

dateValues <- sample(dateValues)

newDate <- paste0("\t\t\t", dateHap, "=", dateValues)

newFile[datePositions] <- paste0(newDate, ",")

if(lastLine == 0){newFile[(matchTraitDateLine + numberDates)] <- 
  paste0 (newDate[numberDates])}

if(lastLine == 1){newFile[(matchTraitDateLine + numberDates)] <-
  paste0 (newDate[numberDates], "\t\t\t\t<taxa id=", date[numberDates* 2 + 1],
  "=", date[numberDates * 2 + 2])}

fileRep <- paste0(" fileName=\"Rep", i, ".")

newFile <- gsub(" fileName=\"", fileRep, newFile)
  
outName <- paste0(name, "_Rep", i) 

out <- paste0(name, ".Rep", i, ".xml")

cat (newFile, file = out, sep = "\n")

}
cat ("Replicates done:", i,"\n")
}
