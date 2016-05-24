traverseTreeHelp <-
function(data, place, treeLvl=1, maxTaxaDepth, split){
REMOVELBL <- TRUE  #set to false to show all labels, otherwise only non-zero labels will be shown
childStr <- ""
myName <- rownames(data)[place]
myVal <- data[place, 1]
stopTree <- FALSE

if(treeLvl >= maxTaxaDepth)
stopTree <- TRUE

data <- data[-place, , drop=FALSE] 
childRows <- grep(paste(myName, split, sep=""), rownames(data), fixed=TRUE)

noSiblings <- TRUE
if(length(childRows) != 0 && !stopTree){ #we know we have children
newData <- data[childRows,, drop=FALSE]

rownames(newData) <- substring(rownames(newData), nchar(myName)+2, nchar(rownames(newData)))
nameSplit <- strsplit(rownames(newData), split, fixed=TRUE)

for(t in 1:nrow(newData)){ 
if(length(nameSplit[[t]]) == 1){ 
temp_hstr <- traverseTreeHelp(newData, t, treeLvl+1, maxTaxaDepth, split)
if(childStr != ""){
childStr <- paste(childStr, ",", temp_hstr, sep="")
}else{
childStr <- temp_hstr
}

noSiblings <- FALSE
}
}

if(myVal == 0){ #set up string to be returned
retStr <- paste("(", childStr, "):", myVal,  sep="")
}else{
retStr <- paste("(", childStr, ")", myName, ":", myVal,  sep="")
}

if(noSiblings) #adds an extra 0 value branch if there are no siblings
retStr <- paste(retStr, ",:0.0", sep="")
}else{ #no childern found
if(myVal == 0 && REMOVELBL == TRUE){
retStr <- paste(":", myVal, sep="")
}else{
retStr <- paste(myName, ":", myVal, sep="")
}

retStr <- paste(retStr, ",:0.0", sep="")
}

return(retStr)
}
