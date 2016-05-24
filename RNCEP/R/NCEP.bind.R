NCEP.bind <-
function(data.west, data.east){

## Bind the U and V arrays separately ##
whole.data <- abind(data.west, data.east, along=2)

## Rename the dimnames to be negative on the west side of the prime meridian ##
dimnames(whole.data)[[2]] <- ifelse(as.numeric(dimnames(whole.data)[[2]]) > 180, as.character(-1*(360 - as.numeric(dimnames(whole.data)[[2]]))), dimnames(whole.data)[[2]])

## Return the merged dataset ##
return(whole.data)
}

