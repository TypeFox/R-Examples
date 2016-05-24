getBreaks <- function (colNames, breaks=NULL) 
{
    popColumns <- grep(
        "^([[:space:]]+)?(m|f|male|female)(_|\\.)?[[:digit:]]+(_|-|to|plus|\\+)?[[:digit:]]*$", 
        colNames, value = TRUE, ignore.case = TRUE)
    ageGroups <- gsub("^([[:space:]]+)?(m|f|male|female)(\\.|_)?", "", popColumns, ignore.case = TRUE)
    ageGroups <- gsub("(\\+|plus)", "_Inf", ageGroups, ignore.case = TRUE)
    ageLower <- as.numeric(gsub("(_|-|to|\\.)([[:digit:]]+|Inf)$", "", 
        ageGroups))
    ageUpper <- as.numeric(gsub("^[[:digit:]]+(_|-|to)", "", ageGroups))
    currentbreaks <- sort(unique(c(ageLower, max(ageUpper))))

	# check upper limit of breaks 
	if(max(currentbreaks)<=max(ageLower))
		currentbreaks = c(currentbreaks, Inf)
	
    sex <- toupper(substr(popColumns, 1, 1))

    result = list(breaks = currentbreaks, age = ageLower, 
        sex = sex, oldNames = popColumns, newNames = paste(toupper(sex), 
            ageLower, sep = "."))

if(length(breaks) ) {
 # check that the breaks are valid

# make sure that the right end of the last bin is Inf
	if(breaks[length(breaks)] != Inf) {
	  if(breaks[length(breaks)] %in% currentbreaks) {
		 breaks = c(breaks, Inf)
	  } else {
		breaks[length(breaks)] = Inf
	  }
	}

# if the breaks don't include some younger age groups, keep those age groups as-is
result$breaks = c(currentbreaks[currentbreaks < min(breaks)], breaks)

# check that population breaks are nested within breaks
if (!all(breaks %in% currentbreaks)) {
	warning("population data doesn't appear to be nested within breaks, ignoring breaks")
}
}
result
}

