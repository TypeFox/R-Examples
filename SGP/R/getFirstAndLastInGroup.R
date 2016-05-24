`getFirstAndLastInGroup` <- 
function(DT) {
	VALID_CASE <- YEAR <- CONTENT_AREA <- FIRST_OBSERVATION <- LAST_OBSERVATION <- ID <- GRADE <- YEAR_WITHIN <- NULL
	setkey(DT, VALID_CASE, CONTENT_AREA, YEAR, GRADE, ID, YEAR_WITHIN)
	setkey(DT, VALID_CASE, CONTENT_AREA, YEAR, GRADE, ID)
	DT[DT[unique(DT["VALID_CASE"]),,mult="first", which=TRUE], FIRST_OBSERVATION:=1L]
	DT[DT[unique(DT["VALID_CASE"]),,mult="last", which=TRUE], LAST_OBSERVATION:=1L]
	return(DT)
} ### END getFirstAndLastInGroup
