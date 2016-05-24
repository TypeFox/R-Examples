naWindow <- function(x, union = F, zero = F){
  if(zero)
    invalid <- function(x){ is.na(x) | as.numeric(x) == 0 }
  else
    invalid <- is.na
  
  if(NCOL(x) == 1)
	valid <- !invalid(x)
  else {
	if(union) valid <- !apply(invalid(x), 1, all)
	else      valid <- !apply(invalid(x), 1, any)
  }
  tsx <- is.ts(x)
  if(tsx) x <- as.tis(x)
  validDates <- ti(x)[valid]
  if(length(validDates) > 0)
    z <- window(x, start = min(validDates), end = max(validDates))
  else z <- x
  if(tsx) as.ts(z)
  else z
}

