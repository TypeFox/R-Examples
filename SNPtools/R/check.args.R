check.args <-
function(chr, start, end, strains, available.strains) {

  # chr cannot be empty.
  if(missing(chr)) {
    stop(paste("chr cannot be empty. Please enter at least one chromosome."))
  } # if(is.null(chr))

  # start cannot be empty.
  if(missing(start)) {
    stop(paste("start cannot be empty. Please enter at least one start value."))
  } # if(is.null(start))

  # end cannot be empty.
  if(missing(end)) {
    stop(paste("end cannot be empty. Please enter at least one end value."))
  } # if(is.null(end))

  # The length of chr, start and end must be the same.
  if(length(chr) != length(start) | length(chr) != length(end)) {
    stop(paste("chr, start and end are not the same lengths. All three",
         "must be the same lenth. len(chr) =", length(chr), "; len(start) =",
         length(start), "; len(end) =", length(end), "."))
  } # if(length(chr) != length(start) |...

  # All start values must be less than or equal the end values.
  if(!all(start <= end)) {
    stop(paste("The start values must all be less than or equal to the end",
         "value."))
  } # if(end < start)

  # Requested strains in the list of available strains.
  if(!all(strains %in% available.strains)) {
    stop(paste("Strains not found in file:", 
         paste(strains[!strains %in% available.strains], collapse = ",")))
  } # if(!all(strains %in% available.strains))

}
