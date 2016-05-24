uniqueID <- function(ped.infile){
##
## Checks if ped file contains unique individual IDs. Returns TRUE if ped file contains unique individual IDs, FALSE otherwise.
##
## ped.infile is a character string giving the name and path of the standard ped file to be modified.
##
#
## Extract ID column from ped file using DatABEL's extract_text_file_columns
.id <- extract_text_file_columns(file = ped.infile, whichcols = 2)
#
## Check if the vector of individual IDs contains duplicated values 
.unique.id <- !any(duplicated(.id))
#
return(.unique.id)
}
