get.chrom.num <- function(file.name, prefix) {

# Given a list of file names, and their common prefix, extracts the chromosome
# number that immediately follows the prefix. The number can be 1 or 2 digits.
# Returns the array of chromosome numbers; or NA if no number was obtained.

# Extract the 2 letters containing chromosome number:
#   expected to be right after <prefix>, if prefix="any.chr.":
#   examples: "any.chr.12.whatever.txt" = "12"
#             "any.chr.3stuffy.txt" = "3s"

len.p <- nchar(prefix)
two.letters <- substr(file.name, len.p+1, len.p+2)

# Next try to convert the chromosome number into numeric.
# Disable warnings - if some can not be converted, then
# it is a one-digit value, so remove the second letter
# symbol and convert again.
oldwarn <- getOption("warn")
options(warn=-1) # disable warning messages to avoid showing user expected NAs.
case.numeric <- as.numeric(two.letters) # convert to numeric 2-digit values.
options(warn=oldwarn) # enable warning options back to their original state

case.tmp <- file.name[is.na(case.numeric)] # extract all values that could not be converted to numeric
case.onedigit <- substr(case.tmp, len.p+1, len.p+1) # take the 1st letter from each (hopefully the 1-digit numeric)
case.numeric[is.na(case.numeric)] <- as.numeric(case.onedigit) # replace the NAs with converted numeric 1-digits

return(case.numeric)

}
