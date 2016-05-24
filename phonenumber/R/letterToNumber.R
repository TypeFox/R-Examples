#' 
#' Convert letters to numbers as on a telephone's keypad
#' 
#' Take a character vector and convert it to the equivalent number sequence
#' from a telephone's key pad
#' 
#' @param value An input value as a character vector with one element (a string)
#' @param qz Whether to assign q and z to zero (qz = 0) or not (any other value)
#' 
#' @export
#' 
#' @examples
#' # Convert an alphabetic string can be converted directly (with 
#' # non-alphanumeric characters replaced by dashes)
#' 
#' letterToNumber("R functions") # returns "7-386284667"
#' 
#' # Of course, vectors containing strings can also be converted
#' 
#' string <- "Phone Number"
#' letterToNumber(string) # returns "74663-686237"
#' 
#' # Alphanumeric strings can also be converted with numbers being returned as
#' # is
#' 
#' letterToNumber("Jenny's number is 867-5309") # returns "53669-7-686237-47-867-5309"
#' 
#' # Specifying qz = 0 maps "q" and "z" to 0 instead of 7 and 9
#' 
#' letterToNumber("qz") # returns "79"
#' letterToNumber("qz", qz = 0) # returns ("00")
#' 
#' @return A character vector of numbers and dashes based on value

letterToNumber <- function(value, qz = 1) {
    value <- as.character(value)
    value <- gsub("[^A-Za-z0-9]", "-", value)
    value <- toupper(value)
    valueSplit <- strsplit(value, "")[[1]]
    numString <- as.character()
    mapphone <- function(char) { 
        if (qz == 0 && (char == LETTERS[17] | char == LETTERS[26])) { "0" } 
        else { ifelse(is.element(char, LETTERS[1:3]), "2",
               ifelse(is.element(char, LETTERS[4:6]), "3",
               ifelse(is.element(char, LETTERS[7:9]), "4",
               ifelse(is.element(char, LETTERS[10:12]), "5",
               ifelse(is.element(char, LETTERS[13:15]), "6",
               ifelse(is.element(char, LETTERS[16:19]), "7",
               ifelse(is.element(char, LETTERS[20:22]), "8",
               ifelse(is.element(char, "-") | 
                      suppressWarnings(!is.na(as.numeric(char))), char, "9")
               ))))))) 
        } 
        } 
    numString <- lapply(valueSplit, mapphone)
    return(paste0(numString, collapse = ""))
}