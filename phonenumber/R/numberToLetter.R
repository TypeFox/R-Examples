#' 
#' Convert numbers to letters as on a telephone's keypad
#' 
#' Take a character vector (i.e., a telephone number) and convert it to all
#' all possible letter combinations as on from a telephone's key pad
#' 
#' @param value An input value as a character vector with one element (a string)
#' @param decreasing Whether to sort the results in alphabetical order or not
#' @param qz Whether to assign q and z to zero (qz = 0) or not (any other value)
#' 
#' @export
#' 
#' @examples
#' # Convert a string or a vector of numeric characters
#'
#' numberToLetter("911") # returns "W11" "X11" "Y11" "Z11"
#' x <- "911"
#' numberToLetter(x) # also returns "W11" "X11" "Y11" "Z11"
#' 
#' # Convert a number directly
#' 
#' numberToLetter(911) # also returns "W11" "X11" "Y11" "Z11"
#' 
#' # Convert an alphanumeric string (letters are returned as is and
#' # non-alphanumeric characters are returned as dashes)
#' 
#' numberToLetter("dial 911!") # returns "DIAL-W11-" "DIAL-X11-" "DIAL-Y11-" "DIAL-Z11-"
#'
#' # Specifying qz = 0 maps "q" and "z" to 0 instead of 7 and 9
#' 
#' numberToLetter("000") # returns "000" 
#' numberToLetter("000", qz = 0) # returns "QQQ" "QQZ" "QZQ" "QZZ" "ZQQ" "ZQZ" "ZZQ" "ZZZ"
#' 
#' @return A character vector of letters and dashes based on value

numberToLetter <- function(value, decreasing = FALSE, qz = 1) {
    value <- as.character(value)
    value <- gsub("[^A-Za-z0-9]", "-", value)
    value <- toupper(value)
    valueSplit <- strsplit(value, "")[[1]]    
    strList <- list()
    phonemap <- function(char) { if (char == "2") { strsplit(LETTERS[1:3], "") 
        } else { if (char == "3") { strsplit(LETTERS[4:6], "") 
        } else { if (char == "4") { strsplit(LETTERS[7:9], "") 
        } else { if (char == "5") { strsplit(LETTERS[10:12], "") 
        } else { if (char == "6") { strsplit(LETTERS[13:15], "") 
        } else { if (char == "7") { strsplit(LETTERS[16:19], "") 
        } else { if (char == "8") { strsplit(LETTERS[20:22], "") 
        } else { if (char == "9") { strsplit(LETTERS[23:26], "") } else { char } 
        }}}}}}}
        }
    letterSplit <- function(char) { if (qz == 0) {if (char == "0") { 
        list(LETTERS[17], LETTERS[26]) } else { if (char == "7") { 
            c(LETTERS[16], strsplit(LETTERS[18:19], "")) } else { 
                if (char == "9") { strsplit(LETTERS[23:25], "") 
                } else { phonemap(char) }}}} else { phonemap(char) } 
        }
    strList <- lapply(valueSplit, letterSplit)
    strDF <- expand.grid(strList)
    lettString <- do.call(paste0, strDF[1:ncol(strDF)])
    lettString <- sort(lettString, decreasing)
    return(lettString)
}