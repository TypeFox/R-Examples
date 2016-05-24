##' Function to translate multipliers
##'
##' This function translates number to character name or vice versa
##'
##' @param vec The vector containing name or number to be translated
##' @export
##'
##' @examples
##' ## Create numeric vector
##' myUnit = c(1000, 1e6, 1000, 1e9, 1e9, 1e12)
##'
##' ## Translate numeric to character
##' myUnit2 = translateUnit(myUnit)
##' myUnit2
##'
##' ## Now translate back
##' translateUnit(myUnit2)
##'

translateUnit = function(vec){
  type = mode(vec)
  if(type == "character"){
    if(!(all(unique(vec) %in% c("hundred", "thousand", "million",
                                "billion", "trillion", NA))))
      stop("Unrecognised name")
    transVec = double(length(vec))
    transVec[which(vec == "hundred")] = 100
    transVec[which(vec == "thousand")] = 1000
    transVec[which(vec == "million")] = 1e6
    transVec[which(vec == "billion")] = 1e9
    transVec[which(vec == "trillion")] = 1e12
    transVec[is.na(vec)] = NA
  } else if(type == "numeric"){
    if(!(all(unique(vec) %in% c(100, 1000, 1e6, 1e9, 1e12, NA))))
      stop("The unit does not have a character name available to translate")
    transVec = character(length(vec))
    transVec[which(vec == 100)] = "hundred"
    transVec[which(vec == 1000)] = "thousand"
    transVec[which(vec == 1e6)] = "million"
    transVec[which(vec == 1e9)] = "billion"
    transVec[which(vec == 1e12)] = "trillion"
    transVec[is.na(vec)] = NA
  } else {
    stop("The type of vector can not be translated")
  }
  names(transVec) = names(vec)
  transVec
}
