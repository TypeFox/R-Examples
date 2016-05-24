rawToDisplay <- function( raws )
{
  raws[ raws <= charToRaw("\037") ] <- charToRaw(".")
  raws[ raws >= charToRaw("\177") ] <- charToRaw(".")
  rawToChar(raws)
}
