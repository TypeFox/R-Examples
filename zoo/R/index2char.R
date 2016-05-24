index2char <- function(x, ...) UseMethod("index2char")

index2char.default <- function(x, ...) as.character(x)

index2char.numeric <- function(x, frequency = NULL, digits = getOption("digits") - 3, ...)
{
  freq <- frequency
  if(is.null(freq)) return(as.character(round(x, digits = digits)))
  if(length(x) < 1) return(character(0))
  d <- diff(x)
  if(freq > 1 && identical(all.equal(freq, round(freq)), TRUE)) freq <- round(freq)
  if(identical(all.equal(freq*d, round(freq*d)), TRUE)) {
    if(freq == 1) return(as.character(round(x)))
      else return(paste(floor(x), "(", round((x - floor(x)) * freq) + 1, ")", sep = ""))
  } else
    return(as.character(round(x, digits = digits)))
}
