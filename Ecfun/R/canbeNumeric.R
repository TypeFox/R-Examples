canbeNumeric <- function(x){
  ((mode(x) == 'numeric') & 
    (!('levels' %in% names(attributes(x)))) )
}
