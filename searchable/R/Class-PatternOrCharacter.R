#' @include Class-Pattern.R

NULL 

# Class: PatternOrCharacter 

  setClassUnion( 'PatternOrCharacter', c('Pattern','character'))
