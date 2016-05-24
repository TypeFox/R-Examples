lengths <- function(design, ...){
 UseMethod("lengths", design)
}

lengths.default <- function(design, ...){
 ## function lengths was added to base R with version 3.2.0
 if (getRversion() < "3.2.0")
 lengths.design(design, ...)
 else
 base::lengths(design)
}

lengths.design <- function (design, ...){
c('2' = length2(design, ...), 
  '3' = length3(design, ...), 
  '4' = length4(design, ...), 
  '5' = length5(design, ...))
}

lengths.matrix <- function (design, ...){
c('2' = length2(design, ...), 
  '3' = length3(design, ...), 
  '4' = length4(design, ...), 
  '5' = length5(design, ...))
}
