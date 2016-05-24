"multcompLetters2" <- 
function (formula, x, data, ...) {
  #Convert formula to character, get rid of "~"
  fm <- as.character(formula)
  fm <- fm[-1]
  #Split char vector with ":" as this points an
  #interaction and is not included in the data
  #per se
  fm <- strsplit(fm, ":", fixed = TRUE)
  y.z  <- tapply(data[,fm[[1]]], data[,fm[[2]]], 
    function(x) do.call(mean, list(x=x)))
  oz <- order(y.z, decreasing= T )
  #This is to handle interactions
  if (length(fm[[2]] > 1)) {
    Lvls <- levels(interaction(data[,fm[[2]]], sep = ":"))[oz]
  } else {
    Lvls <- levels(data[,fm[[2]]])[oz]
  }
  value <- vec2mat(x)
  value <- value[Lvls, Lvls]
  multcompLetters(value, ...)
}

