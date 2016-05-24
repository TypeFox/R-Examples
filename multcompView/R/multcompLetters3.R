"multcompLetters3" <- 
function (z, y , x, data, ...) {
  y.z  <- tapply(data[, y], data[, z], 
    function(x) do.call(mean, list(x=x)))
  oz <- order(y.z, decreasing= T )
  #This is to handle interactions
  if (length(z > 1)) {
    Lvls <- levels(interaction(data[, z], sep = ":"))[oz]
  } else {
    Lvls <- levels(data[, z])[oz]
  }
  value <- vec2mat(x)
  value <- value[Lvls, Lvls]
  multcompLetters(value, ...)
}