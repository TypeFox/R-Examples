corba <-
function(xba, x) {

  cors <- mapply(function(xvar, yvar) cor(xvar, yvar, method = "pearson"), data.frame(x), data.frame(xba))
  return(mean(cors))

}
