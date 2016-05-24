IQrange <- function(x, na.rm = FALSE, type = 7){
  diff(quantile(as.numeric(x), c(0.25, 0.75), na.rm = na.rm, 
                names = FALSE, type = type))
}
