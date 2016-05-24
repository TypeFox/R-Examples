leapYear <- 
function(x) {
  if (!is.numeric(x))
    stop('x must be numeric')
  x <- floor(x)
  x%%4 == 0 & (x%%100 != 0 | x%%400 == 0)
}