xNum <- function(x) {
 
  x <- as.integer(round(x,0))

       if (x == 0) num.c <- "none"
  else if (x == 1) num.c <- "one"
  else if (x == 2) num.c <- "two"
  else if (x == 3) num.c <- "three"
  else if (x == 4) num.c <- "four"
  else if (x == 5) num.c <- "five"
  else if (x == 6) num.c <- "six"
  else if (x == 7) num.c <- "seven"
  else if (x == 8) num.c <- "eight"
  else if (x == 9) num.c <- "nine"
  else if (x == 10) num.c <- "ten"
  else if (x == 11) num.c <- "eleven"
  else if (x == 12) num.c <- "twelve"
  else num.c <- .fmti(x)

  return(num.c)
}

