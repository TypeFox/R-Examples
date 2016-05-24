library(functional)

double <- Curry(`*`, e1=2)
stopifnot(double(4) == 8)

is.even <- function(a) a%%2 == 0
is.odd <- Negate(is.even)
stopifnot(Reduce(`&&`, Map(is.odd, c(1, 3, 5))))

car <- function(list) list[[1]]
cdr <- function(list) list[2:length(list)]
cadr <- Compose(cdr, car)
stopifnot(cadr(c(1,2,3)) == 2)

list.copy <- function(list)
  Reduce(Identity, list)

list <- c(1, 2, 3)
stopifnot(list.copy(list) == list)
