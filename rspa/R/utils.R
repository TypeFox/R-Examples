
all_finite <- function(x){
   stopifnot(is.numeric(x))
   storage.mode(x) <- "double"
   .Call("all_finite_double",x)
}

# code copied from print.proc.time...
getDuration <- function(x){
    y <- x
    if (!is.na(y[4L]))
        y[1L] <- y[1L] + y[4L]
    if (!is.na(y[5L]))
        y[2L] <- y[2L] + y[5L]
    y <- y[1L:3L]
    names(y) <- c("user", "system", "elapsed")
    y
}

