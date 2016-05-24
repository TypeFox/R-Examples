library("partykit")
set.seed(1)

mysplit <- function(x) partysplit(1L, breaks = as.double(x))
x <- vector(mode = "list", length = 5)
x[[1]] <- list(id = 1L, split = mysplit(1 / 3), kids = 2:3, info = "one")
x[[2]] <- list(id = 2L, info = "two")
x[[3]] <- list(id = 3L, split = mysplit(2 / 3), kids = 4:5, info = "three")
x[[4]] <- list(id = 4L, info = "four")
x[[5]] <- list(id = 5L, info = "five")

rx <- as.partynode(x)
stopifnot(identical(as.list(rx), x))

dat <- data.frame(x = runif(100))
kidids_node(rx, dat)
