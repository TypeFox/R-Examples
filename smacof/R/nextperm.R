## Utility functions for unidimensional scaling
next.perm <- function(x) .C("permNext",as.double(x),as.integer(length(x)))[[1]]
are.monotone <- function(x,y) as.logical(.C("isMon",as.double(x),as.double(y),as.integer(length(x)),as.integer(1))[[4]])