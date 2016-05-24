TV1 <- function(x) {
 # Part of R1Magic by mehmet.suzen@physics.org
 n <-length(x)
 return( sum(abs(x[1:n-1]-x[2:n]) ) )
}
