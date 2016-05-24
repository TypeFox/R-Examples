#  Allan Strand 9/27/01
#
#

#
#default seed is same as calling environment.  Zero or any positive integer
#is assigned to R's random number seed.  The type of RNG is inherited from the
#calling environment
#
test.landscape.function <- function(m1,m2)
  {
    .Call("test",m1,m2,PACKAGE = "rmetasim")
  }


