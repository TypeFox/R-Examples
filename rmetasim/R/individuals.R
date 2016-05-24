#Allan Strand 9/29/01
#
#returns a vector of populations assignments for each individual in Rland
landscape.populations <- function(Rland)
  {
    if (is.landscape(Rland))
      {
        (Rland$individuals[,1]%/%Rland$intparam$stages)+1
      }
  }


