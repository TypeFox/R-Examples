setMethod("dim", signature(x = "Speclib"), 
                 definition = function(x)
{
  dimX <- c(nrow(x@spectra), ncol(x@spectra))
  return(dimX)
}
)

nspectra <- function(x)
  return(dim(x)[1])


nbands <- function(x)
  return(dim(x)[2])


is.speclib <- function(x)
  any(c(class(x) == "Speclib",
        class(x) == "Specfeat"))