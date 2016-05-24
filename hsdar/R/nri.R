nri <- function(
                x,
                b1,
                b2,
                recursive = FALSE,
                bywavelength = TRUE
               )
{
if (!is.speclib(x))
  stop("x must be of class 'Speclib'")

range.of.wavelength <- x$fwhm

reflectance <- spectra(x)
wavelength <- wavelength(x)


if (recursive)
{
 
  if (inherits(nrow(reflectance) * (sum(1:length(wavelength))-length(wavelength)), "error"))
  {
      stop("Number of Samples*(number of wavelengths^2) exceeds maximum
      vector size of 2^31-1")
  }
  
  nri_dat <- single(length = nrow(reflectance) * (sum(1:length(wavelength))-length(wavelength)))
  result <- .Fortran("recursive_nri",
                     nwl = as.integer(length(wavelength)),
                     nspec = as.integer(nrow(reflectance)),
                     reflectance = as.single(as.matrix(reflectance)),
                     nri = nri_dat,
                     nri_length = as.integer(nrow(reflectance) *
                                  (sum(1:length(wavelength))-length(wavelength)))#,
#                      PACKAGE = "hsdar"
                    )
  
  result <- distMat3D(as.numeric(result$nri), length(wavelength), nrow(reflectance))
  result <- new("Nri", nri = result, fwhm = range.of.wavelength,
                wavelength = wavelength,
                dimnames = list(Band_1 = paste("B_", wavelength, sep = ""),
                                Band_2 = paste("B_", wavelength, sep = ""),
                                Sample = idSpeclib(x)),
                attributes = attribute(x)                
               )
  if (!is.null(attr(x, "caretParameters")))
    attr(result, "caretParameters") <- attr(x, "caretParameters")
  result@usagehistory <- c(x@usagehistory, "NRI values calculated") 
} else {
  b1 <- as.vector(unlist(b1))
  b2 <- as.vector(unlist(b2))
  
  stopifnot(length(b1) == length(b2))
  
  if (length(b1) > 1)
  {
    res <- apply(matrix(1:length(b1), ncol = 1), 1, 
                 FUN = function(i, x, b1, b2, bywavelength)
                 {
                   index <- nri(x, b1 = b1[i], b2 = b2[i], bywavelength = bywavelength)
                   return(index)
                 }, x, b1, b2, bywavelength)
    colnames(res) <- paste("B", b1, "B", b2, sep = "_")
    rownames(res) <- idSpeclib(x)
    return(res)
  }                 
  if (bywavelength)
  {
    posb1 <- which(wavelength==b1)
    posb2 <- which(wavelength==b2)
  } else {
    posb1 <- b1
    posb2 <- b2
  }
  result <- (reflectance[,posb1]-reflectance[,posb2])/(reflectance[,posb1]+reflectance[,posb2])
  if (class(result)=="data.frame")
    names(result)<-"NRI"
}
return(result)  
}

