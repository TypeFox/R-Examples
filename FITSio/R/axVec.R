`axVec` <-
function (nax = 1, axDat) 
{
### Function makes axis vector from FITS file data
###
### Takes:
  ## nax is axis number
  ## axDat is axis data from data frame produced by readFITSim.r
### Returns:
  ## Vector with axis values
### Requires/Used by:
  ## 
###
### Refs: http://fits.gsfc.nasa.gov/
###       Hanisch et al., Astr. Ap. 376, 359-380 (2001) 
###
### A. Harris, Univ. MD Astronomy, 7/19/08
###
    axVec = seq(from = (1 - axDat$crpix[nax]) * axDat$cdelt[nax] + 
        axDat$crval[nax], by = axDat$cdelt[nax], length = axDat$len[nax])
    return(axVec)
}

