## ---- echo=FALSE, warning=FALSE, error=FALSE-----------------------------
knitr::opts_chunk$set(error=TRUE)
suppressPackageStartupMessages(library(ArgumentCheck))

## ---- eval=FALSE---------------------------------------------------------
#  cylinder.volume <- function(height, radius)
#  {
#    pi * radius^2 * height
#  }

## ------------------------------------------------------------------------
cylinder.volume <- function(height, radius)
{
  if (height < 0) stop("'height' must be >= 0")
  if (radius < 0) stop("'radius' must be >= 0")
  pi * radius^2 * height  
}

## ------------------------------------------------------------------------
cylinder.volume(height = -3, 
                radius = -4)

## ------------------------------------------------------------------------
cylinder.volume(height = 3, 
                radius = -4)

## ------------------------------------------------------------------------
cylinder.volume <- function(height, radius)
{
  #* Establish a new 'ArgCheck' object
  Check <- ArgumentCheck::newArgCheck()
  
  #* Add an error if height < 0
  if (height < 0) 
    ArgumentCheck::addError(
      msg = "'height' must be >= 0",
      argcheck = Check
    )
  
  #* Add an error if radius < 0
  if (radius < 0)
    ArgumentCheck::addError(
      msg = "'radius' must be >= 0",
      argcheck = Check
    )
  
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)
  
  pi * radius^2 * height 
}

## ------------------------------------------------------------------------
cylinder.volume(height = -3,
                radius = -4)

## ------------------------------------------------------------------------
cylinder.volume(height = c(-3, 3),
                radius = -4)
cylinder.volume(height = c(3, 3),
                radius = 8)
cylinder.volume(height = c(3, -3),
                radius = c(8, 4))

## ---- eval=FALSE---------------------------------------------------------
#  if (any(height < 0))
#    ArgumentCheck::addError(
#      msg = "'height' must be >= 0",
#      argcheck = Check
#    )

## ------------------------------------------------------------------------
cylinder.volume <- function(height, radius)
{
  #* Establish a new 'ArgCheck' object
  Check <- ArgumentCheck::newArgCheck()
  
  #* Add an warning if height < 0
  if (any(height < 0)){
    ArgumentCheck::addWarning(
      msg = "'height' must be >= 0. Negative values have been set to NA",
      argcheck = Check
    )
    
    height[height < 0] <- NA
  }
  
  #* Add an error if radius < 0
  if (any(radius < 0)){
    ArgumentCheck::addWarning(
      msg = "'radius' must be >= 0. Negative values have been set to NA",
      argcheck = Check
    )
    
    radius[radius < 0] <- NA
  }
  
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)
  
  pi * radius^2 * height 
}

cylinder.volume(height = c(3, -3, 8, -1),
                radius = c(4, -4, -2, 3))

