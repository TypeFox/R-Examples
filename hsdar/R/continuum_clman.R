setMethod("initialize", signature(.Object = "Clman"),
          function(.Object, ...)
{ 
  copyFromSpeclib <- 0
  dots <- list(...)
  for (i in 1:length(dots))
    if (is.speclib(dots[[i]]))
      copyFromSpeclib <- i
  if (copyFromSpeclib == 0)
  {
    spec <- initialize(new("Speclib"), ...)
  } else {
    spec <- dots[[copyFromSpeclib]]
  }
  
  for (i in names(attributes(spec)))
  {
    if (i != "class")
    {
      if (i %in% slotNames(spec))
      {
        slot(.Object, i) <- slot(spec, i)
      } else {
        attr(.Object, i) <- attr(spec, i)
      }
    }  
  } 
  
  
  if (any(names(dots) == "hull"))
  {
    .Object@hull <- as.matrix(dots$hull)
  } else {
    .Object@hull <- matrix()
  }

  if (any(names(dots) == "cp"))
  {
    .Object@cp <- as.matrix(dots$cp)
  } else {
    .Object@cp <- matrix()
  }
  
  attr(.Object,"outdatedhull") <- if (!any(names(dots) == "outdatedhull")) c(1:nrow(.Object@cp))*0 else dots$outdatedhull
  
  if (ncol(.Object@cp) != length(.Object@wavelength))
    stop("Number of bands in continuum points and length of wavelength differ")
  if (ncol(.Object@spectra) != length(.Object@wavelength))
    stop("Number of bands in spectra and length of wavelength differ")
  if (nrow(.Object@spectra) != nrow(.Object@cp))
    stop("Number of samples in spectra and continuum points differ")
  return(.Object)
}
)

setMethod("plot", signature(x = "Clman"),
          function(
            x,
            ispec,
            subset = NULL,
            numeratepoints = TRUE, 
            hull.style = NULL,
            points.style = list(),
            ...
          )
{
  
  result <- getcp(x, ispec, subset = subset)
    
  ispec <- result$ispec
  result <- result$ptscon
  
  Reflectance  <- spectra(x)[ispec,]
  outdatedhull <- attr(x,"outdatedhull")
  if (outdatedhull[ispec]==1)
  {
    ptscon <- result$Wavelength

    hull <- Reflectance*0
    ref <- Reflectance
    storage.mode(ptscon)      <- "integer"
    storage.mode(ref)         <- "double"
    storage.mode(hull)        <- "double"
    x@hull[ispec,] <- .Fortran("makehull",
                       ncp     = as.integer(length(ptscon)), 
                       n       = as.integer(length(ref)),
                       ptscon  = ptscon, 
                       y       = ref,
                       offset  = as.integer(x$wavelength[1]-1),
                       hull    = hull,
                       PACKAGE = "hsdar"
                      )$hull
  }
  hull <- x@hull[ispec,]
  Wavelength  <- x@wavelength
    
  callNextMethod(x, FUN = ispec, ...)
  
  if (is.null(hull.style))
  {
    lty <- "dashed"
    col <- "black"
    lwd <- par()$lwd
    pch <- par()$pch
    type <- "l"
  } else {
    if (any(names(hull.style) == "lty"))
    {
      lty <- hull.style$lty
    } else {
      lty <- "dashed"
    }
    
    if (any(names(hull.style) == "col"))
    {
      col <- hull.style$col
    } else {
      col <- "black"
    }
    
    if (any(names(hull.style) == "lwd"))
    {
      lwd <- hull.style$lwd
    } else {
      lwd <- par()$lwd
    }
    
    if (any(names(hull.style) == "pch"))
    {
      pch <- hull.style$pch
    } else {
      pch <- par()$pch
    }
    
    if (any(names(hull.style) == "type"))
    {
      type <- hull.style$type
    } else {
      type <- "l"
    }
  }
  
  lines(Wavelength, hull, lty = lty, type = type, pch = pch,
        lwd = lwd, col = col)
  
  if (!is.null(points.style))
  {    
    if (any(names(points.style) == "col"))
    {
      col <- points.style$col
    } else {
      col <- "black"
    }
    
    if (any(names(points.style) == "pch"))
    {
      pch <- points.style$pch
    } else {
      pch <- par()$pch
    }
    
    if (any(names(points.style) == "lwd"))
    {
      lwd <- points.style$lwd
    } else {
      lwd <- par()$lwd
    }
    
    if (any(names(points.style) == "cex"))
    {
      cex <- points.style$cex
    } else {
      cex <- par()$cex
    }
    
    for (i in 1:nrow(result))
    {
      points(result$Wavelength[i],result$Reflectance[i], 
             cex = cex, pch = pch, lwd = lwd, col = col)
      if (numeratepoints)
        text(result$Wavelength[i],result$Reflectance[i],i, pos = 3)
    }
  }
  invisible(result)
}
)

setMethod("spectra", signature(object = "Clman"), 
          function(object, ...)
  return(callNextMethod(object, ...))
)

setReplaceMethod("spectra", signature(object = "Clman", value = "matrix"), 
                 function(object, value)
{
  return(callNextMethod(object, value))
}
)

setReplaceMethod("spectra", signature(object = "Clman", value = "data.frame"),
                 function(object, value)
{
  return(callNextMethod(object, value))
}
)

setReplaceMethod("spectra", signature(object = "Clman", value = "numeric"),
                 function(object, value)
{
  return(callNextMethod(object, value))
}
)

getcp <- function(
                  x,
                  ispec,
                  subset=NULL                 
                 )
{
  if (class(x)!="Clman")
    stop("x must be of class 'Clman'")
  
  if (length(ispec)!=1)
    stop("Multiple spectra selected")
  
  if (mode(ispec)!="numeric")
  {
    ispec <- match(ispec, idSpeclib(x), nomatch=0)
    if (ispec==0)
      stop("Unknown id of spectrum. Cannot select spectrum")
  }
  if (ispec > dim(x)[1])
    stop("ispec out of range")

  Reflectance  <- spectra(x)[ispec,]
  cp <- x@cp[ispec,]
  hull <- x@hull[ispec,]
  Wavelength  <- x@wavelength
  
  if (mode(subset)!="NULL")
  {
    if (length(subset)!=2)
      stop("subset must be a vector of length 2 giving starting and stopping wavelength")
    
    subset <- Wavelength>=subset[1] & Wavelength<=subset[2]
    Reflectance  <- Reflectance[subset]
    cp <- cp[subset]
    hull <- hull[subset]
    Wavelength  <- Wavelength[subset]
  }
  
  result <- data.frame(Wavelength=Wavelength[cp>0],
                       Reflectance=Reflectance[cp>0])
  return(list(ptscon=result, ispec=ispec))
}


deletecp <- function (
                      x,   
                      ispec,
                      cpdelete
                     )
{
  if (class(x)!="Clman")
    stop("x must be of class 'Clman'")
    
  ptscon <- getcp(x,ispec)
  ispec  <- ptscon$ispec
  outdatedhull <- attr(x,"outdatedhull")
  outdatedhull[ispec] <- 1
  attr(x,"outdatedhull") <- outdatedhull
  ptscon <- ptscon$ptscon$Wavelength
  x@cp[ispec,match(cpdelete,x@wavelength,nomatch=0)] <- 0
  return(x)
}

addcp <- function (
                   x,   
                   ispec,
                   cpadd
                  )
{
  if (class(x)!="Clman")
    stop("x must be of class 'Clman'")
    
  ptscon <- getcp(x,ispec)
  ispec  <- ptscon$ispec
  outdatedhull <- attr(x,"outdatedhull")
  outdatedhull[ispec] <- 1
  attr(x,"outdatedhull") <- outdatedhull
  ptscon <- ptscon$ptscon$Wavelength
  x@cp[ispec,match(cpadd,x@wavelength,nomatch=0)] <- cpadd
  return(x)
}