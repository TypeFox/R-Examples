##  Copyright 2005, Ross Ihaka. All Rights Reserved.
##  
##  Redistribution and use in source and binary forms, with or without
##  modification, are permitted provided that the following conditions
##  are met:
##  
##     1. Redistributions of source code must retain the above copyright notice,
##        this list of conditions and the following disclaimer.
##  
##     2. Redistributions in binary form must reproduce the above copyright
##        notice, this list of conditions and the following disclaimer in the
##        documentation and/or other materials provided with the distribution.
##  
##     3. The name of the Ross Ihaka may not be used to endorse or promote
##        products derived from this software without specific prior written
##        permission.
##  
##  THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS''
##  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
##  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
##  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ROSS IHAKA BE LIABLE FOR
##  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
##  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
##  OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
##  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
##  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
##  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
##  POSSIBILITY OF SUCH DAMAGE.

##  ----------------------------------------------------------------------------

##  An S4 Color Space Class
##
##  The following color spaces are available:
##
##     RGB       linearized sRGB space
##     sRGB      Gamma-corrected sRGB space
##     XYZ       CIE-XYZ space
##     LAB       CIE-L*a*b* space
##     polarLAB  CIE-L*a*b* space in polar coordinates
##     HSV       Hue-Saturation-Value
##     LUV       CIE-L*u*v*
##     polarLUV  CIE-L*u*v* in polar coordinates
##     HLS       Hue-Lightness-Saturation
##
##  The ``canonical'' space here is really CIE-XYZ, but in this
##  implementation all spaces are treated equally, because
##  they are all useful.
##

##  ----------------------------------------------------------------------------

## The Abstract Color Class

setClass("color", representation=list(coords="matrix"))

## Subclasses corresponding to various color spaces

setClass("RGB", contains="color")
setClass("sRGB", contains="color")
setClass("XYZ", contains="color")
setClass("LAB", contains="color")
setClass("polarLAB", contains="color")
setClass("HSV", contains="color")
setClass("HLS", contains="color")
setClass("LUV", contains="color")
setClass("polarLUV", contains="color")

setGeneric("coords", function(color) standardGeneric("coords"))

setMethod("coords", "color", function(color) {color@coords})
setMethod("show", "color", function(object) show(coords(object)))

setMethod("[", "color",
          function(x, i, j, drop=FALSE) {
            do.call(class(x), list(coords(x)[i,,drop=FALSE]))
          })

setMethod("plot", signature("color"),
          function(x, y, pch=20, cex=3)
          pairs(coords(x), col=hex(x,fix=TRUE), pch=pch, cex=cex))

CheckMatrix =
  function(x)
  if (!is.numeric(x) || length(x) < 1
      || length(dim(x)) != 2 || dim(x)[2] != 3)
  stop("invalid color matrix")

CheckBounds =
  function(x, lower, upper)
  if (any(x < lower | x > upper, na.rm=TRUE))
  warning("out of gammut color")

RGB =
  function(R, G, B, names)
  {
    if (missing(R)) return(new("RGB"))
    if (missing(names)) names = dimnames(R)[[1]]
    coords = cbind(R, if (missing(G)) NULL else G,
                      if (missing(B)) NULL else B)
    CheckMatrix(coords)
    dimnames(coords) = list(names, c("R", "G", "B"))
    new("RGB", coords = coords)
  }

sRGB =
  function(R, G, B, names)
  {
    if (missing(R)) return(new("sRGB"))
    if (missing(names)) names = dimnames(R)[[1]]
    coords = cbind(R, if (missing(G)) NULL else G,
                      if (missing(B)) NULL else B)
    CheckMatrix(coords)
    dimnames(coords) = list(names, c("R", "G", "B"))
    new("sRGB", coords = coords)
  }

XYZ =
  function(X, Y, Z, names)
  {
    if (missing(X)) return(new("XYZ"))
    if (missing(names)) names = dimnames(X)[[1]]
    coords = cbind(X, if (missing(Y)) NULL else Y,
                      if (missing(Z)) NULL else Z)
    CheckMatrix(coords)
    ## CheckBounds(coords, 0, Inf)
    dimnames(coords) = list(names, c("X", "Y", "Z"))
    new("XYZ", coords = coords)
  }

LAB =
  function(L, A, B, names)
  {
    if (missing(L)) return(new("LAB"))
    if (missing(names)) names = dimnames(L)[[1]]
    coords = cbind(L, if (missing(A)) NULL else A,
                      if (missing(B)) NULL else B)
    CheckMatrix(coords)
    ## CheckBounds(coords[,1], 0, 100)
    ## CheckBounds(coords[,2], -500, 500)
    ## CheckBounds(coords[,3], -200, 200)
    dimnames(coords) = list(names, c("L", "A", "B"))
    new("LAB", coords = coords)
  }

polarLAB =
  function(L, C, H, names)
  {
    if (missing(L)) return(new("polarLAB"))
    if (missing(names)) names = dimnames(L)[[1]]
    coords = cbind(L, if (missing(C)) NULL else C,
                      if (missing(H)) NULL else H)
    CheckMatrix(coords)
    ## CheckBounds(coords[,1], 0, Inf)
    ## CheckBounds(coords[,2], 0, Inf)
    ## CheckBounds(coords[,3], 0, 360)
    dimnames(coords) = list(names, c("L", "C", "H"))
    new("polarLAB", coords = coords)
  }

HSV =
  function(H, S, V, names)
  {
    if (missing(H)) return(new("HSV"))
    if (missing(names)) names = dimnames(H)[[1]]
    coords = cbind(H, if (missing(S)) NULL else S,
                      if (missing(V)) NULL else V)
    CheckMatrix(coords)
    ## CheckBounds(coords[,1], 0, 360)
    ## CheckBounds(coords[,2], 0, 1)
    ## CheckBounds(coords[,3], 0, 1)
    dimnames(coords) = list(names, c("H", "S", "V"))
    new("HSV", coords = coords)
  }

HLS =
  function(H, L, S, names)
  {
    if (missing(H)) return(new("HLS"))
    if (missing(names)) names = dimnames(H)[[1]]
    coords = cbind(H, if (missing(L)) NULL else L,
                      if (missing(S)) NULL else S)
    CheckMatrix(coords)
    ## CheckBounds(coords[,1], 0, 360)
    ## CheckBounds(coords[,2], 0, 1)
    ## CheckBounds(coords[,3], 0, 1)
    dimnames(coords) = list(names, c("H", "L", "S"))
    new("HLS", coords = coords)
  }

LUV =
  function(L, U, V, names)
  {
    if (missing(L)) return(new("LUV"))
    if (missing(names)) names = dimnames(L)[[1]]
    coords = cbind(L, if (missing(U)) NULL else U,
                      if (missing(V)) NULL else V)
    ## CheckBounds(coords[,1], 0, 100)
    ## CheckBounds(coords[,2], -Inf, Inf)
    ## CheckBounds(coords[,3], -Inf, Inf)
    dimnames(coords) = list(names, c("L", "U", "V"))
    new("LUV", coords = coords)
  }

polarLUV =
  function(L, C, H, names)
  {
    if (missing(L)) return(new("polarLUV"))
    if (missing(names)) names = dimnames(L)[[1]]
    coords = cbind(L, if (missing(C)) NULL else C,
                      if (missing(H)) NULL else H)
    dimnames(coords) = list(names, c("L", "C", "H"))
    new("polarLUV", coords = coords)
  }

setAs("color", "RGB", function(from)
      RGB(.Call("as_RGB", from@coords, class(from), .WhitePoint, PACKAGE = "colorspace"),
          names = dimnames(from@coords)[[1]]))

setAs("color", "sRGB", function(from)
      sRGB(.Call("as_sRGB", from@coords, class(from), .WhitePoint, PACKAGE = "colorspace"),
           names = dimnames(from@coords)[[1]]))

setAs("color", "XYZ", function(from)
      XYZ(.Call("as_XYZ", from@coords, class(from), .WhitePoint, PACKAGE = "colorspace"),
          names = dimnames(from@coords)[[1]]))

setAs("color", "LAB", function(from)
      LAB(.Call("as_LAB", from@coords, class(from), .WhitePoint, PACKAGE = "colorspace"),
          names = dimnames(from@coords)[[1]]))

setAs("color", "polarLAB", function(from)
      polarLAB(.Call("as_polarLAB", from@coords, class(from), .WhitePoint, PACKAGE = "colorspace"),
               names = dimnames(from@coords)[[1]]))

setAs("color", "HSV", function(from)
      HSV(.Call("as_HSV", from@coords, class(from), .WhitePoint, PACKAGE = "colorspace"),
          names = dimnames(from@coords)[[1]]))

setAs("color", "HLS", function(from)
      HLS(.Call("as_HLS", from@coords, class(from), .WhitePoint, PACKAGE = "colorspace"),
          names = dimnames(from@coords)[[1]]))

setAs("color", "LUV", function(from)
      LUV(.Call("as_LUV", from@coords, class(from), .WhitePoint, PACKAGE = "colorspace"),
          names = dimnames(from@coords)[[1]]))

setAs("color", "polarLUV", function(from)
      polarLUV(.Call("as_polarLUV", from@coords, class(from), .WhitePoint, PACKAGE = "colorspace"),
               names = dimnames(from@coords)[[1]]))

hex =
  function(from, gamma = NULL, fixup=FALSE)
  {
      if (!is.null(gamma))
          warning("'gamma' is deprecated and has no effect")
      coords <- as(from, "sRGB")@coords
      rval <- .Call("sRGB_to_RColor", coords, fixup, PACKAGE = "colorspace")
      if(!is.null(dnam <- attr(coords, "dimnames"))) names(rval) <- dnam[[1L]]
      return(rval)
  }

hex2RGB =
  function(x, gamma = FALSE) {
      rval <- if (gamma)
          RGB(.Call("hex_to_RGB", x, gamma, PACKAGE = "colorspace"))
      else
          sRGB(.Call("hex_to_RGB", x, gamma, PACKAGE = "colorspace"))
      if(!is.null(nam <- names(x))) attr(rval@coords, "dimnames")[[1L]] <- nam
      return(rval)
  }
  

readRGB =
  function(file, class="RGB")
  {
      x = scan(file, what=list(R = 0, G = 0, B = 0, name = ""))
      as(RGB(R = x$R/255, G = x$G/255, B = x$B/255), Class=class)
  }

readhex =
  function(file="", class="RGB")
  as(hex2RGB(scan(file, what = "")), Class=class)

writehex =
  function(x, file="")
  {
      cat(hex(x), sep="\n", file=file)
      file
  }

.WhitePoint = NULL

mixcolor = 
  function(alpha, color1, color2, where = class(color1))
  {
    alpha = as.numeric(alpha)
    c1 = coords(as(color1, where))
    c2 = coords(as(color2, where))
    na = length(alpha)
    n1 = nrow(c1)
    n2 = nrow(c2)
    n = max(na, n1, n2)
    if (na < n) alpha = rep(alpha, length = n)
    if (n1 < n) c1 = c1[rep(1:n1, length = n),]
    if (n2 < n) c2 = c2[rep(1:n2, length = n),]
    get(where)((1 - alpha) * c1 + alpha * c2)
  }
