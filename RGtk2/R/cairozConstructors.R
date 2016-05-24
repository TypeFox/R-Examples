cairo <- cairoCreate

cairoFontFace <-
function(family, slant, weight)
{
  if (!missing(family)) {
    cairoToyFontFaceCreate(family, slant, weight)
  }
  else {
    cairoUserFontFaceCreate()
  }
}

cairoPattern <-
function(red, green, blue, alpha, surface, x0, y0, x1, y1, cx0, cy0, radius0, cx1, cy1, radius1)
{
  if (!missing(red)) {
    if (!missing(alpha)) {
      cairoPatternCreateRgba(red, green, blue, alpha)
    }
    else {
      cairoPatternCreateRgb(red, green, blue)
    }
  }
  else {
    if (!missing(surface)) {
      cairoPatternCreateForSurface(surface)
    }
    else {
      if (!missing(x0)) {
        cairoPatternCreateLinear(x0, y0, x1, y1)
      }
      else {
        cairoPatternCreateRadial(cx0, cy0, radius0, cx1, cy1, radius1)
      }
    }
  }
}

cairoSurface <-
function(width, height, format, other, content, data, stride, filename, con)
{
  if (!missing(other)) {
    cairoSurfaceCreateSimilar(other, content, width, height)
  }
  else {
    if (!missing(width)) {
      if (!missing(data)) {
        cairoImageSurfaceCreateForData(data, format, width, height, stride)
      }
      else {
        cairoImageSurfaceCreate(format, width, height)
      }
    }
    else {
      if (!missing(filename)) {
        cairoImageSurfaceCreateFromPng(filename)
      }
      else {
        cairoImageSurfaceCreateFromPngStream(con)
      }
    }
  }
}

cairoScaledFont <- cairoScaledFontCreate

cairoFontOptions <- cairoFontOptionsCreate

