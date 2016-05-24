ColorSet <- function(nc, ReferenceZero=NULL) {

  if (is.null(ReferenceZero)) ReferenceZero <- (nc+1)/2
  if (ReferenceZero < 1) ReferenceZero <- .5
  if (ReferenceZero > nc) ReferenceZero <- nc + .5

  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

  maxcolorset <- if (is.wholenumber(ReferenceZero)) (-nc):nc else c(-rev(1:nc), 1:nc)
  start.position <- nc - ceiling(ReferenceZero) + 2
  maxcolorset[seq(start.position, length=nc)]
}


likertColorBrewer <- function(nc, ReferenceZero=NULL,
                        BrewerPaletteName="RdBu", middle.color="gray90") {
  ## These are the diverging palettes in RColorBrewer
  ## c("RdBu", "BrBG", "PiYG", "PRGn", "PuOr",
  ## "RdGy", "RdYlBu", "RdYlGn", "Spectral"),
  ## "#F7F7F7" is the RColorBrewer default for middle.color in the RdBu scheme
  colorset <- ColorSet(nc, ReferenceZero)
  ncolors <- max(abs(colorset))*2 + (0 %in% colorset)
  oneN2 <- (1:max(abs(colorset)))
  which <- c(-rev(oneN2), 0[0 %in% colorset], oneN2) %in% colorset
  brewer.pal.likert(ncolors, BrewerPaletteName, middle.color)[which]
}


likertColor <- function(nc, ReferenceZero=NULL,
                        colorFunction=c("diverge_hcl","sequential_hcl"),
                        colorFunctionOption=c("lighter","flatter","default"),
                        colorFunctionArgs=likertColorFunctionArgs[[colorFunctionOption, colorFunction]],
                        ...) {
  colorFunction <- match.arg(colorFunction)
  colorFunctionOption <- match.arg(colorFunctionOption)
  likertColorFunctionArgs <-
    cbind(
      diverge_hcl=list(
        lighter=list(h=c(260, 0), c=80, l=c(60,90), power=.7),  ## our recommendation
        flatter=list(h=c(260, 0), c=80, l=c(30,90), power=.6),  ## our second choice
        default=list(h=c(260, 0), c=80, l=c(30,90), power=1.5)) ## colorspace default
      ,
      sequential_hcl=list(
        lighter=list( h=260, c=c(80,1), l=c(60,90), power=.7),  ## our recommendation
        flatter=list( h=260, c=c(80,1), l=c(30,90), power=.6),  ## our second choice
        default=list( h=260, c=c(80,1), l=c(30,90), power=1.5)) ## colorspace default
      )

  colorset <- ColorSet(nc, ReferenceZero)
  ncolors <- max(abs(colorset))*2 + (0 %in% colorset)
  oneN2 <- (1:max(abs(colorset)))
  which <- c(-rev(oneN2), 0[0 %in% colorset], oneN2) %in% colorset

  colorspaceFunction <- getExportedValue("colorspace", colorFunction)
   if (nc == 1 && (is.null(ReferenceZero) || ReferenceZero==1))
      do.call(colorspaceFunction, c(n=3, colorFunctionArgs))[2]
  else
    rev(do.call(colorspaceFunction, c(n=ncolors, colorFunctionArgs)))[which]
}

## Older version, selects the palette we now call "lighter".
## likertColor <- function(nc, ReferenceZero=NULL, ...) {
##   colorset <- ColorSet(nc, ReferenceZero)
##   ncolors <- max(abs(colorset))*2 + (0 %in% colorset)
##   oneN2 <- (1:max(abs(colorset)))
##   which <- c(-rev(oneN2), 0[0 %in% colorset], oneN2) %in% colorset
##   if (nc == 1 && (is.null(ReferenceZero) || ReferenceZero==1))
##       diverge_hcl(3)[2]
##   else
##     rev(diverge_hcl(ncolors))[which]
## }


brewer.pal.likert <- function(n, name,  middle.color) {

  is.odd <- function(x)  x%%2 == 1

  palette <-
    if (n <= 2) {
      bp <- RColorBrewer::brewer.pal(n=3, name=name)
      if (n==1) bp[2] else bp[-2]
    }
    else {
      if (n <= 11)
        RColorBrewer::brewer.pal(n=n, name=name)
      else {
        if (is.odd(n))
          colorRampPalette(RColorBrewer::brewer.pal(n=11, name=name))(n)
        else
          colorRampPalette(RColorBrewer::brewer.pal(n=10, name=name))(n)
      }
    }
  if (is.odd(n) && !missing(middle.color)) {
    middle <- (n %/% 2) + 1
    palette[middle] <- middle.color
  }
  palette
}

