## test return value of "try"
isError <- function(x) {
  return(class(x) == "try-error")
}

## fetch current limits
.lim <- function(usr) {
  d <- diff(usr)/1.08
  return(usr[1:2]+d*c(0.04, -0.04))
}

## fetch current xlim
.xlim <- function() {
  return(.lim(par("usr")[1:2]))
}

## fetch current ylim
.ylim <- function() {
  return(.lim(par("usr")[3:4]))
}

## move xlim limits:
.movelim <- function(lim, width) {
  step <- diff(lim)*width
  return(lim+step)
}

## move xlim limits
##  width > 0 => to the right
##  width < 0 => to the left
.moveXlim <- function(width) {
  return(.movelim(.xlim(), width=width))
}

## move ylim limits
##  width > 0 => up
##  width < 0 => down
.moveYlim <- function(width) {
  return(.movelim(.ylim(), width=width))
}

## zoom limits
.zoomLim <- function(lim, width) {
  d <- diff(lim)/2/width
  m <- sum(lim)/2
  return(m+d*c(-1, 1))
}

## zoom only x axis
##  width < 1 => zoom out
##  width > 1 => zoom in
.zoomXlim <- function(width) {
  return(.zoomLim(.xlim(), width))
}

## zoom only y axis
##  width < 1 => zoom out
##  width > 1 => zoom in
.zoomYlim <- function(width) {
  return(.zoomLim(.ylim(), width))
}

