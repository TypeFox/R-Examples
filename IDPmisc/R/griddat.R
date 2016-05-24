## griddat.R

griddat <- function(rho,
                    cyclVar,
                    circle,
                    vp,
                    grid,
                    title)
  ## Author: Rene Locher
  ## Version: 2009-04-08
  ## helper function for plot.rose
  {
    cyclVar.lab.n <- length(grid$cyclVar$lab)

    ## Identifying position of labels on rose
    East <- 0
    South <- 0
    West <- 0
    east <- 0
    west <- 0

    if (cyclVar.lab.n>=3) {
      if (cyclVar.lab.n %%2 ==0) {
        South <- cyclVar.lab.n/2+1
        east <- 2:(cyclVar.lab.n %/% 2)
      } else east <- 2:(cyclVar.lab.n %/% 2 + 1)
      East <- which.min(abs(seq(0,by=360/cyclVar.lab.n,
                                length.out=cyclVar.lab.n)-90))
      west <- (cyclVar.lab.n %/% 2 + 2):cyclVar.lab.n
      West <- which.min(abs(seq(0,by=360/cyclVar.lab.n,
                                length.out=cyclVar.lab.n)-270))
    } else if (cyclVar.lab.n == 2) {
      South <- 2
    } else stop("'cyclVar.lab' must consist of at least 2 elements!")


    pushViewport(vp, recording=FALSE)

    circ.lab.r <- unit(grid$circ$r,"native") +
      stringWidth(grid$circ$value)* grid$circ$cex*0.5 +
      unit(grid$circ$cex,"char") * grid$circ$between
    circ.lab.r <- convertWidth(circ.lab.r,"native",valueOnly=TRUE)

    cyclVar.lab.dr <-
      if (grid$cyclVar$centered) {
        if (cyclVar.lab.n>2) {
          max(convertWidth(stringWidth(grid$cyclVar$lab[c(east,west)]),
              "native", valueOnly = TRUE) * grid$cyclVar$cex/2/
                 abs(sin(seq(0, by = 2*pi/cyclVar.lab.n,
                             length.out = cyclVar.lab.n)
                         )[c(east, west)]))
        } else
        max(convertHeight(stringHeight(grid$cyclVar$lab), "native",
                          valueOnly = TRUE)) * grid$cyclVar$cex/2
      } else 0

    each.cyclVar.lab <- grid$ray$n %/% length(grid$cyclVar$lab)
    cyclVar.lab.r <- grid$circ$r[length(grid$circ$r)] + cyclVar.lab.dr +
      convertWidth(stringWidth("X"),"native", valueOnly = TRUE)*
                   grid$cyclVar$cex*grid$cyclVar$between*0.5


    ##----------------------------------------
    ## calculating x & y coordinates

    ## for labels on circles
    circ.lab.x <- circ.lab.r*sin(grid$circ$dir)
    circ.lab.y <- circ.lab.r*cos(grid$circ$dir)

    ## for labels for cyclVar
    cyclVar.lab.x <- cyclVar.lab.r * sin(2*pi/grid$ray$n*
                       seq(0,  by = each.cyclVar.lab,
                           to = grid$ray$n-1))
    cyclVar.lab.y <- cyclVar.lab.r * cos(2*pi/grid$ray$n*
                       seq(0, by = each.cyclVar.lab,
                           to = grid$ray$n-1))

    ## adjusting coordinates when labels are not centered
    ## and calculating dimension of rose

    ## label North

    if (!grid$cyclVar$centered) cyclVar.lab.y[1] <- cyclVar.lab.y[1] +
      0.5*grid$cyclVar$cex*
        convertHeight(stringHeight(grid$cyclVar$lab[1]),
                      "native", valueOnly = TRUE)

    ymax <- convertY(unit(cyclVar.lab.y[1],"native"),"mm") +
            0.5*grid$cyclVar$cex*
              convertHeight(stringHeight(grid$cyclVar$lab[1]),"mm")

    y.title <- unit(cyclVar.lab.y[1],"native") + convertHeight(
      0.5*grid$cyclVar$cex*stringHeight(grid$cyclVar$lab[1]) +
        unit(title$between, "char")*title$cex,"mm")

    if (!is.null(title$text))
      ymax <-ymax + title$cex*stringHeight(title$text) +
        title$cex*unit(title$between, "char")

    top <- ymax -  convertY(unit(1,"npc"),"mm")

    ## label east
    if (!grid$cyclVar$centered && East>0)
      cyclVar.lab.x[east] <- cyclVar.lab.x[east] +
        0.5*grid$cyclVar$cex*
          convertWidth(stringWidth(grid$cyclVar$lab[east]),
                       "native", valueOnly = TRUE)

    if (East>0){
      xmax <-
      convertX(unit(cyclVar.lab.x[east],"native"),"mm") +
      0.5*grid$cyclVar$cex*
        convertHeight(stringWidth(grid$cyclVar$lab[east]),"mm")
      xmax <- max(xmax)
    } else
    xmax <- convertX(unit(cyclVar.lab.r,"native"), "mm")

    right <- xmax-convertX(unit(1,"npc"),"mm")

    ## label South
    if (!grid$cyclVar$centered && South>0)
      cyclVar.lab.y[South] <- cyclVar.lab.y[South] -
        0.5*grid$cyclVar$cex*convertHeight(stringHeight(
            grid$cyclVar$lab[South]), "native", valueOnly = TRUE)

    if (South>0) {
      ymin <-
        convertY(unit(cyclVar.lab.y[South],"native"),"mm") -
          0.5*grid$cyclVar$cex*
          convertHeight(stringHeight(grid$cyclVar$lab[South]),"mm")
    } else
    ymin <- convertY(unit(-cyclVar.lab.r,"native"), "mm")

    bottom <- convertY(unit(0,"npc"),"mm") - ymin

    ## label west
    if (!grid$cyclVar$centered &&West>0)
      cyclVar.lab.x[west] <- cyclVar.lab.x[west] - 0.5*grid$cyclVar$cex*
        convertWidth(stringWidth(grid$cyclVar$lab[west]),
                                 "native", valueOnly = TRUE)
    if (West>0) {
      xmin <-
        convertX(unit(cyclVar.lab.x[west],"native"),"mm") -
        0.5*grid$cyclVar$cex*
          convertWidth(stringWidth(grid$cyclVar$lab[west]), "mm")
      xmin <- min(xmin)
    } else
    xmin <- convertX(unit(-cyclVar.lab.r,"native"), "mm")

    left <- convertX(unit(0,"npc"),"mm") - xmin

    labSpace = unit.c(bottom,left,top,right)

    popViewport()

    return(list(labSpace = labSpace,
                circ = list(
                  lab = list(
                    x = circ.lab.x,
                    y = circ.lab.y)),
                cyclVar = list(
                  lab = list(
                    x = cyclVar.lab.x,
                    y = cyclVar.lab.y)),
                title = list(
                  y = y.title)))
  } ## griddat


