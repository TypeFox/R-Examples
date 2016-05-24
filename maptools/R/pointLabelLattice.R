# Author: Oscar Perpinan Lamigueiro (oscar.perpinan@gmail.com)
# using original code from maptools::pointLabel
# Date :  October 2012
# Version 0.10
# Licence GPL v3


drawDetails.labels <- function(x, ..., recording) {

  ##------------------------------------------------------------------##
  ## Functions SANN, GA and auxiliars (extracted from maptools::pointLabel)
  ##------------------------------------------------------------------##

  gen_offset <- function(code) c(-1, -1, -1, 0, 0, 1, 1, 1)[code] *
    (width/2) + (0+1i) * c(-1, 0, 1, -1, 1, -1, 0, 1)[code] *
      (height/2)

  rect_intersect <- function(xy1, offset1, xy2, offset2) {
    w <- pmin(Re(xy1 + offset1/2), Re(xy2 + offset2/2)) - 
      pmax(Re(xy1 - offset1/2), Re(xy2 - offset2/2))
    h <- pmin(Im(xy1 + offset1/2), Im(xy2 + offset2/2)) - 
      pmax(Im(xy1 - offset1/2), Im(xy2 - offset2/2))
    w[w <= 0] <- 0
    h[h <= 0] <- 0
    w * h
  }

  nudge <- function(offset) {
    doesIntersect <- rect_intersect(xy[rectidx1] + offset[rectidx1], 
                                    rectv[rectidx1], xy[rectidx2] + offset[rectidx2], 
                                    rectv[rectidx2]) > 0
    pyth <- abs(xy[rectidx1] + offset[rectidx1] - xy[rectidx2] - 
                offset[rectidx2])/nudgeFactor
    eps <- 1e-10
    for (i in which(doesIntersect & pyth > eps)) {
      idx1 <- rectidx1[i]
      idx2 <- rectidx2[i]
      vect <- (xy[idx1] + offset[idx1] - xy[idx2] - offset[idx2])/pyth[idx1]
      offset[idx1] <- offset[idx1] + vect
      offset[idx2] <- offset[idx2] - vect
    }
    offset
  }

  objective <- function(gene) {
    offset <- gen_offset(gene)
    if (allowSmallOverlap) offset <- nudge(offset)
    if (!is.null(rectidx1)) 
      area <- sum(rect_intersect(xy[rectidx1] + offset[rectidx1], 
                                 rectv[rectidx1], xy[rectidx2] + offset[rectidx2], 
                                 rectv[rectidx2]))
    else area <- 0
    n_outside <- sum(Re(xy + offset - rectv/2) < 0 |
                     Re(xy + offset + rectv/2) > xyAspect |
                     Im(xy + offset - rectv/2) < 0 |
                     Im(xy + offset + rectv/2) > 1/xyAspect)
    res <- 1000 * area + n_outside
    res
  }

  GA <- function() {
    n_startgenes <- 1000
    n_bestgenes <- 30
    prob <- 0.2
    mutate <- function(gene) {
      offset <- gen_offset(gene)
      doesIntersect <- rect_intersect(xy[rectidx1] + offset[rectidx1], 
                                      rectv[rectidx1], xy[rectidx2] + offset[rectidx2], 
                                      rectv[rectidx2]) > 0
      for (i in which(doesIntersect)) {
        gene[rectidx1[i]] <- sample(1:8, 1)
      }
      for (i in seq(along = gene)) if (runif(1) <= prob) 
        gene[i] <- sample(1:8, 1)
      gene
    }
    crossbreed <- function(g1, g2) ifelse(sample(c(0, 1), 
                                                 length(g1), replace = TRUE) > 0.5, g1, g2)
    genes <- matrix(sample(1:8, n_labels * n_startgenes, 
                           replace = TRUE), n_startgenes, n_labels)
    for (i in 1:10) {
      scores <- array(0, NROW(genes))
      for (j in 1:NROW(genes)) scores[j] <- objective(genes[j, 
                                                            ])
      rankings <- order(scores)
      genes <- genes[rankings, ]
      bestgenes <- genes[1:n_bestgenes, ]
      bestscore <- scores[rankings][1]
      if (bestscore == 0) {
        ## if (trace) 
        ##   cat("overlap area =", bestscore, "\n")
        break
      }
      genes <- matrix(0, n_bestgenes^2, n_labels)
      for (j in 1:n_bestgenes) {
        for (k in 1:n_bestgenes) {
          genes[n_bestgenes * (j - 1) + k, ] <-
            mutate(crossbreed(bestgenes[j,], bestgenes[k, ]))
        }}
      genes <- rbind(bestgenes, genes)
      ## if (trace) 
      ##   cat("overlap area =", bestscore, "\n")
    }
    nx <- Re(xy + gen_offset(bestgenes[1, ]))
    ny <- Im(xy + gen_offset(bestgenes[1, ]))
    list(x = nx, y = ny)
  }

  SANN <- function() {
    gene <- rep(8, n_labels)
    score <- objective(gene)
    bestgene <- gene
    bestscore <- score
    T <- 2.5
    for (i in 1:50) {
      k <- 1
      for (j in 1:50) {
        newgene <- gene
        newgene[sample(1:n_labels, 1)] <- sample(1:8, 
                                                 1)
        newscore <- objective(newgene)
        if (newscore <= score || runif(1) < exp((score - 
              newscore)/T)) {
          k <- k + 1
          score <- newscore
          gene <- newgene
        }
        if (score <= bestscore) {
          bestscore <- score
          bestgene <- gene
        }
        if (bestscore == 0 || k == 10) 
          break
      }
      if (bestscore == 0) 
        break
      ## if (trace) 
      ##   cat("overlap area =", bestscore, "\n")
      T <- 0.9 * T
    }
    ## if (trace) 
    ##   cat("overlap area =", bestscore, "\n")
    nx <- Re(xy + gen_offset(bestgene))
    ny <- Im(xy + gen_offset(bestgene))
    list(x = nx, y = ny)
  }

  ## ------------------------------------------------------------------ ##
  ## Extraction of information from the supplied grob (x)
  gl <- x
  x <- gl$x
  y <- gl$y
  labels <- gl$labels
  
  gp <- gl$gp

  allowSmallOverlap <- gl$allowSmallOverlap
  if (allowSmallOverlap) 
    nudgeFactor <- 0.02

  method <- gl$method

  ## ------------------------------------------------------------------ ##
  ## Coordinates, labels and graphical parameters
  z <- xy.coords(x, y, recycle = TRUE)
  x <- convertX(unit(z$x, 'native'), 'npc', valueOnly=TRUE)
  y <- convertY(unit(z$y, 'native'), 'npc', valueOnly=TRUE)

  if (length(labels) < length(x)) 
    labels <- rep(labels, length(x))
  n_labels <- length(x)

  vp <- current.vpTree()$parent
  windowWidth <- convertWidth(vp$width, 'cm', valueOnly=TRUE)
  windowHeight <- convertHeight(vp$height, 'cm', valueOnly=TRUE)
  xyAspect <- windowWidth/windowHeight ##windowHeight/windowWidth
  
  ## widthOriginal <- convertX(stringWidth(labels), 'npc',
  ##                           valueOnly=TRUE)
  ## heightOriginal <- convertX(stringHeight(labels), 'npc',
  ##                            valueOnly=TRUE)

  widthOriginal <-  sapply(labels, function(s){
    tg <- textGrob(s, gp=gp)
    gw <- grobWidth(tg);
    convertX(gw, 'npc', valueOnly=TRUE)
  })

  heightOriginal <-  sapply(labels, function(s){
    tg <- textGrob(s, gp=gp)
    gh <- grobHeight(tg);
    convertY(gh, 'npc', valueOnly=TRUE)
  })
  
  width <- widthOriginal + 0.015 ##xyAspect * widthOriginal  ## en pointLabel agregan 0.015
  height <- heightOriginal + 0.015 ##1/xyAspect * (heightOriginal + 1e-2)
  
  ## ------------------------------------------------------------------ ##

  xy <- x + (0+1i) * y
  rectv <- width + (0+1i) * height
  rectidx1 <- rectidx2 <- array(0, (length(x)^2 - length(x))/2)
  k <- 0
  for (i in 1:length(x)) {
    for (j in seq(len = (i - 1))) {
      k <- k + 1
      rectidx1[k] <- i
      rectidx2[k] <- j
    }}
    
  canIntersect <- rect_intersect(xy[rectidx1], 2 * rectv[rectidx1], 
                                 xy[rectidx2], 2 * rectv[rectidx2]) > 0
  rectidx1 <- rectidx1[canIntersect]
  rectidx2 <- rectidx2[canIntersect]

  ## if (trace) 
  ##   cat("possible intersects =", length(rectidx1), "\n")
  ## if (trace) 
  ##   cat("portion covered =", sum(rect_intersect(xy, rectv, 
  ##                                               xy, rectv)), "\n")
  if (method == "SANN") 
    xy <- SANN()
  else xy <- GA()

  ## ------------------------------------------------------------------ ##
  dots = list(...)
#  if (hasArg(group.number)) {
  if (!is.null(dots$group.number)) {
    group <- dots$group.number
  } else group <- 0

  idRect <- trellis.grobname('rect', type='panel', group=group)
  
  rg <- rectGrob(x=xy$x, y=xy$y,
                 width=widthOriginal, height=heightOriginal,
                 gp=gpar(fill=gp$fill, col='transparent', alpha=gp$alpha),
                 name=idRect)
  grid.draw(rg)
  
  idLabel <- trellis.grobname('label', type='panel', group=group)
  tg <- textGrob(x=xy$x, y=xy$y, label=labels,
                 gp=gp,
                 name=idLabel)
  grid.draw(tg)
}

panel.pointLabel <- function (x, y = NULL,
                              labels = seq(along = x),
                              method = c("SANN", "GA"),
                              allowSmallOverlap = FALSE,
                              col = add.text$col,
                              alpha = add.text$alpha,
                              cex = add.text$cex,
                              lineheight = add.text$lineheight,
                              font = add.text$font,
                              fontfamily = add.text$fontfamily,
                              fontface = add.text$fontface,
                              fill='transparent',
                              ...){

  add.text <- trellis.par.get("add.text")
  
  if (!missing(y) && (is.character(y) || is.expression(y))) {
    labels <- y
    y <- NULL
  }

  labels <- as.graphicsAnnot(labels)
  ## if (length(labels) < length(x)) 
  ##   labels <- rep(labels, length(x))

  method <- match.arg(method)
   
  labelGrob <- grob(x=x, y=y,
                    labels=labels,
                    gp=gpar(col = col,
                      alpha = alpha,
                      cex = cex,
                      lineheight = lineheight,
                      ##font = add.text$font,
                      fontfamily = fontfamily,
                      fontface = fontface,
                      fill = fill), ##rectangle color
                    method=method,
                    allowSmallOverlap=allowSmallOverlap,
                    cl='labels')
  grid.draw(labelGrob)
}

setGeneric('sp.pointLabel', function(object, labels, ...){standardGeneric('sp.pointLabel')})

setMethod('sp.pointLabel',
          signature=(object='SpatialPoints'),
          definition=function(object, labels, ...){
            xy = coordinates(object)
            if (missing(labels)) labels <- row.names(object)
            panel.pointLabel(xy[,1], xy[,2], labels, ...)
          })




