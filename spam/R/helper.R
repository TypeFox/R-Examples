# This is file ../spam/R/helper.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     



########################################################################
########################################################################
# a few nice helper functions:


bandwidth <- function(A) {
  if (!is.spam(A)) {
    warning("Matrix not 'spam' object. Coerced to one")
    A <- as.spam(A)
  }
  ret <- .Fortran("getbwd",A@dimension[1],A@entries,A@colindices,
                  A@rowpointers,low=integer(1),upp=integer(1),
                 NAOK = .Spam$NAOK, PACKAGE = "spam")
  return(c(ret$low,ret$upp))
}
                  
  

bdiag.spam <- function(...){
  nargs <- nargs()
  if (nargs == 0)     return( NULL)
  args <- list(...)
  args[which(sapply(args, is.null))] <- NULL

  if (nargs == 1)     return( args[[1]])
  if (nargs == 2) {
    # Classical case, concatenate two matrices
    A <- args[[1]]
    B <- args[[2]]
    if(!is.spam(A))
      A <- as.spam(A)
    if(!is.spam(B))
      B <- as.spam(B)
    dimA <- A@dimension
    lenA <- length(A@entries)

    A@entries <- c(A@entries,B@entries)
    A@colindices <- c(A@colindices,B@colindices+dimA[2])
    A@rowpointers <- c(A@rowpointers,B@rowpointers[-1]+lenA) 
    A@dimension <-  dimA+B@dimension
    return(A)
  } else {
    # "recursive" approach only, e.g. no checking
    tmp <- bdiag.spam( args[[1]], args[[2]])
    for ( i in 3:nargs)
      tmp <- bdiag.spam( tmp, args[[i]])
    return( tmp)
  }
}


adjacency.landkreis <- function(loc)
  # this reads the germany graph file provide by
  # loc <- "http://www.math.ntnu.no/~hrue/GMRF-book/germany.graph"
  # or
  # loc <- system.file("demodata/germany.graph", package="INLA")
  # 
  {
    n <- as.numeric( readLines(loc, n=1))

    nnodes <- nodes <- numeric( n)

    adj <- list()
    for (i in 1:n) {
      tmp <- as.numeric(scan(loc, skip=i, nlines=1, quiet=T, what=list(rep("",13)))[[1]])
      nodes[i] <- tmp[1]
      nnodes[i] <- tmp[2]
      adj[[i]] <- tmp[-c(1:2)]
    }
    
    adj <- adj[ order(nodes)]
    nnodes <- nnodes[ order(nodes)]
    A <- spam(0)
    A@colindices <- as.integer( unlist(adj)+1)
    A@rowpointers <- as.integer( c(1,cumsum(lapply(adj, length))+1))
    A@entries <- rep(1, length(unlist(adj)))
    A@dimension <- as.integer( c(n, n))
    return(A)
  }

map.landkreis <- function(data, col=NULL, zlim=range(data), add=FALSE, legendpos=c( 0.88,0.9,0.05,0.4))
# This is a stripped-down version of the function provided by the INLA package.
# Added color argument, changed 'append' to 'add'.
# Legend is tuned for a mai=rep(0,4) call 
{
  npoly <- length(spam::germany)
  ymax <- ymin <- xmax <- xmin <- 1:npoly

  if (length(data)!=npoly)
    stop('data has wrong length')
  
  if (is.null(col)) {
      if (requireNamespace("fields", quietly = TRUE)) {
          col <- fields::tim.colors(64)
      } else {
          col <- gray(seq(.05,to=0.95,length=64))
      }
  }
  ncol <- length(col)
  polycol <- col[round(((data-zlim[1])/diff(zlim)+1e-6)*(ncol-1))+1]
  
  for(i in 1:length(spam::germany)) {
    xmin[i] <- min(spam::germany[[i]][,2],na.rm=T)
    xmax[i] <- max(spam::germany[[i]][,2],na.rm=T)
    ymin[i] <- min(spam::germany[[i]][,3],na.rm=T)
    ymax[i] <- max(spam::germany[[i]][,3],na.rm=T)
  }


  if (!add)
    plot(c(min(xmin),max(xmax)),c(min(ymin),max(ymax)), type="n", axes=F, xlab="", ylab="")
  for(k in npoly:1)
    polygon(spam::germany[[k]][,2],spam::germany[[k]][,3],col=polycol[k])
  if (requireNamespace("fields", quietly = TRUE)) 
    fields::image.plot(as.matrix(data), zlim=zlim, legend.only=T, smallplot=legendpos, cex=.2, col=col)

  invisible()
}



germany.plot <- function(vect,  col=NULL, zlim=range(vect), legend=TRUE, 
             main=NULL, cex.axis=1, add=FALSE, ... )
{
  if (length(vect) != spam::germany.info$n) 
        stop("data has wrong length")

  if (!add) {
    par(mai=c(.1,.1,.1,.3))
    plot(0,0, xlim=spam::germany.info$xlim, ylim=spam::germany.info$ylim,    
         type = "n", axes = F, xlab = "", ylab = "")
  }
  if (is.null(col)) {
      ## from: require(RColorBrewer); col <- colorRampPalette(brewer.pal(9,"Blues"))(100)
    col <- c("#F7FBFF", "#F4F9FE", "#F2F8FD", "#F0F7FD", "#EEF5FC", "#ECF4FB", "#EAF3FB", "#E8F1FA", "#E6F0F9", "#E4EFF9", "#E2EEF8", "#E0ECF7", "#DEEBF7", "#DCEAF6", "#DAE8F5", "#D8E7F5", "#D6E6F4", "#D5E5F4", "#D3E3F3", "#D1E2F2", "#CFE1F2", "#CDDFF1", "#CBDEF0", "#C9DDF0", "#C7DBEF", "#C5DAEE", "#C1D9ED", "#BED7EC", "#BBD6EB", "#B8D5EA", "#B5D3E9", "#B1D2E7", "#AED1E6", "#ABCFE5", "#A8CEE4", "#A4CCE3", "#A1CBE2", "#9ECAE1", "#9AC8E0", "#96C5DF", "#92C3DE", "#8EC1DD", "#89BEDC", "#85BCDB", "#81BADA", "#7DB8DA", "#79B5D9", "#75B3D8", "#71B1D7", "#6DAFD6", "#69ACD5", "#66AAD4", "#62A8D2", "#5FA6D1", "#5CA3D0", "#58A1CE", "#559FCD", "#529DCC", "#4E9ACB", "#4B98C9", "#4896C8", "#4493C7", "#4191C5", "#3E8EC4", "#3C8CC3", "#3989C1", "#3686C0", "#3484BE", "#3181BD", "#2E7EBC", "#2C7CBA", "#2979B9", "#2776B8", "#2474B6", "#2171B5", "#1F6FB3", "#1D6CB1", "#1B69AF", "#1967AD", "#1764AB", "#1562A9", "#135FA7", "#115CA5", "#0F5AA3", "#0D57A1", "#0B559F", "#09529D", "#084F9A", "#084D96", "#084A92", "#08478E", "#08458A", "#084286", "#083F82", "#083D7E", "#083A7A", "#083776", "#083572", "#08326E", "#08306B")
  }
  ncol <- length(col)
  polycol <- (col)[round((((vect) - zlim[1])/diff(zlim) + 1e-06) * 
                            (ncol - 1)) + 1]
  
  polygon( spam::germany.poly[17965L:1L,],
          col = (polycol[spam::germany.info$polyid])[599L:1L], ...) 

  if (legend&&requireNamespace("fields", quietly = TRUE)){
    legendpos <- c(0.845, 0.89, 0.05, 0.4)
    fields::image.plot(as.matrix(vect), zlim = zlim, legend.only = TRUE, 
               smallplot = legendpos, axis.args=list(cex.axis=cex.axis,lwd=0, lwd.ticks=1.3),
               col = col)
  }
  
  if(!is.null(main))
    text( min(spam::germany.info$xlim), max(spam::germany.info$ylim), main, cex=1.5, adj=c(0,1))
  
  invisible()
}



grid_zoom <- function(inputGrob = pointsGrob(runif(200),runif(200)),
                      inputViewport = viewport(name='main'),
                      x = 'topleft', y, just,
                      ratio = c(.3,.4), zoom_xlim, zoom_ylim,
                      rect = TRUE, rect_lwd = 1, rect_fill = 'gray92',
                      draw =TRUE, zoom_fill = 'white', zoom_frame_gp = gpar(lwd = 1),
                      zoom_gp = NULL, zoom_xaxis = xaxisGrob(main = FALSE), zoom_yaxis = NULL) {

## inputGrob <- pointsGrob(runif(100), runif(100), pch='.', gp=gpar(cex=4),default.units='native',name='cc') 
## inputViewort <- viewport(name='main')
## x <- 'topleft'
## ratio <- unit(c(.3,.3), 'npc')
## zoom_xlim <- c(0.1,.5)
## zoom_ylim <- c(0.1,.5)
## rect <- TRUE
## rect_lwd = 1
## rect_fill = 'gray92'
## zoom_fill = 'white'
## zoom_frame_gp = gpar(lwd = 1)
## draw = TRUE
## zoom_gp = NULL

#  cat('may not work if other units than \'native\' and if the inputGrob has a complex structure. \n')
  
  if (!identical(length(ratio), 1)) 
    ratio <- c(ratio, ratio)
  if(class(x) == 'character')
    switch(x,
           topleft = {x = 0; y = 1; just = c(0, 1)},
           topright = {x = 1; y = 1 ; just = c(1, 1)},
           bottomright = {x = 1; y = 0; just = c(1, 0)},
           bottomleft = {x = 0; y = 0; just = c(0, 0)})

  inputViewport$name <- 'main'
  vp <- vpStack(inputViewport,
                vpList(viewport(name='zoom', x = unit(x,'npc'), y = unit(y,'npc'), width=unit(ratio[1],'npc'), height=unit(ratio[2],'npc'), just = just, xscale=zoom_xlim, yscale=zoom_ylim, default.units='native', clip = 'on'),
                       viewport(name='zoom_noClip', x = unit(x,'npc'), y = unit(y,'npc'), width=unit(ratio[1],'npc'), height=unit(ratio[2],'npc'), just = just, xscale=zoom_xlim, yscale=zoom_ylim, default.units='native', clip = 'off')))
  
  inputGrob <- editGrob(inputGrob, name='main', vp='main')
  zoom_grob <- editGrob(inputGrob, name='zoom', vp='main::zoom')

  if(!is.null(zoom_gp))
    zoom_grob <- editGrob(inputGrob, name='zoom', vp='main::zoom', gp=zoom_gp)
  
  if(!is.null(zoom_xaxis))
    zoom_xaxis <- editGrob(zoom_xaxis, vp='main::zoom_noClip', name = 'xaxis')

  if(!is.null(zoom_yaxis))
    zoom_yaxis <- editGrob(zoom_yaxis, vp='main::zoom_noClip', name = 'yaxis')

  
  rect <- rectGrob(x = zoom_xlim[1], y = zoom_ylim[1], just = c(0,0),
                    width = diff(zoom_xlim), height = diff(zoom_ylim),
                    default.units = 'native', vp = 'main', gp = gpar(lwd = rect_lwd, fill=rect_fill))

  rect_frame <- rectGrob(x = zoom_xlim[1], y = zoom_ylim[1], just = c(0,0),
                    width = diff(zoom_xlim), height = diff(zoom_ylim),
                    default.units = 'native', vp = 'main', gp = gpar(lwd = rect_lwd))
  
  
  gr <- gList(rect, inputGrob, rectGrob(vp='main::zoom', gp=gpar(fill=zoom_fill)), zoom_grob, rectGrob(vp='main::zoom_noClip', gp=zoom_frame_gp), zoom_xaxis, zoom_yaxis, rect_frame)

  out <- gTree(children=gr, childrenvp = vp)

  if (draw)
    grid.draw(out)
  invisible(out)
}


grid_trace2 <- function (chain1, chain2 = NULL,
                         xlim = NULL, ylim1 = NULL, ylim2=NULL,
                         chain1_lab = "", chain2_lab = "", main = "",
                         chain1_yaxis_at = NULL, chain2_yaxis_at = NULL,
                         log = FALSE,
                         cex_points = unit(.2, "mm"),
                         cex_main = unit(1.2, "mm"),
                         lwd_lines = unit(.1, "mm"),
                         lwd_frame = unit(.8, "mm"),
                         draw = TRUE) 
{
    if (is.null(chain2)) {
        chain2 <- chain1[, 2]
        chain1 <- chain1[, 1]
    }
    if (log) {
        chain1 <- log(chain1)
        chain2 <- log(chain2)
    }
    stopifnot(identical(length(chain1), length(chain2)))
    
    n <- length(chain1)
    if(!is.null(xlim)){
      stopifnot(length(xlim)==2)
      chain1.sub <- chain1[xlim[1]:xlim[2]]
      chain2.sub <- chain2[xlim[1]:xlim[2]]
    }else{
      chain1.sub <- chain1
      chain2.sub <- chain2
    }
    if(!is.null(ylim1))
      stopifnot(length(ylim1)==2)
    if(!is.null(ylim2))
      stopifnot(length(ylim2)==2)
       
    vp1 <- plotViewport(unit(c(2.5, 3, 2.5, 2), "cm"), name = "frame")
    vp2 <- viewport(layout = grid.layout(nrow = 1, ncol = 3, 
        respect = rbind(c(0, 0, 1)), widths = unit(c(1, 0.3, 
            1), c("null", "cm", "null"))), name = "lay1")
    vp3 <- viewport(layout.pos.col = 1, name = "traces")
    vp4 <- viewport(layout = grid.layout(nrow = 2, ncol = 1), 
        name = "lay2")
    vp5 <- viewport(layout.pos.row = 1, name = "trace1")
    vp5data <- dataViewport(xData = 1L:n, yData = chain1,
                            xscale = xlim, yscale = ylim1,
                            extension = c(0.02, 
        0.03), name = "trace1data", clip="off")
    vp5data_clip <- dataViewport(xData = 1L:n, yData = chain1,
                            xscale = xlim, yscale = ylim1,
                            extension = c(0.02, 
                              0.03), name = "trace1data_clip", clip="on")

    vp6 <- viewport(layout.pos.row = 2, name = "trace2")
    vp6data_clip <- dataViewport(xData = 1L:n, yData = chain2,
                            xscale = xlim, yscale = ylim2,
                            extension = c(0.02, 
        0.03), name = "trace2data_clip", clip="on")
    vp6data <- dataViewport(xData = 1L:n, yData = chain2,
                            xscale = xlim, yscale = ylim2,
                            extension = c(0.02, 
                              0.03), name = "trace2data", clip="off")

    vp7 <- viewport(layout.pos.col = 3, name = "scatter")
    vp7data_clip <- dataViewport(xData = chain1, yData = chain2,
                            xscale = ylim1, yscale = ylim2,
                            extension = 0.03, 
                                 name = "scatterData_clip", clip="on")
    vp7data <- dataViewport(xData = chain1, yData = chain2,
                            xscale = ylim1, yscale = ylim2,
                            extension = 0.03, 
        name = "scatterData", clip="off")

    vps <- vpStack(vp1, vp2, vpList(vpStack(vp3, vp4, vpList(vpStack(vp5, 
        vp5data, vp5data_clip), vpStack(vp6, vp6data, vp6data_clip))), vpStack(vp7, vp7data, vp7data_clip)))

   
    grs <- gList(rectGrob(vp = "frame::lay1::traces::lay2::trace1", 
        gp = gpar(lwd = lwd_frame), name = "rect_trace1"),
                 linesGrob(x = 1L:n, y = chain1, gp = gpar(lty = 1, lwd = lwd_lines), default.units = "native",  vp = "frame::lay1::traces::lay2::trace1::trace1data::trace1data_clip",  name = "lines_chain1"),
                 yaxisGrob(at = chain1_yaxis_at,
                           vp = "frame::lay1::traces::lay2::trace1::trace1data", 
                           name = "yaxis_chain1"),
                 rectGrob(vp = "frame::lay1::traces::lay2::trace2", 
                          gp = gpar(lwd = lwd_frame), name = "rect_trace2"),
                 linesGrob(x = 1L:n, 
                           y = chain2, gp = gpar(lty = 1, lwd = lwd_lines), default.units = "native", 
                           vp = "frame::lay1::traces::lay2::trace2::trace2data::trace2data_clip", 
                           name = "lines_chain2"),
                 yaxisGrob(at = chain2_yaxis_at, 
                           vp = "frame::lay1::traces::lay2::trace2::trace2data", 
                           name = "yaxis_chain2"),
                 xaxisGrob(vp = "frame::lay1::traces::lay2::trace2::trace2data", 
                           name = "xaxis_chains"),
                 pointsGrob(x = chain1.sub, y = chain2.sub, 
                            pch = 20, gp = gpar(cex = cex_points),
                            default.units = "native", 
                            vp = "frame::lay1::scatter::scatterData::scatterData_clip",
                            name = "points_scatter"), 
                 rectGrob(vp = "frame::lay1::scatter::scatterData",
                          gp = gpar(lwd = lwd_frame),  name = "rect_scatter"),
                 textGrob(chain1_lab, y = unit(-1, "line") - unit(0.2, "cm"),
                          vp = "frame::lay1::scatter", 
                          name = "text_scatter_lab1"),
                 textGrob(chain2_lab, 
                          x = unit(1, "npc") + unit(0.5, "cm"), rot = 90,
                          vp = "frame::lay1::scatter", 
                          name = "text_scatter_lab2"),
                 textGrob(main, 
                          x = unit(.5, "npc") + grobHeight(rectGrob())*.5,
                          y = unit(1, "npc") + max(stringHeight("Fg"),unit(.6, "cm")),
                          vp = "frame::lay1::traces", 
                          name = "title",
                          just=c(.5, .5),
                          gp=gpar(cex=cex_main))
                 )
    out <- gTree(childrenvp = vps, children = grs)
    if (draw) 
        grid.draw(out)
    invisible(out)
}

