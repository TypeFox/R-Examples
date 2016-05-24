plot.graphlist <- function(x, names = NULL, i2 = NULL, layout = NULL, 
                           plot.i1 = TRUE, max.thickness=15, circle.diameter=40, ...) {
    oldPar <- par(no.readonly = TRUE)
      ## igraph changes settings at least for abline
    graphlist <- x
    
    #######################################################################################################
    ### new circle type for vertex frame width
    mycircle <- function(coords, v=NULL, params) {
      vertex.color <- params("vertex", "color")
      if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- vertex.color[v]
      }
      vertex.size  <- 1/200 * params("vertex", "size")
      if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
      }
      vertex.frame.color <- params("vertex", "frame.color")
      if (length(vertex.frame.color) != 1 && !is.null(v)) {
        vertex.frame.color <- vertex.frame.color[v]
      }
      vertex.frame.width <- params("vertex", "frame.width")
      if (length(vertex.frame.width) != 1 && !is.null(v)) {
        vertex.frame.width <- vertex.frame.width[v]
      }
      
      mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
             vertex.size, vertex.frame.width,
             FUN=function(x, y, bg, fg, size, lwd) {
               symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                       circles=size, add=TRUE, inches=FALSE)
             })
    }
    
    add.vertex.shape("fcircle", clip=igraph.shape.noclip,
                     plot=mycircle, parameters=list(vertex.frame.color=1,
                                                    vertex.frame.width=1))
    #######################################################################################################
    
    d <- graphlist$d
    V <- graphlist$V
    i1 <- graphlist$i1[,1]
    tii <- graphlist$tii[,1]
    active <- which(tii > 0)
    tii <- tii[active]
    E <- t(combn(d,2))
    E <- E[active,]
    
    g <- graph(as.vector(t(E)) , n = d, directed = FALSE)
    if (is.null(layout)) {
        layout <- layout.fruchterman.reingold(g)
    }
    
    if (is.null(names)) {
      names <- 1:d
    }
    max.frame.width <- max.thickness
      ## thickness of the greates value
    V(g)$size <- circle.diameter
      ## vertex size (circle diameter)
  ###################################
    if(plot.i1 == FALSE){
      i1 <- rep(0,d)
    }
  ###################################
    ## scaling
    max <- max(c(tii,i1,0.000001))  # 0.000001 for the case everything is zero
    p.i1<-pmax(i1,0)
    vertex.weight.scale <- p.i1 * max.frame.width/max
    v.col <- "darkgreen"
    vertex.weight.scale <- vertex.weight.scale + 0.0001 # to avoid zero values
    edge.weight.scale <- tii * max.frame.width/max
    e.col <- "darkgreen"
    
    if (plot.i1 == FALSE){
      vertex.weight.scale = 4
      v.col <- "black"
    }
    
    #### plotting
    plot(g, layout = layout, vertex.shape = "fcircle", vertex.frame.width = vertex.weight.scale, 
        edge.width = edge.weight.scale, vertex.frame.color = v.col,
        vertex.color = "white", edge.color = e.col, vertex.label = names,
        vertex.label.color = "black")
    
    if (!is.null(i2)) {
        if (any(i2 > tii)) {
            warning("some interaction indices larger then corresponding total interaction indices, omitted")
            i2 <- pmin(i2, tii)
        }
        edge.weight.scale2 <- i2 * max.frame.width/max
        par(new = TRUE)
        plot(g, layout = layout, vertex.shape = "fcircle", vertex.frame.width = vertex.weight.scale,
            edge.width = edge.weight.scale2, vertex.frame.color = v.col,
            vertex.color = "white", edge.color = "lightgreen", vertex.label = names,
            vertex.label.color = "black")
    }
    par(oldPar)
}

