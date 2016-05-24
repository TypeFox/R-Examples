
dg.prcomp <- function(data, check.plot = FALSE, UserMenus = NULL,
                      Object = new("your.Model", name = "AnModelObject")) {

  Ncases <- dim(data)[1]
  N <- dim(data)[2]

  m <- apply(data, 2, mean)
  v <- apply(data, 2, var)
  data.scaled <- t(t(data) - rep(m, dim(data)[1]))
  data.scaled <- t(t(data.scaled) / sqrt(rep(v, dim(data)[1])))

  Loadings <- prcomp(data.scaled)$rotation

  if (check.plot) {
    print(apply(data.scaled, 2, mean))
    print(apply(data.scaled, 2, var))
    par(mfrow=c(2,2))
    pc.cr <- princomp(data, cor = TRUE)
    print(summary(pc.cr))
    print(loadings(pc.cr))
    biplot(pc.cr)
    biplot(prcomp(data.scaled))
    plot(data[,1], -data[,2])
    text(data[,1], -data[,2], label=dimnames(data)[[1]])
    plot(Loadings[,1], Loadings[,2])
    text(Loadings[,1], Loadings[,2], label=dimnames(Loadings)[[1]])
  }

  Vertices <- returnVertexList(dimnames(data)[[1]], N = dim(data)[2], 
                               color = "darkred")

  ExtraVertices <- returnVertexList(c("0", dimnames(data)[[2]]),
                                    types = rep("TextVertex", dim(data)[2]+1),
                                    line = TRUE, N = dim(data)[2],
                                    color = "blue", 
                                    vertexClasses = vertexClasses())


  two.to.pairs <- function(from, to) { 
      edge.list <- vector("list", length(to))
      for (j in seq(along = to)) edge.list[[j]] <- c(from[j], to[j])
      return(edge.list) }
  extra.edge.list <- two.to.pairs(from = rep(-1, dim(data)[2]), 
                                  to = -(1:dim(data)[2])-1)
  ExtraEdges <- returnExtraEdgeList(extra.edge.list, Vertices, ExtraVertices, 
                                    color = "DarkSlateGray")

  setMethod("draw", "dg.TextVertex",
            function(object, canvas, position,
                     x = position[1], y = position[2], stratum = 0,
                     w = 2, color = "black", background = "white")
            {
              s <- w * sqrt(4 / pi) * 2
              p <- tkcreate(canvas, "oval", x - s, y - s, x + s, y + s, 
                            fill = color(object), activefill = "IndianRed")
            # l <- tkcreate(canvas, "line", 200, 200, x, y, width = w * 2, 
            #               arrow = "last", dash = "",
            #               arrowshape = paste(c(2, 5, 3) * w, collapse = " "),
            #               fill = color(object), activefill = "DarkSlateGray")
              return(list(dynamic = list(p) # , fixed = list(l)
                          )) })

  Positions(Vertices) <- 10 * matrix(unlist(c(data.scaled)), 
                                     ncol = N) %*% Loadings
  Positions(ExtraVertices) <- rbind(rep(0, N), Ncases * Loadings)

  if (FALSE) {

    Z <- DynamicGraph(Vertices = Vertices, ExtraVertices = ExtraVertices,
                      ExtraEdges = ExtraEdges, 
                      diagonal = FALSE, title = "Pre-rotated", 
                      object = Object, w = 1, UserMenus = UserMenus)
  
    Positions(Vertices) <- 10 * data.scaled
    Positions(ExtraVertices) <- rbind(rep(0, N), Ncases * diag(N))
  
    U <- DynamicGraph(Vertices = Vertices, ExtraVertices = ExtraVertices,
                      ExtraEdges = ExtraEdges, 
                      returnLink = TRUE,
                      diagonal = FALSE, title = "Master", 
                      object = Object, w = 1, UserMenus = UserMenus)
  
    W <- DynamicGraph(Vertices = Vertices, ExtraVertices = ExtraVertices,
                      ExtraEdges = ExtraEdges, 
                      transformation = t(Loadings), 
                      addModel = TRUE, frameModels = U,                    
                      diagonal = FALSE, title = "With transformation", 
                      object = Object, w = 1, UserMenus = UserMenus)

  } else {

    control <- dg.control(title = "<< Pre-rotated >>",
                          w = 1, UserMenus = UserMenus)
  
    dg <- new("dg.graph", vertexList = Vertices, 
                          extraList = ExtraVertices, 
                          extraEdgeList = ExtraEdges)
  
    Z <- dg(dg, modelObject = Object, control = control, title = "Pre-rotated")
  
    Positions(Vertices) <- 10 * data.scaled
    Positions(ExtraVertices) <- rbind(rep(0, N), Ncases * diag(N))
  
    dg2 <- new("dg.graph", vertexList = Vertices, 
                           visibleVertices = 1:Ncases,
                           extraList = ExtraVertices, 
                           extraEdgeList = ExtraEdges)
  
    control <- dg.control(title = "<< Master >>",
                          w = 1, UserMenus = UserMenus)
  
    U <- dg(dg2, modelObject = Object, control = control, title = "Master")
  
    control <- dg.control(title = "<< With transformation >>",
                          transformation = t(Loadings), 
                          w = 1, UserMenus = UserMenus)
  
    W <- addModel(dg2, frameModels = U, modelObject = Object, 
                  control = control, title = "With transformation")

  }

  invisible(U)

}
