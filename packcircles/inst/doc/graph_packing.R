## ---- eval=FALSE---------------------------------------------------------
#  
#  library(packcircles)
#  library(ggplot2)
#  
#  ## List of tangencies. Vector elements are circle IDs.
#  ## The first element in each vector is an internal circle
#  ## and the subsequent elements are its neighbours.
#  internal <- list(
#    c(9, 4, 5, 6, 10, 11),
#    c(5, 4, 8, 6, 9),
#    c(8, 3, 2, 7, 6, 5, 4),
#    c(7, 8, 2, 1, 6)
#  )
#  
#  ## Specification of external circle radii.
#  external <- data.frame(id = c(1, 2, 3, 4, 6, 10, 11), radius=10.0)
#  
#  ## The circleGraphLayout function is used to find an arrangement
#  ## of circles corresponding to the tangencies specified by `internal`
#  ## and the external circle sizes specified by `external`. The
#  ## result is a four-column data.frame: id, x, y, radius.
#  ##
#  layout <- circleGraphLayout(internal, external)
#  
#  ## Get data for circle vertices.
#  plotdat <- circlePlotData(layout, xyr.cols = 2:4, id.col = 1)
#  
#  ## Draw circles annotated with their IDs.
#  ggplot() +
#    geom_polygon(data=plotdat, aes(x, y, group=id), fill=NA, colour="black") +
#    geom_text(data=layout, aes(x, y, label=id)) +
#    coord_equal() +
#    theme_bw()
#  

