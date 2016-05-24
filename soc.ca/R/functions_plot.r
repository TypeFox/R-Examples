### Plot for all modalitites
#' Map all modalities
#' 
#' Creates a map of all active and supplementary modalities on two selected 
#' dimension.
#' @param object a soc.ca class object as created by \link{soc.mca} and 
#'   \link{soc.csa}
#' @param dim the dimensions in the order they are to be plotted. The first 
#'   number defines the horizontal axis and the second number defines the 
#'   vertical axis.
#' @param map.title the title of the map. If set to its default the standard 
#'   title is used.
#' @param labelx the label of the horizontal axis. If set to NULL a standard 
#'   label is used.
#' @param labely the label of the vertical axis. If set to NULL a standard label
#'   is used.
#' @param point.shape a numerical value defining the shape of the points. If set
#'   to its default, the default scale is used. It may be mapped to a variable 
#'   with a suitable length and order.
#' @param point.alpha defines the alpha of the points. Values range from 0 to 1.
#'   It may be mapped to a variable with a suitable length and order.
#' @param point.fill defines the fill color of the points. It may be mapped to a
#'   variable with a suitable length and order.
#' @param point.color defines the color of the points. It may be mapped to a
#'   variable with a suitable length and order. See \link{colors} for some of
#'   the valid values.
#' @param point.size a numerical value defining the size of the points. If set 
#'   to its default, the size is determined by the frequency of each modality. 
#'   It may be defined by a variable with a suitable length.
#' @param label if TRUE each point is assigned its label, defined in the soc.ca 
#'   object. See \link{assign.label} and \link{add.to.label} for ways to alter 
#'   the labels.
#' @param label.repel if TRUE overlapping labels are rearranged, see \link{geom_text_repel} or \link{geom_label_repel}.
#' @param label.alpha defines the alpha of the labels. Values range from 0 to 1.
#'   It may be mapped to a variable with a suitable length and order.
#' @param label.color defines the color of the labels. It may be mapped to a
#'   variable with a suitable length and order. See \link{colors} for some of
#'   the valid values.
#' @param label.size defines the size of the labels. It may be mapped to a
#'   variable with a suitable length and order.
#' @param label.fill defines the color of the box behind the labels. It may be mapped to a
#'   variable with a suitable length and order. This only works if label.repel is TRUE. See \link{geom_label_repel}.
#' @param legend if set to TRUE a legend is provided. Change the legend with the
#'   \link{guides}, \link{theme} and link{guide_legend} functions from the
#'   ggplot2 package.
#' @examples
#' example(soc.ca)
#' map.mod(result)
#' map.mod(result, dim = c(3, 2), point.size = 2)
#' @export
map.mod         <- function(object, dim = c(1, 2),
                            point.shape = "variable", point.alpha = 0.8,
                            point.fill = "whitesmoke", point.color = "black", point.size = "freq",
                            label = TRUE, label.repel = FALSE, label.alpha = 0.8, label.color = "black", label.size = 4, label.fill = NULL,
                            map.title = "mod", labelx = "default", labely = "default", legend = NULL){
  
  plot.type   <- "mod"

  plot.flow(object,
            dim         = dim,
            point.shape = point.shape,
            point.alpha = point.alpha,
            point.fill  = point.fill,
            point.color = point.color,
            point.size  = point.size,
            label       = label,
            label.repel = label.repel,
            label.alpha = label.alpha,
            label.color = label.color, 
            label.size  = label.size,
            label.fill  = label.fill,
            map.title   = map.title,
            labelx      = labelx,
            labely      = labely,
            legend      = legend,
            plot.type   = plot.type)
}

#' Map the most contributing modalities
#' 
#' Creates a map of the modalities contributing above average to one or more 
#' dimensions on two selected dimension.
#' @param object a soc.ca class object as created by \link{soc.mca} and 
#'   \link{soc.csa}
#' @param dim the dimensions in the order they are to be plotted. The first 
#'   number defines the horizontal axis and the second number defines the 
#'   vertical axis.
#' @param ctr.dim the dimensions of the contribution values
#' @param map.title the title of the map. If set to its default the standard 
#'   title is used.
#' @param labelx the label of the horizontal axis. If set to NULL a standard 
#'   label is used.
#' @param labely the label of the vertical axis. If set to NULL a standard label
#'   is used.
#' @param point.shape a numerical value defining the shape of the points. If set
#'   to its default, the default scale is used. It may be mapped to a variable 
#'   with a suitable length and order.
#' @param point.alpha defines the alpha of the points. Values range from 0 to 1.
#'   It may be mapped to a variable with a suitable length and order.
#' @param point.fill defines the fill color of the points. It may be mapped to a
#'   variable with a suitable length and order.
#' @param point.color defines the color of the points. It may be mapped to a
#'   variable with a suitable length and order. See \link{colors} for some of
#'   the valid values.
#' @param point.size a numerical value defining the size of the points. If set 
#'   to its default, the size is determined by the frequency of each modality. 
#'   It may be defined by a variable with a suitable length.
#' @param label if TRUE each point is assigned its label, defined in the soc.ca 
#'   object. See \link{assign.label} and \link{add.to.label} for ways to alter 
#'   the labels.
#' @param label.repel if TRUE overlapping labels are rearranged, see \link{geom_text_repel} or \link{geom_label_repel}.
#' @param label.alpha defines the alpha of the labels. Values range from 0 to 1.
#'   It may be mapped to a variable with a suitable length and order.
#' @param label.color defines the color of the labels. It may be mapped to a
#'   variable with a suitable length and order. See \link{colors} for some of
#'   the valid values.
#' @param label.size defines the size of the labels. It may be mapped to a
#'   variable with a suitable length and order.
#' @param label.fill defines the color of the box behind the labels. It may be mapped to a
#'   variable with a suitable length and order. This only works if label.repel is TRUE. See \link{geom_label_repel}.   
#' @param legend if set to TRUE a legend is provided. Change the legend with the
#'   \link{guides}, \link{theme} and link{guide_legend} functions from the
#'   ggplot2 package.
#' @examples
#' example(soc.ca)
#' map.ctr(result)
#' map.ctr(result, ctr.dim = c(1, 2))
#' @export

map.ctr         <- function(object, dim = c(1, 2), ctr.dim = 1,
                            point.shape = "variable", point.alpha = 0.8, point.fill = "whitesmoke",
                            point.color = "black", point.size = "freq",
                            label = TRUE, label.repel = TRUE, label.alpha = 0.8, label.color = "black", label.size = 4, label.fill = NULL,
                            map.title = "ctr", labelx = "default", labely = "default", legend = NULL){
  
  plot.type   <- "ctr"
  
  
  plot.flow(object,
            dim         = dim,
            ctr.dim     = ctr.dim,
            point.shape = point.shape,
            point.alpha = point.alpha,
            point.fill  = point.fill,
            point.color = point.color,
            point.size  = point.size,
            label       = label,
            label.repel = label.repel,
            label.alpha = label.alpha,
            label.color = label.color, 
            label.size  = label.size,
            label.fill  = label.fill,
            map.title   = map.title,
            labelx      = labelx,
            labely      = labely,
            legend      = legend,
            plot.type   = plot.type)
}

#' Map the active modalities
#' 
#' Creates a map of the active modalities on two selected dimensions.
#' @param object a soc.ca class object as created by \link{soc.mca} and 
#'   \link{soc.csa}
#' @param dim the dimensions in the order they are to be plotted. The first 
#'   number defines the horizontal axis and the second number defines the 
#'   vertical axis.
#' @param map.title the title of the map. If set to its default the standard 
#'   title is used.
#' @param labelx the label of the horizontal axis. If set to NULL a standard 
#'   label is used.
#' @param labely the label of the vertical axis. If set to NULL a standard label
#'   is used.
#' @param point.shape a numerical value defining the shape of the points. If set
#'   to its default, the default scale is used. It may be mapped to a variable 
#'   with a suitable length and order.
#' @param point.alpha defines the alpha of the points. Values range from 0 to 1.
#'   It may be mapped to a variable with a suitable length and order.
#' @param point.fill defines the fill color of the points. It may be mapped to a
#'   variable with a suitable length and order.
#' @param point.color defines the color of the points. It may be mapped to a
#'   variable with a suitable length and order. See \link{colors} for some of
#'   the valid values.
#' @param point.size a numerical value defining the size of the points. If set 
#'   to its default, the size is determined by the frequency of each modality. 
#'   It may be defined by a variable with a suitable length.
#' @param label if TRUE each point is assigned its label, defined in the soc.ca 
#'   object. See \link{assign.label} and \link{add.to.label} for ways to alter 
#'   the labels.
#' @param label.repel if TRUE overlapping labels are rearranged, see \link{geom_text_repel} or \link{geom_label_repel}.
#' @param label.alpha defines the alpha of the labels. Values range from 0 to 1.
#'   It may be mapped to a variable with a suitable length and order.
#' @param label.color defines the color of the labels. It may be mapped to a
#'   variable with a suitable length and order. See \link{colors} for some of
#'   the valid values.
#' @param label.size defines the size of the labels. It may be mapped to a
#'   variable with a suitable length and order.
#' @param label.fill defines the color of the box behind the labels. It may be mapped to a
#'   variable with a suitable length and order. This only works if label.repel is TRUE. See \link{geom_label_repel}.
#' @param legend if set to TRUE a legend is provided. Change the legend with the
#'   \link{guides}, \link{theme} and link{guide_legend} functions from the
#'   ggplot2 package.
#' @examples
#' example(soc.ca)
#' map.active(result)
#' map.active(result, dim = c(2, 1))
#' map.active(result, point.size = result$ctr.mod[, 1],
#'  map.title = "All active modalities with size according to contribution")
#' @export

map.active         <- function(object, dim = c(1, 2),
                               point.shape = "variable", point.alpha = 0.8, point.fill = "whitesmoke",
                               point.color = "black", point.size = "freq",
                               label = TRUE, label.repel = FALSE, label.alpha = 0.8, label.color = "black", label.size = 4, label.fill = NULL,
                               map.title = "active", labelx = "default", labely = "default", legend = NULL){
  
  plot.type   <- "active"
  
  plot.flow(object,
            dim         = dim,
            point.shape = point.shape,
            point.alpha = point.alpha,
            point.fill  = point.fill,
            point.color = point.color,
            point.size  = point.size,
            label       = label,
            label.repel = label.repel,
            label.alpha = label.alpha,
            label.color = label.color, 
            label.size  = label.size,
            label.fill  = label.fill,
            map.title   = map.title,
            labelx      = labelx,
            labely      = labely,
            legend      = legend,
            plot.type   = plot.type)
}


#' Map the supplementary modalities
#' 
#' Creates a map of the supplementary modalities on two selected dimension.
#' @param object a soc.ca class object as created by \link{soc.mca} and 
#'   \link{soc.csa}
#' @param dim the dimensions in the order they are to be plotted. The first 
#'   number defines the horizontal axis and the second number defines the 
#'   vertical axis.
#' @param map.title the title of the map. If set to its default the standard 
#'   title is used.
#' @param labelx the label of the horizontal axis. If set to NULL a standard 
#'   label is used.
#' @param labely the label of the vertical axis. If set to NULL a standard label
#'   is used.
#' @param point.shape a numerical value defining the shape of the points. If set
#'   to its default, the default scale is used. It may be mapped to a variable 
#'   with a suitable length and order.
#' @param point.alpha defines the alpha of the points. Values range from 0 to 1.
#'   It may be mapped to a variable with a suitable length and order.
#' @param point.fill defines the fill color of the points. It may be mapped to a
#'   variable with a suitable length and order.
#' @param point.color defines the color of the points. It may be mapped to a
#'   variable with a suitable length and order. See \link{colors} for some of
#'   the valid values.
#' @param point.size a numerical value defining the size of the points. If set 
#'   to its default, the size is determined by the frequency of each modality. 
#'   It may be defined by a variable with a suitable length.
#' @param label if TRUE each point is assigned its label, defined in the soc.ca 
#'   object. See \link{assign.label} and \link{add.to.label} for ways to alter 
#'   the labels.
#' @param label.repel if TRUE overlapping labels are rearranged, see \link{geom_text_repel} or \link{geom_label_repel}.
#' @param label.alpha defines the alpha of the labels. Values range from 0 to 1.
#'   It may be mapped to a variable with a suitable length and order.
#' @param label.color defines the color of the labels. It may be mapped to a
#'   variable with a suitable length and order. See \link{colors} for some of
#'   the valid values.
#' @param label.size defines the size of the labels. It may be mapped to a
#'   variable with a suitable length and order.
#' @param label.fill defines the color of the box behind the labels. It may be mapped to a
#'   variable with a suitable length and order. This only works if label.repel is TRUE. See \link{geom_label_repel}.
#' @param legend if set to TRUE a legend is provided. Change the legend with the
#'   \link{guides}, \link{theme} and link{guide_legend} functions from the
#'   ggplot2 package.
#' @examples
#' example(soc.ca)
#' map.sup(result)
#' map.sup(result, dim = c(2, 1))
#' map.sup(result, point.size = result$coord.sup[, 4],
#'  map.title = "All supplementary modalities with size according to coordinate on the 4th dimension")
#' @export

map.sup         <- function(object, dim = c(1, 2),
                            point.shape = "variable", point.alpha = 0.8, point.fill = "whitesmoke",
                            point.color = "black", point.size = "freq",
                            label = TRUE, label.repel = TRUE, label.alpha = 0.8, label.color = "black", label.size = 4, label.fill = NULL,
                            map.title = "sup", labelx = "default", labely = "default", legend = NULL){
  
  plot.type     <- "sup"
    
  plot.flow(object,
            dim         = dim,
            point.shape = point.shape,
            point.alpha = point.alpha,
            point.fill  = point.fill,
            point.color = point.color,
            point.size  = point.size,
            label       = label,
            label.repel = label.repel,
            label.alpha = label.alpha,
            label.color = label.color, 
            label.size  = label.size,
            label.fill  = label.fill,
            map.title   = map.title,
            labelx      = labelx,
            labely      = labely,
            legend      = legend,
            plot.type   = plot.type)
}

#' Map the individuals of a soc.ca analysis
#' 
#' Creates a map of the individuals on two selected dimension.
#' @param object a soc.ca class object as created by \link{soc.mca} and 
#'   \link{soc.csa}
#' @param dim the dimensions in the order they are to be plotted. The first 
#'   number defines the horizontal axis and the second number defines the 
#'   vertical axis.
#' @param map.title the title of the map. If set to its default the standard 
#'   title is used.
#' @param labelx the label of the horizontal axis. If set to NULL a standard 
#'   label is used.
#' @param labely the label of the vertical axis. If set to NULL a standard label
#'   is used.
#' @param point.shape a numerical value defining the shape of the points. It may
#'   be mapped to a variable with a suitable length and order.
#' @param point.alpha defines the alpha of the points. Values range from 0 to 1.
#'   It may be mapped to a variable with a suitable length and order.
#' @param point.fill defines the fill color of the points. It may be mapped to a
#'   variable with a suitable length and order.
#' @param point.color defines the color of the points. It may be mapped to a
#'   variable with a suitable length and order. See \link{colors} for some of
#'   the valid values.
#' @param point.size a numerical value defining the size of the points. It may
#'   be defined by a variable with a suitable length.
#' @param label if TRUE each point is assigned its label, defined in the soc.ca 
#'   object. See \link{assign.label} and \link{add.to.label} for ways to alter 
#'   the labels.
#' @param label.repel if TRUE overlapping labels are rearranged, see \link{geom_text_repel} or \link{geom_label_repel}.
#' @param label.alpha defines the alpha of the labels. Values range from 0 to 1.
#'   It may be mapped to a variable with a suitable length and order.
#' @param label.color defines the color of the labels. It may be mapped to a
#'   variable with a suitable length and order. See \link{colors} for some of
#'   the valid values.
#' @param label.size defines the size of the labels. It may be mapped to a
#'   variable with a suitable length and order.
#' @param label.fill defines the color of the box behind the labels. It may be mapped to a
#'   variable with a suitable length and order. This only works if label.repel is TRUE. See \link{geom_label_repel}.
#' @param legend if set to TRUE a legend is provided. Change the legend with the
#'   \link{guides}, \link{theme} and link{guide_legend} functions from the
#'   ggplot2 package.
#' @examples
#' example(soc.ca)
#' map.ind(result)
#' map.ind(result, map.title = "Each individual is given its shape according to a value in a factor",
#'  point.shape = active[, 1], legend = TRUE)
#' map  <- map.ind(result, map.title = "The contribution of the individuals with new scale",
#'  point.color = result$ctr.ind[, 1], point.shape = 18) 
#' map + scale_color_continuous(low = "white", high = "red")
#' quad   <- create.quadrant(result)
#' map.ind(result, map.title = "Individuals in the space given shape and color by their quadrant",
#'  point.shape = quad, point.color = quad)
#' @export

map.ind         <- function(object, dim = c(1, 2),
                            point.shape = 21, point.alpha = 0.8, point.fill = "whitesmoke",
                            point.color = "black", point.size = 3,
                            label = FALSE, label.repel = FALSE, label.alpha = 0.8, label.color = "black", label.size = 4, label.fill = NULL,
                            map.title = "ind", labelx = "default", labely = "default", legend = NULL){
  
  plot.type   <- "ind"
  
  plot.flow(object,
            dim         = dim,
            point.shape = point.shape,
            point.alpha = point.alpha,
            point.fill  = point.fill,
            point.color = point.color,
            point.size  = point.size,
            label       = label,
            label.repel = label.repel,
            label.alpha = label.alpha,
            label.color = label.color, 
            label.size  = label.size,
            label.fill  = label.fill,
            map.title   = map.title,
            labelx      = labelx,
            labely      = labely,
            legend      = legend,
            plot.type   = plot.type)
}

#' Map select modalities and individuals
#' 
#' Creates a map of selected modalities or individuals
#' @param object a soc.ca class object as created by \link{soc.mca} and 
#'   \link{soc.csa}
#' @param dim the dimensions in the order they are to be plotted. The first 
#'   number defines the horizontal axis and the second number defines the 
#'   vertical axis.
#' @param list.mod a numerical vector indicating which active modalities to 
#'   plot. It may also be a logical vector of the same length and order as the 
#'   modalities in object$names.mod.
#' @param list.sup a numerical vector indicating which supplementary modalities 
#'   to plot. It may also be a logical vector of the same length and order as 
#'   the modalities in object$names.sup.
#' @param list.ind a numerical vector indicating which individuals to plot. It 
#'   may also be a logical vector of the same length and order as the modalities
#'   in object$names.ind.
#' @param ctr.dim the dimensions of the contribution values
#' @param map.title the title of the map. If set to its default the standard 
#'   title is used.
#' @param labelx the label of the horizontal axis. If set to NULL a standard 
#'   label is used.
#' @param labely the label of the vertical axis. If set to NULL a standard label
#'   is used.
#' @param point.shape a numerical value defining the shape of the points. If set
#'   to its default, the default scale is used. It may be mapped to a variable 
#'   with a suitable length and order.
#' @param point.alpha defines the alpha of the points. Values range from 0 to 1.
#'   It may be mapped to a variable with a suitable length and order.
#' @param point.fill defines the fill color of the points. It may be mapped to a
#'   variable with a suitable length and order.
#' @param point.color defines the color of the points. It may be mapped to a
#'   variable with a suitable length and order. See \link{colors} for some of
#'   the valid values.
#' @param point.size a numerical value defining the size of the points. If set 
#'   to its default, the size is determined by the frequency of each modality. 
#'   It may be defined by a variable with a suitable length.
#' @param label if TRUE each point is assigned its label, defined in the soc.ca 
#'   object. See \link{assign.label} and \link{add.to.label} for ways to alter 
#'   the labels.
#' @param label.repel if TRUE overlapping labels are rearranged, see \link{geom_text_repel} or \link{geom_label_repel}.
#' @param label.alpha defines the alpha of the labels. Values range from 0 to 1.
#'   It may be mapped to a variable with a suitable length and order.
#' @param label.color defines the color of the labels. It may be mapped to a
#'   variable with a suitable length and order. See \link{colors} for some of
#'   the valid values.
#' @param label.size defines the size of the labels. It may be mapped to a
#'   variable with a suitable length and order.
#' @param label.fill defines the color of the box behind the labels. It may be mapped to a
#'   variable with a suitable length and order. This only works if label.repel is TRUE. See \link{geom_label_repel}.
#' @param legend if set to TRUE a legend is provided. Change the legend with the
#'   \link{guides}, \link{theme} and \link{guide_legend} functions from the
#'   ggplot2 package.
#' @param ... further arguments are currently ignored.
#' @examples
#' example(soc.ca)
#' map.select(result, map.title = "Map of the first ten modalities", list.mod = 1:10)
#' select   <- active[, 3]
#' select   <- select == levels(select)[2]
#' map.select(result, map.title = "Map of all individuals sharing a particular value",
#'  list.ind = select, point.size = 3)
#' map.select(result, map.title = "Map of both select individuals and modalities",
#'  list.ind = select, list.mod = 1:10)
#' @export

map.select         <- function(object, dim = c(1, 2), ctr.dim = 1,
                               list.mod = NULL, list.sup = NULL, list.ind = NULL,
                               point.shape = "variable", point.alpha = 0.8, point.fill = "whitesmoke",
                               point.color = "black", point.size = "freq",
                               label = TRUE, label.repel = FALSE, label.alpha = 0.8, label.color = "black", label.size = 4, label.fill = NULL,
                               map.title = "select", labelx = "default", labely = "default", legend = NULL, ...){
  
  modal.list  <- list(list.mod = list.mod, list.sup = list.sup, list.ind = list.ind)
  
  plot.type   <- "select"
  
  plot.flow(object,
            dim         = dim,
            ctr.dim     = ctr.dim,
            modal.list = modal.list,
            point.shape = point.shape,
            point.alpha = point.alpha,
            point.fill  = point.fill,
            point.color = point.color,
            point.size  = point.size,
            label       = label,
            label.repel = label.repel,
            label.alpha = label.alpha,
            label.color = label.color, 
            label.size  = label.size,
            label.fill  = label.fill,
            map.title   = map.title,
            labelx      = labelx,
            labely      = labely,
            legend      = legend,
            plot.type   = plot.type)
}


#' Add points to an existing map created by one of the soc.ca mapping functions.

#' @param object a soc.ca class object as created by \link{soc.mca} and 
#'   \link{soc.csa}
#' @param ca.map a map created using one of the soc.ca map functions
#' @param plot.type defines which type of points to add to the map. Accepted 
#'   values are: "mod", "sup", "ind", "ctr". These values correspond to the 
#'   different forms of
#' @param list.mod a numerical vector indicating which active modalities to 
#'   plot. It may also be a logical vector of the same length and order as the 
#'   modalities in object$names.mod.
#' @param list.sup a numerical vector indicating which supplementary modalities 
#'   to plot. It may also be a logical vector of the same length and order as 
#'   the modalities in object$names.sup.
#' @param list.ind a numerical vector indicating which individuals to plot. It 
#'   may also be a logical vector of the same length and order as the modalities
#'   in object$names.ind.
#' @param dim the dimensions in the order they are to be plotted. The first 
#'   number defines the horizontal axis and the second number defines the 
#'   vertical axis.
#' @param ctr.dim the dimensions of the contribution values
#' @param labelx the label of the horizontal axis. If set to NULL a standard 
#'   label is used.
#' @param labely the label of the vertical axis. If set to NULL a standard label
#'   is used.
#' @param point.shape a numerical value defining the shape of the points. If set
#'   to its default, the default scale is used. It may be mapped to a variable 
#'   with a suitable length and order.
#' @param point.alpha defines the alpha of the points. Values range from 0 to 1.
#'   It may be mapped to a variable with a suitable length and order.
#' @param point.fill defines the fill color of the points. It may be mapped to a
#'   variable with a suitable length and order.
#' @param point.color defines the color of the points. It may be mapped to a
#'   variable with a suitable length and order. See \link{colors} for some of
#'   the valid values.
#' @param point.size a numerical value defining the size of the points. If set 
#'   to its default, the size is determined by the frequency of each modality. 
#'   It may be defined by a variable with a suitable length.
#' @param label if TRUE each point is assigned its label, defined in the soc.ca 
#'   object. See \link{assign.label} and \link{add.to.label} for ways to alter 
#'   the labels.
#' @param label.repel if TRUE overlapping labels are rearranged, see \link{geom_text_repel} or \link{geom_label_repel}.
#' @param label.alpha defines the alpha of the labels. Values range from 0 to 1.
#'   It may be mapped to a variable with a suitable length and order.
#' @param label.color defines the color of the labels. It may be mapped to a
#'   variable with a suitable length and order. See \link{colors} for some of
#'   the valid values.
#' @param label.size defines the size of the labels. It may be mapped to a
#'   variable with a suitable length and order.
#' @param label.fill defines the color of the box behind the labels. It may be mapped to a
#'   variable with a suitable length and order. This only works if label.repel is TRUE. See \link{geom_label_repel}.
#' @param legend if set to TRUE a legend is provided. Change the legend with the
#'   \link{guides}, \link{theme} and link{guide_legend} functions from the
#'   ggplot2 package.
#' @examples
#' example(soc.ca)
#' original.map    <- map.sup(result)
#' map.add(result, original.map, plot.type = "ctr", ctr.dim = 2)
#' map.add(result, map.ind(result), plot.type = "select",list.ind = 1:50,
#'  point.color = "red", label = FALSE, point.size = result$ctr.ind[1:50, 1]*2000)
#' @export

map.add         <- function(object, ca.map, plot.type = NULL,
                            ctr.dim = 1, list.mod = NULL, list.sup = NULL, list.ind = NULL,
                            point.shape = "variable", point.alpha = 0.8, point.fill = "whitesmoke",
                            point.color = "black", point.size = "freq",
                            label = TRUE, label.repel = TRUE, label.alpha = 0.8, label.color = "black", label.size = 4, label.fill = NULL,
                            labelx = "default", labely = "default", legend = NULL){
  p           <- ca.map
  dim         <- ca.map$dimensions
  org.scales  <- ca.map$ca.scales
  modal.list  <- list(list.mod = list.mod, list.sup = list.sup, list.ind = list.ind)
  
  gg.data     <- data.plot(object,
                           plot.type      = plot.type,
                           dim            = dim,
                           ctr.dim        = ctr.dim,
                           modal.list     = modal.list,
                           point.size     = point.size,
                           point.variable = point.shape)
  
  if(identical(point.shape,"variable")) point.shape <- gg.data$variable
  if(identical(point.size,"freq"))      point.size  <- gg.data$freq
  
  point.l            <- list(x    = gg.data$x,
                            y     = gg.data$y,
                            shape = point.shape,
                            alpha = point.alpha,
                            fill  = point.fill,
                            color = point.color,
                            size  = point.size)
  
  p.i                <- unlist(lapply(point.l, length)) == 1
  point.attributes   <- point.l[p.i == TRUE]
  point.aes          <- point.l[p.i == FALSE]
  
  label.l            <- list(x     = gg.data$x,
                             y     = gg.data$y,
                             label = gg.data$names,
                             alpha = label.alpha,
                             color = label.color,
                             fill  = label.fill,
                             size  = label.size)
  
  t.i                <- unlist(lapply(label.l, length)) == 1
  label.attributes   <- label.l[t.i == TRUE]
  label.aes          <- label.l[t.i == FALSE]
  
  gg.input           <- list(gg.data     = gg.data,
                             label       = label,
                             repel         = label.repel,
                             point.aes   = point.aes,
                             point.attributes = point.attributes,
                             label.aes    = label.aes,
                             label.attributes = label.attributes)
  
  # Plot
  # Points
  point.attributes             <- gg.input$point.attributes
  point.attributes$mapping     <- do.call("aes", gg.input$point.aes)
  p                            <- p + do.call("geom_point", point.attributes, quote = TRUE)
  shapes                       <- c(21, 22, 23, 24, 25, 6, 0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 15, 16, 17, 18,
                                    42, 45, 61, 48, 50:120)   
  p                            <- p + scale_shape_manual(values = shapes)   
  
  # label
  if (gg.input$label == TRUE){
    label.attributes             <- gg.input$label.attributes
    label.attributes$vjust       <- 1.8
    label.attributes$family      <- "sans"
    label.attributes$lineheight  <- 0.9
    label.attributes$mapping     <- do.call("aes", gg.input$label.aes)
    
    if(is.null(gg.input$label.aes$fill) & gg.input$repel == TRUE & is.null(label.attributes$fill)){
      label.attributes$vjust       <- NULL
      p                            <- p + do.call("geom_text_repel", label.attributes, quote = TRUE)
    }
    
    if(is.null(gg.input$label.aes$fill) == FALSE | is.null(label.attributes$fill) == FALSE){
      label.attributes$vjust       <- NULL
      p                            <- p + do.call("geom_label_repel", label.attributes, quote = TRUE)
    }
    
    if(is.null(gg.input$label.aes$fill) & gg.input$repel == FALSE){
      p                            <- p + do.call("geom_text", label.attributes, quote = TRUE)
    }
    
  }
  
  # Scales and labels
  
  org.max     <- org.scales$lim.max
  org.min     <- org.scales$lim.min
  n.max       <- max(c(max(gg.data[,1:2]), org.max))
  n.min       <- (c(min(gg.data[,1:2]), org.min))
  tscales     <- gg.data
  tscales[1, 1:2] <- n.max
  tscales[2, 1:2] <- n.min
  scales      <- breaksandscales(tscales)
  
  p           <- p + scale_x_continuous(breaks = scales$scalebreaks, labels = scales$breaklabel)
  p           <- p + scale_y_continuous(breaks = scales$scalebreaks, labels = scales$breaklabel)
  p$ca.scales <- scales
  
  #### Legend
  if(is.null(legend))            p  <- p + theme(legend.position = "none")
  if(identical(legend,"bottom")) p  <- p + theme(legend.position = 'bottom',
                                                 legend.direction = "horizontal", legend.box = "horizontal")
  p
}


#' Density plot for the cloud of individuals
#' 
#' Draws a 2d density plot on top of an existing soc.ca map. The density is 
#' calculated by the \link[MASS]{kde2d} function from MASS and plotted by
#' \link[ggplot2]{geom_density2d} from \code{ggplot2} \code{map.density} uses the
#' coordinates of the individuals as a basis for the density calculation. 
#' Borders are arbitrary.
#' 
#' @param object a soc.ca class object
#' @param map a soc.ca map object created by one of the soc.ca mapping functions
#' @param group a factor determining group membership. Density is mapped for
#'   each group individually.
#' @param color a single value or vector determining the color. See the scale
#'   functions of \code{ggplot2} for ways to alter the scales.
#' @param alpha a single value or vector determining the alpha.
#' @param size a single value or vector determining the size of the lines.
#' @param linetype a single value or vector determining the linetype
#' @export
#' @examples
#' example(soc.ca)
#' map.density(result, map.ind(result, dim = 2:3, point.alpha = 0.2))
#' map.density(result, map.ind(result, legend = TRUE, point.alpha = 0.2),
#'  group = duplicated(active), color = duplicated(active),
#'  linetype = duplicated(active))
#' map.density(result, map.ctr(result))

map.density  <- function(object, map = map.ind(object), group = NULL,
                         color = "red", alpha = 0.8, size = 0.5, linetype = "solid"){
  
  dim                        <- map$dimensions 
  dens.data                  <- as.data.frame(object$coord.ind[, dim])
  colnames(dens.data)        <- c("x", "y")
  
  density.l                  <- list(color = color,
                                     alpha = alpha,
                                     size  = size,
                                     linetype = linetype,
                                     n  = 100,
                                     group = group,
                                     na.rm = TRUE)
  
  d.i                        <- unlist(lapply(density.l, length)) == 1
  density.attributes         <- density.l[d.i]
  density.aes                <- density.l[unlist(lapply(density.l, length)) > 1]
  density.aes$x              <- dens.data$x
  density.aes$y              <- dens.data$y
  
  
  density.attributes$mapping <- do.call("aes", density.aes)
  p     <- map + do.call("geom_density2d", density.attributes, quote = TRUE)
  p
}





############################## Plot delfunktioner ########################################

plot.flow   <- function(object, dim = c(1, 2), ctr.dim = NULL, modal.list = NULL,
                        point.shape = 21, point.alpha = 0.8, point.fill = "whitesmoke",
                        point.color = "black", point.size = 3,
                        label = FALSE, label.repel = FALSE, label.alpha = 0.8, label.color = "black", label.size = 4, label.fill = NULL,
                        map.title = map.title, labelx = "default", labely = "default", legend = NULL,
                        plot.type = plot.type){
  
  
  gg.proc     <- round(object$adj.inertia[,4])                        # Adjusted inertia
  gg.data     <- data.plot(object,
                           plot.type   = plot.type,
                           dim,
                           ctr.dim     = ctr.dim,
                           modal.list  = modal.list,
                           point.size  = point.size,
                           point.color = point.color)                 # Data selection
  
  axis.labels <- plot.axis(labelx = labelx, labely = labely, gg.proc = gg.proc, dim = dim) # Axis labels
  map.title   <- plot.title(map.title = map.title, ctr.dim = ctr.dim) # Plot title
  scales      <- breaksandscales(gg.data)                             # Scales og breaks
  
  if (identical(point.shape,"variable")) point.shape <- gg.data$variable
  if (identical(point.size,"freq"))      point.size  <- gg.data$freq
  
  point.l            <- list(x    = gg.data$x,
                            y     = gg.data$y,
                            shape = point.shape,
                            alpha = point.alpha,
                            fill  = point.fill,
                            color = point.color,
                            size  = point.size)
  
  p.i                <- unlist(lapply(point.l, length)) == 1
  point.attributes   <- point.l[p.i == TRUE]
  point.aes          <- point.l[p.i == FALSE]
  
  label.l            <- list(x     = gg.data$x,
                             y     = gg.data$y,
                             label = gg.data$names,
                             alpha = label.alpha,
                             color = label.color,
                             fill  = label.fill,
                             size  = label.size)
  
  t.i                <- unlist(lapply(label.l, length)) == 1
  label.attributes   <- label.l[t.i == TRUE]
  label.aes          <- label.l[t.i == FALSE]
  
  gg.input           <- list(gg.data       = gg.data,
                             axis.inertia  = gg.proc,
                             map.title     = map.title,
                             labelx        = axis.labels$x,
                             labely        = axis.labels$y,
                             scales        = scales,
                             label         = label,
                             repel         = label.repel,
                             point.aes     = point.aes,
                             point.attributes = point.attributes,
                             label.aes     = label.aes,
                             label.attributes = label.attributes)
  
  b.plot              <- basic.plot(gg.input) 
  t.plot              <- b.plot + theme_min() 
  
  if(is.null(legend)) t.plot   <- t.plot + theme(legend.position = "none")
  if(identical(legend,"bottom")) t.plot  <- t.plot + theme(legend.position = 'bottom',
                                                           legend.direction = "horizontal", legend.box = "horizontal")
  t.plot$dimensions <- dim
  return(t.plot)
}

#################################################################################
# Soc.ca.graphics environment

# soc.ca.scales                  <- function(...){
# soc.ca.graphics                <- new.env()
# assign("scale_colour_continuous", function(...) ggplot2::scale_colour_continuous(..., low = getOption("soc.ca.gradient")[1], high = getOption("soc.ca.gradient")[2]), envir = soc.ca.graphics)
# assign("scale_color_continuous",  function(...) ggplot2::scale_colour_continuous(..., low = getOption("soc.ca.gradient")[1], high = getOption("soc.ca.gradient")[2]), envir = soc.ca.graphics)
# assign("scale_fill_continuous",   function(...) ggplot2::scale_fill_continuous(...,   low = getOption("soc.ca.gradient")[1], high = getOption("soc.ca.gradient")[2]), envir = soc.ca.graphics)
# 
# assign("scale_colour_gradient",   function(...) ggplot2::scale_colour_continuous(..., low = getOption("soc.ca.gradient")[1], high = getOption("soc.ca.gradient")[2]), envir = soc.ca.graphics)
# assign("scale_color_gradient",    function(...) ggplot2::scale_colour_continuous(..., low = getOption("soc.ca.gradient")[1], high = getOption("soc.ca.gradient")[2]), envir = soc.ca.graphics)
# assign("scale_fill_gradient",     function(...) ggplot2::scale_fill_continuous(...,   low = getOption("soc.ca.gradient")[1], high = getOption("soc.ca.gradient")[2]), envir = soc.ca.graphics)
# 
# assign("scale_colour_discrete",   function(...) ggplot2::scale_colour_manual(..., values = getOption("soc.ca.colors")), envir = soc.ca.graphics)
# assign("scale_color_discrete",    function(...) ggplot2::scale_colour_manual(..., values = getOption("soc.ca.colors")), envir = soc.ca.graphics)
# assign("scale_fill_discrete",     function(...) ggplot2::scale_fill_manual(...,   values = getOption("soc.ca.colors")), envir = soc.ca.graphics)
# if(("soc.ca.graphics" %in% search()) == FALSE) attach(soc.ca.graphics)
# 
# 
# }

basic.plot <- function(gg.input){ 
  
  #####################
  # Changing default ggplot2 behaviour
 # soc.ca.scales() # When loading the soc.ca.graphics environment, certain functions are masked from ggplot2. This is intended.
  # Her skal der perhaps en on.exit
  # Defining point.size
  p       <- ggplot()   
  # The middle plot axis
  p       <- p + geom_hline(yintercept = 0, color = "grey50", size = 0.5, linetype = "solid")
  p       <- p + geom_vline(xintercept = 0, color = "grey50", size = 0.5, linetype = "solid")
  
  # Points
  point.attributes             <- gg.input$point.attributes
  point.attributes$mapping     <- do.call("aes", gg.input$point.aes)
  p                            <- p + do.call("geom_point", point.attributes, quote = TRUE)
  shapes                       <- getOption("soc.ca.shape")  
  p                            <- p + scale_shape_manual(values = shapes)   
  
  # label
  if (gg.input$label == TRUE){
    label.attributes             <- gg.input$label.attributes
    label.attributes$vjust       <- 1.8
    label.attributes$family      <- "sans"
    label.attributes$lineheight  <- 0.9
    label.attributes$mapping     <- do.call("aes", gg.input$label.aes)
    
    if(is.null(gg.input$label.aes$fill) & gg.input$repel == TRUE & is.null(label.attributes$fill)){
    label.attributes$vjust       <- NULL
    label.attributes$max.iter    <- 2000
    p                            <- p + do.call("geom_text_repel", label.attributes, quote = TRUE)
    }
    
    if(is.null(gg.input$label.aes$fill) == FALSE | is.null(label.attributes$fill) == FALSE){
      label.attributes$vjust       <- NULL
      label.attributes$max.iter    <- 2000
      p                            <- p + do.call("geom_label_repel", label.attributes, quote = TRUE)
    }
    
    if(is.null(gg.input$label.aes$fill) & gg.input$repel == FALSE){
      p                            <- p + do.call("geom_text", label.attributes, quote = TRUE)
    }
    
  }
  
  # Title and axis labels
  p     <- p + ggtitle(label = gg.input$map.title)
  p     <- p + xlab(gg.input$labelx) + ylab(gg.input$labely)
  p     <- p + labs(alpha = "Alpha", shape = "Shape",
                    color = "Color", linetype = "Linetype", size = "Size", fill = "Fill")
  # Scales and breaks
  p     <- p + scale_x_continuous(breaks = gg.input$scales$scalebreaks, labels = gg.input$scales$breaklabel)
  p 		<- p + scale_y_continuous(breaks = gg.input$scales$scalebreaks, labels = gg.input$scales$breaklabel)
  p     <- p + coord_fixed()
  p$ca.scales <- gg.input$scales
  
 # detach(soc.ca.graphics) # Detaches the environment containing the soc.ca version of the scales, see(soc.ca.scales)
  
  p
}




####################### Plot title

plot.title  <- function(map.title = NULL, ctr.dim = NULL){
    # Map title for ctr plot    
    if (identical(map.title, "ctr") == TRUE) {
        map.title     <- paste("The modalities contributing above average to dimension ",
                               paste(ctr.dim, collapse = " & "), sep = "")
    }
    # Map title for all plot
    if (identical(map.title, "mod") == TRUE) {
        map.title     <- "Map of all modalities"
    }
    # Map title for active plot
    if (identical(map.title, "active") == TRUE) {
        map.title     <- "Map of active modalities"
    }
    # Map title for sup plot
    if (identical(map.title, "sup") == TRUE) {
        map.title     <- "Map of supplementary modalities"
    }
    # Map title for id plot
    if (identical(map.title, "ind") == TRUE) {
        map.title     <- "Map of individuals"
    }
    # Map title for list plot
    if (identical(map.title, "select") == TRUE) {
        map.title     <- "Map of selected modalities"
    }
    return(map.title)
}

############# Axis labels
plot.axis       <- function(labelx = "default", labely = "default", gg.proc = NA, dim = NULL){
    if (identical(labelx, "default") == TRUE) {
labelx         <- paste(dim[1], ". Dimension: ", gg.proc[dim[1]], "%", sep = "")
}
    if (identical(labely, "default") == TRUE) {
labely         <- paste(dim[2], ". Dimension: ", gg.proc[dim[2]], "%", sep = "")
}
axis.labels <- list(x = labelx, y = labely)
return(axis.labels)    
}

####################### breaksandscales
breaksandscales <- function(gg.data){
    mround <- function(x, base){
    base*round(x/base)
  }
  
  lim.min.x   	<- min(gg.data[, 1])
  lim.min.y   	<- min(gg.data[, 2])
  lim.max.x   	<- max(gg.data[, 1])
  lim.max.y   	<- max(gg.data[, 2])
  
  scalebreaks   <- seq(-10,10, by = 0.25)
  nolabel       <- seq(-10, 10, by = 0.5)
  inter         <- intersect(scalebreaks, nolabel)
  truelabel     <- is.element(scalebreaks, nolabel)
  breaklabel    <- scalebreaks  
  breaklabel[truelabel == FALSE] <-  ""
    
  length.x      <- sqrt(lim.min.x^2) + sqrt(lim.max.x^2)
  length.y      <- sqrt(lim.min.y^2) + sqrt(lim.max.y^2)
  
  scales        <- list(scalebreaks = scalebreaks,
                        breaklabel  = breaklabel,
                        lim.min.x   = lim.min.x,
                        lim.min.y   = lim.min.y,
                        lim.max.x   = lim.max.x,
                        lim.max.y   = lim.max.y,
                        length.x    = length.x,
                        length.y    = length.y)
  return(scales)
}

#################### data.plot v.2

data.plot   <- function(object, plot.type, dim, ctr.dim = NULL,
                        modal.list = NULL, point.size = "freq", point.variable = NULL, point.color = NULL){
  # Types: mod, active, sup, ctr, ind, list
  
  # mod
  if (identical(plot.type, "mod") == TRUE){
    coord       <- rbind(object$coord.mod, object$coord.sup)
    mnames  	  <- c(object$names.mod, object$names.sup)
    freq        <- c(object$freq.mod, object$freq.sup)
    variable    <- c(object$variable, object$variable.sup)
  }
  
  # ctr
  if (identical(plot.type, "ctr") == TRUE){
    av.ctr      <- contribution(object, ctr.dim, indices = TRUE, mode = "mod")
    coord     	<- object$coord.mod[av.ctr, ]
    mnames		  <- object$names.mod[av.ctr]
    freq        <- object$freq.mod[av.ctr]
    variable    <- object$variable[av.ctr]
  }
  # active
  if (identical(plot.type, "active") == TRUE){
    coord 		  <- object$coord.mod
    mnames		  <- object$names.mod
    freq        <- object$freq.mod
    variable    <- object$variable
  }
  # sup
  if (identical(plot.type, "sup") == TRUE){
    coord 		<- object$coord.sup
    mnames		<- object$names.sup
    freq      <- object$freq.sup
    variable  <- object$variable.sup
  }
  # ind
  if (identical(plot.type, "ind") == TRUE){
    coord 		<- object$coord.ind
    mnames		<- object$names.ind
    freq      <- rep(1, object$n.ind)
     if(identical(point.variable, NULL)){
       variable  <- rep("ind", nrow(object$coord.ind))
     }else{ 
       variable  <- point.variable
    }
    if(identical(point.color, NULL) == FALSE)
    point.color  <- point.color
  }
  
  # select
  if (identical(plot.type, "select") == TRUE){
    coord.mod 		<- object$coord.mod[modal.list$list.mod, ]
    coord.sup     <- object$coord.sup[modal.list$list.sup, ]
    coord.ind     <- object$coord.ind[modal.list$list.ind, ]
    names.mod		  <- object$names.mod[modal.list$list.mod]
    names.sup  	  <- object$names.sup[modal.list$list.sup]
    names.ind     <- object$names.ind[modal.list$list.ind]
    coord         <- rbind(coord.mod, coord.sup, coord.ind)
    rownames(coord) <- NULL
    mnames        <- c(names.mod, names.sup, names.ind)
    freq.mod      <- object$freq.mod[modal.list$list.mod]
    freq.sup      <- object$freq.sup[modal.list$list.sup]
    freq.ind      <- rep(1, object$n.ind)[modal.list$list.ind]
    freq          <- c(freq.mod, freq.sup, freq.ind)
    variable.mod  <- object$variable[modal.list$list.mod]
    variable.sup  <- object$variable.sup[modal.list$list.sup]
    variable.ind  <- rep("ind", object$n.ind)[modal.list$list.ind]
    variable      <- c(variable.mod, variable.sup, variable.ind)
  }
  
  if(is.numeric(point.size))  freq         <- rep(point.size, length.out = nrow(coord))
  if(is.null(point.color))    point.color  <- rep("Nothing", length.out = nrow(coord))
  
  
  gg.data             <- data.frame(cbind(coord[, dim[1]], coord[,dim[2]]), mnames, freq, variable, point.color)
  colnames(gg.data)   <- c("x", "y", "names", "freq", "variable", "point.color")
  return(gg.data)
}

#' Add a new layer of points on top of an existing plot with output from the 
#' min_cut function
#' @param x a matrix created by the min_cut function
#' @param p is a ggplot object, preferably from one of the mapping functions in 
#'   soc.ca
#' @param label if TRUE the labels of points will be shown
#' @param ... further arguments are passed on to geom_path, geom_point and 
#'   geom_text

add.count <- function(x, p, label = TRUE, ...){
  p   <- p + geom_point(data = x, x = x$X, y = x$Y, ...) + geom_path(data = x, x = x$X, y = x$Y, ...)
  if (identical(label, TRUE)) p <- p + geom_text(data = x, x = x$X, y = x$Y, label = x$label, vjust = 0.2, ...)
}

#' Map path along an ordered variable
#' 
#' Plot a path along an ordered variable. If the variable is numerical it is cut
#' into groups by the \link{min_cut} function.
#' 
#' @param object is a soc.ca result object
#' @param x is an ordered vector, either numerical or factor
#' @param map is a plot object created with one of the mapping functions in the 
#'   soc.ca package
#' @param dim the dimensions in the order they are to be plotted. The first 
#'   number defines the horizontal axis and the second number defines the 
#'   vertical axis.
#' @param label if TRUE the label of the points are shown
#' @param min.size is the minimum size given to the groups of a numerical 
#'   variable, see \link{min_cut}.
#' @param ... further arguments are passed onto \link{geom_path},
#'   \link{geom_point} and \link{geom_text} from the ggplot2 package
#' @export
#' @examples
#' example(soc.ca)
#' map <- map.ind(result, point.color = as.numeric(sup$Age))
#' map <- map + scale_color_continuous(high = "red", low = "yellow")
#' map.path(result, sup$Age, map)

map.path  <- function(object, x, map = map.ind(object, dim), dim = c(1, 2),
                      label = TRUE, min.size = length(x)/10, ...){
  
  x.c     <- x
  if (is.numeric(x)) x.c     <- min_cut(x, min.size = min.size) 
  
  x.av    <- average.coord(object = object, x = x.c, dim = dim)  
  map.p   <- add.count(x.av, map, label, ...) 
  map.p
}