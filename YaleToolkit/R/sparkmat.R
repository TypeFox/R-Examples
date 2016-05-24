############################################################
#################### 60 Characters Wide ####################
############################################################
# Walton A. Green
# Department of Geology
# Yale Unversity
# P. O. Box 208109 Yale Station
# New Haven, Connecticut 06520
# walton.green@yale.edu
#
# John W. Emerson
#
############################################################
#################### 60 Characters Wide ####################
############################################################
# 
#  plots panels of sparklines at arbitrary
#   locations (all panels must have data of the same
#   shape; multiple calls are needed if the shape of the
#   data varies).
#
############################################################
# sparkmat()
############################################################

sparkmat <- function(x,
                     locs = NULL,
                     w = NULL,
		     h = NULL,
                     lcol = NULL,
                     yscales = NULL,
		     tile.shading = NULL,
                     tile.margin = unit(c(0,0,0,0), 'points'),
                     tile.pars = NULL,
                     just = c('right', 'top'),
                     new = TRUE,
                     ...) {

  if (new) grid.newpage()

### Walton, why would x[[1]] ever be NULL?  Trap elsewhere?
###
  if (!is.null(x[[1]]) && is.null(yscales)) {
    yscales <- vector(mode="list", length=length(x[[1]]))
    for (i in 1:length(x)) { # Over the panels
      for (j in 1:length(x[[1]])) { # Over the sparklines in a panel
        yscales[[j]] <- c(min(yscales[[j]][1], min(x[[i]][,j], na.rm=TRUE)),
                          max(yscales[[j]][2], max(x[[i]][,j], na.rm=TRUE)))
      }
    }
  }

  # A little function that takes two vectors of, for instance, the x and y
  #  values at which a regular grid of points should be drawn and
  #  returns a data frame with the ordered pairs at which to plot all the
  #  points as two columes of a data frame.
  vectorize <- function(x,y){
    x.v <- rep(x, length(y))
    y.v <- as.numeric(matrix(y, nrow = length(x),
                             ncol = length(y), byrow = TRUE))
    return(data.frame(x = x.v, y = y.v))
  }

  # if locs is NULL, set up a regular grid for the sparkmats
  if (is.null(locs)) {
    mats.down <- floor(sqrt(length(x)))
    mats.across <- ceiling(length(x) / mats.down)
    locs <- vectorize(x = (1:mats.across) / mats.across,
	              y = (mats.down:1) / mats.down)
    locs$x <- unit(locs$x, 'npc')
    locs$y <- unit(locs$y, 'npc')
    if (is.null(w)) w <- unit(1/mats.across, 'npc')
    if (is.null(h)) h <- unit(1/mats.down, 'npc')
  } else {       ############ Added by Jay 11/3/06
    if (new) {
      pushViewport(viewport(x=0.15, y=0.1, width=0.75, height=0.75, 
                            just=c("left", "bottom"),
                            xscale=range(pretty(locs[,1])),
                            yscale=range(pretty(locs[,2]))))
      grid.xaxis()
      grid.yaxis()
    }
  }

  if (!is.unit(w)) w <- unit(w, "native")  
  if (!is.unit(h)) h <- unit(h, "native")
  
  for (i in 1:length(x)) {
    if (is.unit(locs[i,1])) xloc <- locs[i,1] 
    else xloc <- unit(locs[i,1], "native")
    if (is.unit(locs[i,2])) yloc <- locs[i,2] 
    else yloc <- unit(locs[i,2], "native")
			
    sparklines.viewport <- viewport(x=xloc, y=yloc,
                                    just=just, width=w, height=h)
    pushViewport(sparklines.viewport)
    if (!is.null(tile.pars)) grid.rect(gp=tile.pars)
    sparklines(x[[i]], new=FALSE, lcol=lcol, yscale=yscales,
               outer.margin=tile.margin,
               outer.margin.pars = gpar(fill = tile.shading[i], 
               col=tile.shading[i]), xaxis = FALSE, yaxis = FALSE)		      
    popViewport(1)
  }

} # End of function
