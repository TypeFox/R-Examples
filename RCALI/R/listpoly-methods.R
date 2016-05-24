# +++++++++++++++++++++++++++++++++++++++++++++++++++++++
# AB: oct 2005, sept 2009
# ++++++++++++++++++++++++++++++++++++++++++++++++++++
# -----------------------------------------------------
# crlistpoly:
# FUNCTION:
  # Create an object of classe listpoly from objects "poly"
# -----------------------------------------------------
crlistpoly <- function(...) {
  retour<- list(...)
  for (i in retour) {
    if (class(i) != "poly")
      stop("Arguments should be objects of class 'poly'")
  }
  class(retour) <- "listpoly"
  return(retour)
} # end crlistpoly
        

# ++++++++++++++++++++++++++++++++++++++++++++++++++++
# METHODS 'listpoly'
# ++++++++++++++++++++++++++++++++++++++++++++++++++++
# Creator: NULL
# -----------------------------------------------------
listpoly <- function() {
  return(NULL)
}

# -----------------------------------------------------
# range.listpoly: 
# FUNCTION:
# Determine the min and max of the coordinates of all the polygons
# of an object of classe listpoly by means of my function 'range.poly'
# ARGUMENTS:
# - x: an object of classe listpoly
# - ...: a variable list of arguments which will be passed as it to
#  'range.poly'
# -----------------------------------------------------
range.listpoly <- function(x, ...) {
  xrange <- yrange <- c(NULL, NULL)

  for (l in x) {
    xrange <- range(xrange, l[,"xcoord"],...)
    yrange <- range(yrange, l[,"ycoord"],...)
  }
  retour <- cbind(xrange,yrange)
  dimnames(retour) <- list(c("lower","upper"),c("xrange","yrange"))
  return(retour)
}

# -----------------------------------------------------
# Plotting functions
# -----------------------------------------------------

# -----------------------------------------------------
# plot.listpoly:
# FUNCTION:
# Plot all the 'poly' of a 'listpoly' on the same frame, 
# each one in a distinct color, by means of my function 'plot.poly'
# ARGUMENTS:
# - x: an object  of classe listpoly
# - add: when TRUE, the frame is not drawn
#  (for use after a zoom, i.e)
# - ...: a variable list of arguments which will be passed as it to
#  'plot.poly'
# NOTE:
# When a graphic device is already opened, the plot is drawn on it;
# Otherwise, X11 device is opened
# Example:
# plot(parc)
# zoom() -> click on 2 points
# plot(parc, add=T)
# -----------------------------------------------------
plot.listpoly <- function(x, add=F, color=T, ...) {
listpoly=x

## Open X11 device, in adding colors,
# when no device is already opened en rajoutant des couleurs

# To center in zero
  l <- range(listpoly)



# 12/11/2007: lapply ne marche plus
# Remplacement de ceci par la boucle qui suit:
# ###
# ###  listpoly <- lapply(x, function(a) {
# ### a[,1]=a[,1]-l["lower","xrange"]; 
# ### a[,2]=a[,2]-l["lower","yrange"]; 
# ### return(a)
# ### })

listpoly <-x
for (xl in 1:length(x))
{
listpoly[[xl]][,1]= listpoly[[xl]][,1]-l["lower","xrange"];
listpoly[[xl]][,2]= listpoly[[xl]][,2]-l["lower","yrange"];
}




  class(listpoly) <- "listpoly"
if (add == FALSE) {
# To have the same scale on both axes:
  l <- range(listpoly)
  etendue <- max(apply(l,2,diff))

# Draw the frame: (pty='s' => a square)
  par(pty="s")
  plot( x=c(l[1,"xrange"], l[1,"xrange"]+etendue),
       y=c(l[1,"yrange"], l[1,"yrange"]+etendue),
       type="n", xlab="", ylab="")
} # end (add == FALSE)

if(color==TRUE) {
# Create the colors:
  lg <- length(listpoly)
  r <- rainbow(lg)
  for (i in 1:lg) 
    attr( listpoly[[i]], "couleur") <- r[i]
}

# Loop over the polygons:
# 12/11/2007: lapply ne marche plus
# Remplacement de ceci par la boucle qui suit:
# ###  lapply(listpoly, plot.poly, add=T, ...)
for (xl in 1:length(listpoly))
    plot.poly(listpoly[[xl]], add=T, ...)


# Add the number of the poly:
  for (i in 1:length(listpoly)) {
    text(x=mean(listpoly[[i]][,1]),
	 y=mean(listpoly[[i]][,2]),
	 labels=names(listpoly[i]))
  }

  invisible()
} # end plot
# -----------------------------------------------------
# export:
# FUNCTION:
# Write on file a 'listpoly' according to the requested format.
# For each polygon, 2 lines:
# identificator, followed by the x-coordinates
# identificator, followed by the y-coordinates
# ---------------------------------------------------------

export <- function(x, filename)
  UseMethod("export")

export.default <- function(x,filename)
  {
  if (!inherits(x,"listpoly"))
    stop("This function is valid only on a 'listpoly' object.")
  }
# -----------------------------------
export.listpoly <- function(x, filename) {
  listpoly<-x 
  l <- length(listpoly)
  sink(filename)
  cat(l)
  cat("\n")
  for (i in 1:l) {
    cat(i); cat("  ");cat(listpoly[[i]][,1]); cat("\n");
    cat(i); cat("  ");cat(listpoly[[i]][,2]); cat("\n");
  }
  sink()
  invisible()
}

