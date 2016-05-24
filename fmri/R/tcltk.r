# modified tkrplot from file tkrplot
# the parent of the plot depends on type of the plot (2/3 dimensional) and the situation it is called in (extract images
# or plot)
mytkrplot <- function(parent, fun, hscale=1, vscale=1,slicenr=-1,which="unset",parent2=list(),typeV=2,overview=FALSE) {
    if (!overview) if (typeV == 2) parent = parent2  
    # ACHTUNG aenderung, bisher nur in vied2dNew benutzt!!
    image <- paste("Rplot", .make.tkindex(), sep="")
    .my.tkdev(hscale, vscale)
    try(fun())
   .Tcl(paste("image create Rplot", image))
    lab<-tklabel(parent,image=image) #**** use try, delete image on failure
   tkbind(lab,"<Destroy>", function() .Tcl(paste("image delete", image)))
    lab$image <- image
    lab$fun <- fun
    lab$hscale <- hscale
    lab$vscale <- vscale
    lab
}
helpFunc <- function(a){
  a <- tclVar()
}



