#' GCD sites interactive selection
#' 
#' Interactive selection of GCD sites by drawing a polygon on a map.
#' 
#' 
#' @param addata An optional XY matrix of coordinates to specify a polygon to
#' be drawn on the map.
#' @return An object of the class "pfSiteSel".
#' @author O. Blarquez
#' @seealso \code{\link{pfSiteSel}}
#' @examples
#' 
#' \dontrun{
#' # Type: 
#' ID=pfInteractive()
#' # And follow text instructions}
#' 
#' 
pfInteractive=function(addata=NULL){
  
  # install.packages("Imap");
  ## Avoid no visible binding for global variable
  paleofiresites=NULL; rm(paleofiresites)
  countries=NULL; rm(countries)
  
  if (!requireNamespace("Imap", quietly = TRUE)) {
    install.packages("Imap")
  }
  
  ## Load data  
  data(paleofiresites,envir = environment())
  data(countries,envir = environment())
  ## Define vectors
  yy=cbind(countries$x,countries$y)
  pp=cbind(paleofiresites$long,paleofiresites$lat)
  
  ## Use imap for interactive plot
  cat("You can zoom the map by left clicking in two different locations
on the figure (upper left and bottom right corners) to define a rectangle that will be zoomed into. Left
clicking outside the plot region (but somewhere in the figure region)
will zoom out. Double left clicking on the same spot will reset the plot.
Rignt click or click finish when ready.")
  plot(pp,bg="blue",col = "black",pch = 21,ylim=c(-90,90),xlim=c(-180,180), ylab="Latitude",xlab="Longitude")
  par(bg="white")
  Imap::imap(yy,fill=F,zoom=T,col="black",poly = rgb(238/255,220/255,130/255),add.all=T,grid=T)
  points(pp,bg="blue",col = "black",pch = 21,ylim=c(-90,90))
  if (is.matrix(addata)){
    lines(addata[,1],addata[,2])
  }
  ## Use select.pts for selection
  cat("\n")
  cat("\nInteractively select points by drawing a polygon composed of at least 
three vertices (polygon must be convex ). Right click or click finish when ready.")
  a=Imap::select.pts(pp)
  
  
  ## Plot selected points
  if (length(a)>2){
    plot(pp,bg="blue",col = "black",pch = 21,xlim=c(min(a[,1])-10,max(a[,1])+10),ylim=c(min(a[,2])-10,max(a[,2])+10))
    ## Retrive site IDs ans site names
    IDs=paleofiresites[paleofiresites$lat %in% a[,2] & paleofiresites$long %in% a[,1], 1]
    points(a,bg="red",col = "black",pch = 21)
    
  }
  if (length(a)==2) {plot(pp,bg="blue",col = "black",pch = 21,xlim=c((a[1])-10,(a[1])+10),ylim=c((a[2])-10,(a[2])+10))
                     ## Retrive site IDs ans site names
                     IDs=paleofiresites[paleofiresites$lat %in% a[2] & paleofiresites$long %in% a[1], 1]
                     points(a[1],a[2],bg="red",col = "black",pch = 21)
                     
  }
  lines(yy)
  if (is.matrix(addata)){
    lines(addata[,1],addata[,2])
  }
  
  ## Retrive site IDs ans site names
  
  site_name=as.character(paleofiresites$site_name[paleofiresites$id_site%in% IDs])
  ## Return output ID list
  output=list(id_site=IDs,site_name=site_name)
  class(output)="pfSiteSel"
  ## Remove data
  #rm(paleofiresites,envir = globalenv())
  #rm(countries,envir = globalenv())
  ## Return output
  return(output)
  
}
