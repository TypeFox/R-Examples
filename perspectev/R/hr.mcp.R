

hr.mcp <-
function(x,y=NULL,n.min=50,plot=TRUE,add=FALSE,ID=NULL,...){
  if(is.null(y))xy <- x
  if(!is.null(y))xy <- cbind(x,y)
  if(inherits(x,what="data.frame"))xy <- as(xy,"matrix")
  #Missing values allowed but dropped and trigger a warning 
  xy. <- xy
  xy <- na.omit(xy)
  if(nrow(xy.)!=nrow(xy))warning("Missing values dropped")
  #IF XY DOES NOT SUPPORT ESTIMATION, STOP HERE
  if(nrow(xy)<n.min){
    warning("Fewer than 'n.min' qualifying observations")
    return(NA)
  }
  mcp <- xy[chull(xy),]
  mcp <- rbind(mcp,mcp[1,])
  mcp <- Polygon(mcp,hole=FALSE)
  if(is.null(ID))ID <- as.logical(NA)
  mcp <- Polygons(list(mcp),ID=ID)
  mcp
}