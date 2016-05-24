#'Get Contour Lines (data.frame)
#'
#'The following routine produces contour lines for a set of non-regular \code{x,y} and \code{z} values.
#'via utilizing a Deleaunay Mesh strung between the supplied \code{x,y} coordinates in order to produce 
#'iso-contour data representing the third variable, \code{z}. To this end, by using a Deleaunay mesh, 
#'this routine does not require regular \code{x} and \code{y} data, although it can be expected to yield
#''better' result, with regular / fine-grained data. 
#'
#'@param x,y Numeric data for x and y coordinate, a single matrix or data-frame object can be provided for 
#'\code{x}, which will be used in preference to the y and z arguments. These do \strong{NOT} need to be in any particular order
#'nor do they need to be regular.
#'@param z numeric Data for z coordinate (the coordinate to model)
#'@param nlevels An integer number of bins to split the data into \strong{iff} \code{levels} or \code{binwidth}
#'have not been specified.
#'@param binwidth The desired width of the bins, if specified, will override \code{nlevels}.
#'@param levels A numeric vector of the explicitly specified levels (z values) to contour, by specifying this argument, 
#'it will override \code{nlevels} and/or \code{binwidth}. If this argument is provided, the stacking order of the 
#'contours will be preserved in the order of first occurence within the supplied vector.
#'@param criticalRatio When producing the Delaunay Mesh, the quality of the mesh can be poor in the proximity
#'to the convex hull, Del's that have an aspect ratio greater than this value are not considered when producing the 
#'contours. In this context, the aspect ratio is defined as the circumradius to twice its inradius, 
#'equilateral triangles have an aspect ratio of 1, everything else is larger
#'@return For the function \code{getContourLines(...)}, the return object is a \code{data.frame} obect 
#'representing the contours, assembled in such way to permit easy use within the \code{ggplot2} paradigm.
#'Such data frame contains seven (7) columns: 
#'\item{\code{LID}}{A number representing the level}
#'\item{\code{GID}}{Within each level, a number representing the contour group}
#'\item{\code{PID}}{Within each group, a number representing the path/segment}
#'\item{\code{x}}{The x-coordinates of the contour}
#'\item{\code{y}}{The y-coordinates of the contour}
#'\item{\code{z}}{The z-coordinates of the contour (ie value of the level)}
#'\item{\code{Group}}{The unique identifyer for each independent contour path, calculated as being the 
#'interaction between \code{LID} and \code{GID}}
#'@seealso \code{\link{contourLinesR}}, \code{\link{getDelaunayMesh}} and \code{\link{getConvexHull}}
#'This is a wrapper to the \code{C++} interface function, \code{\link{contourWalker}}.
#'@examples
#'# Contour Lines for Volcano Data
#'library(ggplot2)
#'data(volcano)
#'x = 1:nrow(volcano)
#'y = 1:ncol(volcano)
#'z = expand.grid(x=x,y=y); z$z = apply(z,1,function(xx){ volcano[ xx[1],xx[2] ]} )
#'df = getContourLines(z)
#'ggplot(df,aes(x,y,group=Group,colour=z)) + geom_path()
#'    
#'# Contour Lines for a Function
#'library(ggplot2)
#'a      = -2; b = 2; n = 75
#'x  = y = seq(a,b,by=diff(c(a,b))/(n+1))
#'df     = expand.grid(x=x,y=y)
#'df$z   = with(df,-x*y*exp(-x^2-y^2)) 
#'df.cnt = getContourLines(df)
#'ggplot(data=df.cnt,aes(x,y,group=Group,colour=z)) + geom_path()
getContourLines <- function(x,y,z,nlevels=10,binwidth,levels,criticalRatio=5.0){
  bins = nlevels
  buildFrame = TRUE
  if(class(x) %in% c('data.frame','matrix')){
    if(ncol(x) == 3){
      inputData = data.frame(x); 
      buildFrame = FALSE;
    }
  }
  if(buildFrame && !exists('inputData')){ inputData = data.frame(x,y,z); }
  if(ncol(inputData) != 3){ stop("Problem resolving input xyz data") }
  colnames(inputData) = c("x","y","z")
  
  #Determine the range of z levels
  zRange = range(inputData$z)
  
  ##Determine the levels
  useExplicitLevels = FALSE
  if(!missing(levels)){
    if(is.numeric(levels)){
      if(any(levels >= zRange[1]) & any(levels <= zRange[2])){
        useExplicitLevels = TRUE  
      }else{
        warning(sprintf("Requested levels are outside the range [%E,%E]",zRange[1],zRange[2]))
      }
    }
  }
  
  #If Levels havn't explicitly specified
  if(!useExplicitLevels){
    if(!missing(binwidth)){
      binwidth = binwidth[1]
      if(is.numeric(binwidth) && binwidth != 0){
        bins = diff(zRange)/abs(binwidth)
      }
    }
    if(!is.numeric(bins))
      stop("'bins' is expected to be numeric")
    bins    = max(as.integer(round(bins)),1)
    step    = diff(zRange)/(bins+1)
    #levels  = seq(zRange[1]+step/2,zRange[2]-step/2, by = step)
    levels  = seq(zRange[1],zRange[2], by = step)
  }
  
  #Now Process for the given levels
  dm      = getDelaunayMesh(inputData$x,inputData$y) 
  
  #Solve the Contour
  result  = data.frame(contourWalker(dm,as.matrix(inputData),levels=levels,criticalRatio=criticalRatio));
  colnames(result) = c("LID","GID","PID","x","y","z")
  
  #Add the interaction between LID and GID, to Produce Group
  result$Group = paste(result[,1],result[,2],sep="-")
  result$Group = factor(result$Group,levels=unique(result$Group))
  #result$Group = interaction(result[,1],result[,2],sep="-")
  
  
  #Reorder and return
  return(invisible(result[,c("LID","GID","PID","Group","x","y","z")]));
}


#' Get Contour Lines (list)
#' 
#'A wrapper to the \code{\link{getContourLines}} function, provided to ease the transition from
#'the \code{\link[grDevices]{contourLines}} function as part of \code{\link{grDevices}}
#'
#'This function returns data in the same format/structure as \code{\link[grDevices]{contourLines}}, 
#'ie, list of lists, which is different from the (preferred) \code{dataframe} object returned
#'by the \code{\link{getContourLines}} function as part of the present work.
#' @inheritParams getContourLines
#' @param ... any other parameters passed through to \code{\link{getContourLines}}
#' @seealso \code{\link{getContourLines}}
#' @return
#' A list of contours is returned, Each contour is a list with the elements:
#'\item{\code{level}}{The contour of the level}
#'\item{\code{x}}{The x-coordinates of the contour}
#'\item{\code{y}}{The y-coordinates of the contour}
#' @examples
#' library(contoureR)
#' library(ggplot2)
#' x = runif(100)
#' y = runif(100)
#' df = expand.grid(x=x,y=y)
#' z  = with(df,x+y)
#' result = contourLinesR(df$x,df$y,z)
#' @aliases contourLinesR
contourLinesR <- function(x,y,z,nlevels=10,levels=pretty(range(z, na.rm = TRUE)),...){
  result = getContourLines(x,y,z,nlevels,levels,...)
  result = dlply(result,c("Group"),function(df){
    r = list(); 
    r$level = df$z[1]; r$x = df$x; r$y = df$y
    r
  })
  names(result) <- NULL
  attributes(result) <- NULL
  return(result)
}
