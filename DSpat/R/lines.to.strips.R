lines_to_strips = function(lines, study.area, width=NULL)
################################################################################
#
# Convert lines with transect widths to strips and compute rotation angle
# to be able to rotate it to vertical.
#
# Arguments:
#    lines      - dataframe with names label,x0,x1,y0,y1 and optionally width
#    study.area - owin class giving study area window
#    width      - if all lines have the same width it can be specified here
#
# Value: a list with elements
#    lines     - a psp class of lines with label and angles for rotation added
#    transects - a list of dataframes of polygon coordinates
#
# 7 April 2008
# Jeff Laake
################################################################################
{
# Define function that does the work for a single line
 line_to_strip=function(x,study.area,clip=TRUE)
 {
# Compute length of line, slope and set up x and y values for vertices
   d=sqrt((x[1]-x[2])^2+(x[3]-x[4])^2)
   slope=(x[4]-x[3])/(x[2]-x[1])
#  set up y and x points
   yp=c(rep(x[3],2),rep(x[4],2))
   xp=c(rep(x[1],2),rep(x[2],2))
   w=x[5]/2
#  As long as line is not vertical compute angle and the vertices of the polygon
   if(!is.infinite(slope))
   {
      theta=pi/2-asin((x[4]-x[3])/d)
      deltay=w*sin(theta)
#     Vertices of strip must be listed so they are counter-clockwise
#     This depends on whether the line is going from left to right or right to left
      if(x[1]<x[2])
         y=yp+c(deltay,-deltay,-deltay,deltay)
      else
      {
         y=yp+c(-deltay,deltay,deltay,-deltay)
         theta=-theta
      }
      x=-slope*(y-yp)+xp
   }
# If vertical, the computation is straightforward except that you still
# have the 2 cases depending on direction.
   else
   {
      if(x[3]<x[4])
      {
         theta=0
         x=xp+c(-w,w,w,-w)
      }
      else
      {
         theta=pi
         x=xp+c(w,-w,-w,w)
      }
      y=yp
   }
   if(clip)
   {
      gpc.area=owin.gpc.poly(study.area)
      b=as(data.frame(x=x,y=y),"gpc.poly")
      inside.poly=get.pts(intersect(b,gpc.area))[[1]]
      xdf= data.frame(x=rev(inside.poly$x),y=rev(inside.poly$y))
      return(list(poly=xdf[!duplicated(xdf),],angle=theta))
   }
   else
      return(list(poly=data.frame(x=x,y=y),angle=theta))
 }
# create the lns matrix with line begin/end points and strip width
 if(!is.null(width))
   lns = cbind(lines[,c('x0','x1','y0','y1')],rep(width,dim(lines)[1]))
 else
   if("width" %in% colnames(lines))
     lns = lines[,c('x0','x1','y0','y1','width')]
   else
     stop("Missing either width specification or width column in lines matrix")
# call line_to_strip for each line and return the
 lines.psp <- psp(x0=lines$x0,y0=lines$y0,x1=lines$x1,y1=lines$y1,window=study.area,check=FALSE)
 lines.psp$width=lns[,5]
 lines.psp$label=lines[,1]
 xlns=apply(lns,1,line_to_strip,study.area=study.area)
 angles=unlist(lapply(xlns,function(x)x$angle))
 names(angles)=NULL
 lines.psp$angles=angles
 transects=lapply(xlns,function(x)x$poly)
# get full transects at max width
 lns[,5]=max(lns[,5])
 xlns=apply(lns,1,line_to_strip,study.area=study.area,clip=FALSE)
 full.transects=lapply(xlns,function(x)x$poly)
 return(list(lines=lines.psp,transects=transects,full.transects=full.transects))
}


