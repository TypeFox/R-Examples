"create.lines"<-function(study.area,nlines=NULL,width,spacing=NULL,angle=0)
################################################################################
# Create a systematic set of lines to sample a rectangular grid.  The systematic
# grid can be set at any angle and the number of lines is set by the spacing
# or the spacing is set by width and number of lines.
#
# Arguments:
#
#  study.area - owin class defining area
#  nlines     - number of lines
#  width      - full transect width
#  spacing    - spacing distance between centerlines
#  angle      - angle of rotation in degrees anticlockwise from x-axis
#
# Value: line dataframe with label,x0,y0,x1,y1,width where x0,y0 is beginning
#        and x1,y1 is end of the line
#
# Jeff Laake
# 10 April 2008
################################################################################
{
#  To fit the transects into study.area, a buffer must be created if the lines
#  are anything other than vertical or horizontal
   xr=study.area$xrange
   yr=study.area$yrange
   W=study.area
# Next work out the line spacing
   if(is.null(nlines))
   {
      if(is.null(spacing))
      {
        stop("\n Either nlines or spacing must be specified")
      }
      else
      {
         if(spacing <= width)
            stop("\n spacing must exceed the transect width")
      }
   }
   else
   {
      if((xr[2]-xr[1])/nlines < width)
        stop("\n too many lines for the specified width")
      else
      {
         if(angle==90 | angle ==270 )
            spacing=(xr[2]-xr[1])/nlines
         else
            if(angle==180 | angle==0 )
               spacing=(yr[2]-yr[1])/nlines
            else
               spacing=sqrt(diff(yr)^2+diff(xr)^2)/nlines
      }
        
   }
   xp=rlinegrid(angle=angle,spacing=spacing,win=W)
   xlines=data.frame(label=1:xp$n,x0= xp$ends$x0,x1=xp$ends$x1,y0=xp$ends$y0,
                       y1=xp$ends$y1,width=rep(width,xp$n))
   return(xlines)
}
