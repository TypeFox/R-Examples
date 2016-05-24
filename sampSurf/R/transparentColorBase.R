transparentColorBase = function(color, alphaTrans=alphaTrans) 
{
#---------------------------------------------------------------------------
#
#   Adding transparency in base graphics is a little more difficult than
#   in lattice. This function will take a color and alpha transparency
#   value and return the correct rgb color to be used in plot(),
#   image() and the like.
#
#   Arguments...
#     color = the color desired
#     alphaTrans = alpha transparency value [0,1]
#
#   Returns...
#     the rgb value
#
#Author...									Date: 5-May-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#
#---------------------------------------------------------------------------
#
    if(alphaTrans>1 || alphaTrans<0)
      stop(paste('alphaTrans (',alphaTrans,') must be in [0,1]!', sep=''))
    if( any(is.na(match(color, colors()))) )
      stop(paste('color (',color,') must be legal R colors()!',sep=''))
    
    ac = t(col2rgb(color))/255
    return(rgb(ac[,1], ac[,2], ac[,3], alpha=alphaTrans))
}   #transparentColorBase
