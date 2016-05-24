#################################################################################
#
#  plot methods...
#
#  Right now, because of inheritance, importance and cv will use this one method
#  for cmc defined below.
#
#  There is no plot method for antitheticSampling at the moment.
#
#Author...									Date: 7-May-2013
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#



#================================================================================
#   method for class "crudeMonteCarlo"...
#
setMethod('plot',
          signature(x = 'crudeMonteCarlo', y='missing'),
function(x,
         axes = FALSE,
         renderAs = c('profile', 'crossSection'),
         isHeightColor = .StemEnv$isHeightColor,
         ...                                      #for plot
        )
{
#------------------------------------------------------------------------------
#
#   Things are kept simple here, since we are always dealing with only one
#   "Stem" we center it at (0,0) for drawing. 
#
#------------------------------------------------------------------------------
#
    renderAs = match.arg(renderAs)

    taper = x@stem@taper
    taperNames = names(taper)

#
#   for 'profile' always render as a downLog object regardless of what it really is;
#   for 'crossSection', always render as a standingTree...
#
    if(renderAs == 'profile') {
      names(taper) = c(taperNames[1], 'length')
      ypos = max(taper$length)/2
      stem = downLog(taper, solidType=x@stem@solidType, centerOffset = c(x=0,y=ypos), logAngle = pi/2)
      plot(stem, axes=axes, ...)
      rad = x@diam.s/2
      SL = list()
      for(i in seq_along(rad)) {
        line = cbind(c(-rad[i], rad[i]), c(x@hgt.s[i], x@hgt.s[i]))
        L = Line(line)
        SL = c( SL, Lines(list(L), ID=paste('is',i,sep='.')) )
      }
      IS = SpatialLines(SL)
      suppressWarnings({                               #for non-plot arguments in ...    
        plot(IS, add=TRUE, col=isHeightColor, ...)     #e.g., add
      })   
    }
    else {                                             #cross-sectional
      names(taper) = c(taperNames[1], 'height')
      stem = standingTree(taper, solidType=x@stem@solidType, centerOffset = c(x=0,y=0))
      stem@dbh = stem@buttDiam                         #not kosher, but it makes buttDiam render rather than dbh
      stem@spDBH = spCircle(stem@buttDiam/2)$spCircle  #and it is not used after this
      plot(stem, axes=axes, ...)
      rad = x@diam.s/2
      for(i in seq_along(rad)) {
        cir = spCircle(rad[i])$spCircle
        suppressWarnings({                                #for non-plot arguments in ...    
          plot(cir, add=TRUE, border=isHeightColor, ...)  #e.g., border
        })   
      }
      #top diameter...
      topRad = x@stem@topDiam/2
      if(!identical(topRad, 0)) {
        suppressWarnings({                                #for non-plot arguments in ...    
          plot(spCircle(topRad)$spCircle, add=TRUE, col=.StemEnv$treeColor,
               border = .StemEnv$treeBorderColor, lty='dashed', ...)  #just make it a little darker
        })
      }
    }

        
    return(invisible())

}    #plot for 'crudeMonteCarlo'
) #setMethod
    

