#---------------------------------------------------------------------------
#
#   This file holds the S4 class definitions for the Tract & related classes.
#
#   Classes...
#     1. Tract: 
#     2. bufferedTract: 
#
#Author...									Date: 31-Aug-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#




#=================================================================================================
#
#  define the Tract class...
#
setClass('Tract',
         
#
#  slots for the class and its subclasses...
#
    representation(description = 'character',      #more descriptive name
                   units = 'character',            #English or metric units
                   area = 'numeric'                #area in square units
                  ),
         
    prototype = list(description = '',
                     units = .StemEnv$msrUnits$metric,
                     area = 0
                    ),
         
    contains = 'RasterLayer',                      #subclass of RasterLayer
         
    validity = function(object) {
                 if(!(object@units %in% c('English','metric')))
                   return('units of measure must be "English" or "metric"')

                # if(!is.na(object@spUnits@projargs) && object@spUnits@projargs == '+proj=longlat')
                #   return(paste('spUnits must be commensurate with units,',
                #                'please convert to non-geographic coordinate system!')
                #         )
                 #if(!is.na(projection(object)) && projection(object) == '+proj=longlat')
                 ##gp = grep(projection(object), '+proj=longlat') #projargs can have multiple parts
                 ##if(!is.na(projection(object)) && length(gp) > 0)
                 if(!is.na(projection(object)) && isLonLat(object) > 0)
                   return(paste('CRS must be commensurate with units,',
                                'please convert to non-geographic coordinate system!')
                         )

                # if(object@units=='English' && !is.na(object@spUnits@projargs))
                #   return('English units are not compatible with metric projections!')
                 
                 #if(object@units=='English' && projection(object)!="NA" )
                 if(object@units=='English' && !is.na(projection(object)) )
                   return('English units are not compatible with metric projections!')

                 
                 return(TRUE)
               } #validity check
) #class Tract
 





#=================================================================================================
#
#  the bufferedTract class is a direct descendant of 'Tract'...
#
#
setClass('bufferedTract',
    representation(bufferRect = 'matrix',              #holds the buffer in matrix form
                   spBuffer = 'SpatialPolygons'
                  ),
         
    prototype = list(bufferRect = matrix(NA, nrow=2, ncol=2,
                                         dimnames=list(c('x','y'), c('min','max'))),
                     spBuffer = new('SpatialPolygons')
                    ),
         
    contains = 'Tract',                      #subclass of Tract
         
    validity = function(object) {
                 if(any(is.na(object@bufferRect)))
                   return('bufferRect can not have missing values!')

                 #essentially the same checks as in bboxCheck()...
                 if(!class(object@bufferRect) == 'matrix')
                   return('bufferRect slot must be a 2x2 matrix')
                 bboxNames = match(rownames(object@bufferRect), c('x','y'))
                 if(any(is.na(bboxNames)))
                   return('slot bufferRect rownames must be "x", "y"!')
                 bboxNames = match(colnames(object@bufferRect), c('min','max'))
                 if(any(is.na(bboxNames)))
                   return('slot bufferRect colnames must be "min", "max"!')
                 if(any( apply(object@bufferRect,1,function(x) if(x['min'] >= x['max']) TRUE else FALSE) ))
                   return('in slot bufferRect, "min" must be less than "max" for x and y!')
                 
                 
                 return(TRUE)
               } #validity check
) #class bufferedTract
