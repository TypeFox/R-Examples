#---------------------------------------------------------------------------
#
#   This generic and method will take and "InclusionZoneGrid" object and a
#    'Tract' object and "heap" or accumulate the inclusion zone into the tract.
#
#   There is need for only one method...
#     1. "IncluzionZoneGrid" & "Tract" signatures, or any of their subclasses.
#
#Author...									Date: 9-Sept-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#   generic definition...
#
setGeneric('heapIZ',  
           function(izgObject, tract, ...) standardGeneric('heapIZ'),
           signature = c('izgObject', 'tract')
          )



  


          
#================================================================================
#  method for 'InclusionZoneGrid' and 'Tract' classes...
#
setMethod('heapIZ',
          signature(izgObject = 'InclusionZoneGrid', tract='Tract'),
function(izgObject,
         tract,
         estimate = unlist(c(.StemEnv$puaEstimates, .StemEnv$ppEstimates)),
         ...
        )
{
#---------------------------------------------------------------------------
#
#   get the estimate desired, then fill the grid cells with this attribute
#   and extend/expand the subgrid within the izgObject to the full tract size...
#
    estimate = match.arg(estimate)
    grid = setValues(izgObject@grid, izgObject@data[,estimate]) #set desired attribute
    grid = extend(grid, tract)                                  #new cells are assigned NA values
    gv = getValues(grid)                                        #so we must change these
    k = ifelse(is.na(gv), 0, gv)                                #to zero for subsequent accumulation
    grid = setValues(grid, k)

#    
#   accumulate into the original tract object and return...
#
    tmp = grid + tract                          #adding changes the class to 'RasterLayer' so use tmp
    tract = setValues(tract, getValues(tmp) )   #still a "Tract" subclass

    return(tract)
}   #heapIZ for'InclusionZoneGridIZ'
)   #setMethod




