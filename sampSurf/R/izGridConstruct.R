izGridConstruct = function(izObject,
                           tract,
                           description = 'sausageIZ inclusion zone grid object',
                           wholeIZ = TRUE,           #TRUE: grid the whole object; FALSE: just grid the IZ
                           ...
                          )
{         
#---------------------------------------------------------------------------
#
#   This is the base constructor for the InclusionZoneGrid class that should 
#   be used for the majority of izObject objects (e.g., sausageIZ, standUpIZ,...);
#   the exception being chainSawIZ objects, whose IZs are really just a single grid
#   point.
#
#   Arguments...
#     izObject = an object of class InclusionZone or subclass; but it can not
#                be an object of class chainSawIZ
#     tract = an object of class tract of subclass
#     description = a comment
#     wholeIZ = TRUE: make the bounding grid contain not only the IZ, but the
#               entire object (i.e., including the log); FALSE: the grid
#               will only contain the inclusion zone, the log can lie outside
#               of the bounding grid in this case (i.e., standUpIZ)
#     ... = to be gobbled or passed on
#
#   Returns...
#    An object of class InclusionZoneGrid
#
#Author...									Date: 17-Sept-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   unfortunately, we need these checks since it is not a method...
#
    if(!is(izObject, 'InclusionZone') || !is(tract,'Tract'))
      stop('must pass valid inclusion zone and tract objects!')
    vo = validObject(izObject)
    vo = validObject(tract)
    if(is(izObject, 'chainSawIZ'))
      stop('please use izGrid() method for objects of class chainSawIZ!')
    if(izObject@units != tract@units)
      stop('Inclusion zone and tract units are not compatible!')

    
    if(!wholeIZ) 
      iz.bbox = bbox(perimeter(izObject))              #IZ only
    else
      iz.bbox = bbox(izObject)                         #IZ plus log (always internal)    
    izg = izGrid(iz.bbox, tract)                       #get the minimal bounding grid

#
#   grid with NAs outside (background) the overlay (inclusion zone) region and zeros inside...
#
    #grid = polygonsToRaster(perimeter(izObject), izg, mask=TRUE, silent=TRUE)  ##replaced by rasterize...
    grid = rasterize(perimeter(izObject), izg, mask=TRUE, silent=TRUE)

#
#   a data frame with each pua estimate; now we make the background 0 (rather than NA) and
#   fill the inclusion zone with the appropriate estimate (which can be NA if it does not
#   exist for a particular method)...
#
    nr = ncell(grid)
    npua = length(izObject@puaEstimates)
    df = data.frame(matrix(NA, nrow = nr, ncol = npua))
    colnames(df) = names(izObject@puaEstimates)

    gridVals = getValues(grid)
    for(i in seq_len(npua)) 
      df[,i] = ifelse(is.na(gridVals), 0, izObject@puaEstimates[[i]])

    df$depth = ifelse(is.na(gridVals), 0, 1)   #sample or overlap zone depth   #****

#
#   combine them for the overall bbox...
#
    grid.bbox = bbox(grid)
    iz.bbox = bbox(izObject)                        #must always be whole object extents here
    min = apply(cbind(grid.bbox, iz.bbox), 1, min)
    max = apply(cbind(grid.bbox, iz.bbox), 1, max)
    bbox = matrix(cbind(min,max),ncol=2, dimnames=list(c('x','y'), c('min','max')))
    

    griz = new('InclusionZoneGrid', description = description,
               iz = izObject,
               grid = grid,
               data = df,
               bbox = bbox
              )
    return(griz)
}   #izGridConstruct
