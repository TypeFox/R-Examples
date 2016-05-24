#---------------------------------------------------------------------------
#
#   Contains all of the classes and methods necessary for using the mirage
#   method of boundary overlap correction in sampSurf
#
#   A.1 "mirageTract" class
#   A.2 "mirageTract" generic definition
#   A.3 "mirageTract" constructor method
#
#   B.1 "mirageInclusionZoneGrid" class definition
#   B.2 "mirageInclusionZoneGrid" generic definition
#   B.3 "mirageInclusionZoneGrid" constructor method
#
#   C.1 plot method for "mirageInclusionZoneGrid" objects (20-Nov-2013)
#
#Author...									Date: 12-Sept-2013
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
# A.1 the mirageTract class is a direct descendant of 'Tract'...
#
#
setClass('mirageTract',
          contains = 'Tract'                      #subclass of Tract
) #class mirageTract





          
#================================================================================
#================================================================================
# A.2 mirageTract generic...
#
#
setGeneric('mirageTract',  
           function(tract, ...) standardGeneric('mirageTract'),
           signature = c('tract')
          )

#--------------------------------------------------------------------------------
#
# A.3 mirageTract constructor method for class mirageTract...
#
#

setMethod('mirageTract',
          signature(tract = 'Tract'),
function(tract,
         ...
        )
{
#---------------------------------------------------------------------------
#
#   just cast, no new slot assignments...
#
    mt = as(tract, 'mirageTract')

    if(validObject(mt))
      return(invisible(mt))
    else
      return(NULL)

}   #mirageTract constructor
)   #setMethod
    









#================================================================================
#
#  B.1 "mirageInclusionZoneGrid" class definition as an extension to izGrid...
#
#================================================================================
#
setClassUnion('SPNULL', c('SpatialPolygons','NULL'))
setClassUnion('RLNULL', c('RasterLayer','NULL'))
setClassUnion('izgNULL', c('InclusionZoneGrid','NULL')) #for izGrid.extended

setClass('mirageInclusionZoneGrid',
    representation(slopOver = 'logical',       #TRUE/FALSE in .StemEnv$cardinal directions
                   north.polygon = 'SPNULL',   #SP sliver or NULL
                   north.grid = 'RLNULL',      #RasterLayer grid for sliver or NULL
                   south.polygon = 'SPNULL',
                   south.grid = 'RLNULL',
                   east.polygon = 'SPNULL',
                   east.grid = 'RLNULL',
                   west.polygon = 'SPNULL',
                   west.grid = 'RLNULL',
                   izGrid.extended = 'izgNULL' #holds the extended-tract version 
                  ),
    contains = 'InclusionZoneGrid',
    #below object won't really be valid, but it will assign NULLs as default for, e.g.,
    #IZs that have no overlap...
    prototype = list(description = 'mirageInclusionZoneGrid object',
                     slopOver = vector('logical', length(.StemEnv$cardinal)), #all FALSE by default
                     north.polygon = NULL,
                     north.grid = NULL,
                     south.polygon = NULL,
                     south.grid = NULL,
                     east.polygon = NULL,
                     east.grid = NULL,
                     west.polygon = NULL,
                     west.grid = NULL,
                     izGrid.extended = NULL, 
                     bbox = matrix(rep(0,4), nrow=2, dimnames=list(c('x','y'), c('min','max'))),
                     data = data.frame(matrix(NA,
                                              nrow = 0,
                                              ncol = length(c(.StemEnv$puaEstimates,.StemEnv$ppEstimates)),
                                              dimnames = list(character(0),
                                                         names(c(.StemEnv$puaEstimates,.StemEnv$ppEstimates)))
                                             ) #matrix
                                      ) #df
                    ),
    sealed = TRUE,                           #no further changes or subclasses
    validity = function(object) {
                 #nothing really to check here since all new slots can be NA or NULL
                 return(TRUE)
               } #validity check
) #class mirageInclusionZoneGrid 








#================================================================================
#  mirageInclusionZoneGrid constructor for 'InclusionZone' and 'mirageTract'
#  classes: generic and method...
#--------------------------------------------------------------------------------
#
#  B.2 generic...
#
setGeneric('izGridMirage',  
           function(izObject, tract, ...) standardGeneric('izGridMirage'),
             signature = c('izObject', 'tract')
            )


#--------------------------------------------------------------------------------
#
#  B.3 constructor method...
#
#  Note that the following only applies to the izGrid.extended slot for the object...
#     truncateOverlap = TRUE: assign grid cells within the IZ that are external to the
#                             tract zero values
#                       FALSE: assign grid cells within the IZ, but external to the
#                              tract what would be their estimates if a point fell
#                              fell there--this is useful for illustration and checking
#                              only and must never be used for estimation
#
#  One other thing to be aware of is that if a region folds back outside the original
#  IZ, creating a "phantom" region, the actual grid values (i.e. from using getValues()
#  on the raster grid object) will have associated values for attributes in the
#  data frame, even on some "background" cells. This will only present a potential
#  problem when getValues is used to mask background. It will also only be a problem
#  with some sampling methods. See "The Mirage Method" vignette for details.
#
#
setMethod('izGridMirage',
          signature(izObject = 'InclusionZone', tract='mirageTract'),
function(izObject,
         tract,
         description = 'izGrid object using mirage',
         wholeIZ = TRUE,           #TRUE: grid the whole object; FALSE: just grid the IZ
         truncateOverlap = FALSE,   #by default, see above
         ...
        )
{
#---------------------------------------------------------------------------
#
#   get the perimeter of the IZ object first...
#
    p.iz = perimeter(izObject)
    iz = izObject         
  
#
    require(rgeos)

#
#   create individual SpatialLines objects from each side of the tract...
#
    tr = tract              
    p.tr = perimeter(tract)
    p.trsl = as(p.tr, "SpatialLines")    #cast from SP to SL
    xy = coordinates(p.trsl)[[1]][[1]]   #now the coordinates for each side are rows in a matrix
    west = xy[1:2,]
    West = SpatialLines(list(Lines(Line(west),'west')))
    north = xy[2:3,]
    North = SpatialLines(list(Lines(Line(north),'north')))
    east = xy[3:4,]
    East =  SpatialLines(list(Lines(Line(east),'east')))
    south = xy[4:5,]                                      #tract is a closed polygon, so 5==1
    South = SpatialLines(list(Lines(Line(south),'south')))

    
#
#   see if there are any intersection points of individual bounday legs with IZ--if
#   there are no intersections, call the default constructor as we are done...
#
#   Note that it is not enough to test the iz and tract perimeters for intersection
#   as it will be true for iz completely inside the tract!
#
    gI = rep(FALSE, 4)
    names(gI) = .StemEnv$cardinal
    gI['north'] = gIntersects(North, p.iz)
    gI['south'] = gIntersects(South, p.iz)
    gI['east'] = gIntersects(East, p.iz)
    gI['west'] = gIntersects(West, p.iz)
    #are we are done?...
    if(!any(gI)) {
      griz = izGrid(iz, tract)
      griz = as(griz, 'mirageInclusionZoneGrid')  #cast, all other slots get NULL 
      return(griz)
    }
   
    

#
#   okay, there's at least one intersection...
#
#   first extend the tract so that izGrid will compute what's in the whole iz, including
#   the portion outside the tract...
#
    bbox.iz = bbox(iz)
    bbox.tr = bbox(tr)
    bbs = array(dim=c(2,2,2))         #array of bboxes for overall sum of IZ and tract
    bbs[,,1] = bbox.iz
    bbs[,,2] = bbox.tr
    dimnames(bbs) = dimnames(bbox.tr)
    ext.bb = bboxSum(bbs)             #combine bbox's for the larger tract extents
    ext.tr = extend(tr, ext.bb)       #and extend the tract to a new RasterLayer
    ext.tr = Tract(ext.tr)            #and make it a Tract
    #create the izgrid object on the extended tract...
    izg = izGrid(iz, ext.tr, wholeIZ=wholeIZ, ...)



#
#   The following function is the basic code for miraging each section that falls outside the tract
#   for the cardinal directions...
#

    mirageSection = function(direction = .StemEnv$cardinal) {
      #--------------------------------------------------------------------------------------
      #Note: this function gets any missing info from the enclosure and modifies izg there...
      #--------------------------------------------------------------------------------------
      direction = match.arg(direction)
      #first make an extended bbox rectangle that will encompass the inclusion zone sliver in "direction"...
      extRect = ext.bb                                           #get from enclosure
      if(direction=='north')
        extRect['y','min'] = ymax(tr)
      else if(direction=='south')
        extRect['y','max'] = ymin(tr)
      else if(direction=='west')
        extRect['x','max'] = xmin(tr)
      else if(direction=='east')
        extRect['x','min'] = xmax(tr)
      extRect.poly = bboxToPoly(extRect)                         #make it a SpatialPolygons object
      #now create the sliver of the IZ outside the boundary (within the extended bbox rectangle)...
      sliver.gi = gIntersection(extRect.poly, p.iz)    #returns a SpatialPolygons object of the external IZ sliver
      #check to see if the sliver overlaps any grid cell centers, as crop() will fail if not...
      suppressWarnings({jis = intersect(extent(izg@grid), sliver.gi)})
      if(is.null(jis)) 
        return(list(gi = sliver.gi, grid=NULL))          #sliver too narrow: no grid cell centers inside sliver
      
      sliver.grid = crop(izg@grid, sliver.gi)            #get the grid cells within the polygon as a raster object
      #next transform the xy coordinates of the cells within the sliver: reflect...
      cells = cellsFromExtent(sliver.grid, extent(sliver.grid))  #get the cell indexes from the sliver
      xy = xyFromCell(sliver.grid, cells)                        #and their xy coordinates to match to same ones in izg
      sliver.icx = cellFromXY(izg@grid, xy)                      #the corresponding cell indexes within the larger izg
      #reflect the xy's...
      if(direction=='north') {
        ydiff = xy[,'y'] - ymax(tr)
        xy[,'y'] = ymax(tr) - ydiff
      }
      else if(direction=='south') {
        ydiff = ymin(tr) - xy[,'y']
        xy[,'y'] = ymin(tr) + ydiff
      }
      else if(direction=='west') {
        xdiff = xmin(tr) - xy[,'x']
        xy[,'x'] = xmin(tr) + xdiff
      }
      else if(direction=='east') {
        xdiff = xy[,'x'] - xmax(tr)
        xy[,'x'] = xmax(tr) - xdiff
      }
      #match the transformed xys to the right cell inside the full izg...
      icx = cellFromXY(izg@grid, xy)                         #cell indexes
      #add the contents for all columns in the data frame (except depth)...
      data = izg@data
      attrNames = colnames(data)
      attrNames = attrNames[!attrNames==.StemEnv$ppEstimates$depth]   #don't sum "depth" in reflected areas
      for(attr in attrNames) {
        data[icx, attr] = data[icx, attr] + data[sliver.icx, attr]    #this will work for variable surfaces
        #to truncate the external data, just set to zero for all attributes...
        if(truncateOverlap)                      
          data[sliver.icx, attr] = 0.0           #set to background value
      }
      izg@data <<- data                          #note assignment to parent enclosure

      return(list(gi = sliver.gi, grid=sliver.grid))
    } #mirageSection

#
#   only reflect if there's an intersection: note again that the function adjusts the
#   izg object found here, in its enclosure in each cardinal call below...
#
    if(gI['north']) {
      ms = mirageSection('north')
      north.gi = ms$gi
      north.grid = ms$grid
    }
    else {
      north.gi = NULL
      north.grid = NULL
    }
       
    if(gI['south']) {
      ms = mirageSection('south')
      south.gi = ms$gi
      south.grid = ms$grid
    }
    else {
      south.gi = NULL
      south.grid = NULL
    }
       
    if(gI['east']) {
      ms = mirageSection('east')
      east.gi = ms$gi
      east.grid = ms$grid
    }
    else {
      east.gi = NULL
      east.grid = NULL
    }
       
    if(gI['west']) {
      ms = mirageSection('west')
      west.gi = ms$gi
      west.grid = ms$grid
    }
    else {
      west.gi = NULL
      west.grid = NULL
    }


#
#   truncate/crop the extended izg back to the correct boundary and explicitely crop the
#   data slot since it is not part of the raster image...
#
    grid = crop(izg@grid, tr)                    #original tract interior of izgrid object...
    cells = cellsFromExtent(grid, extent(grid))  #cell indexes from the interior
    xy = xyFromCell(grid, cells)                 #and their xy coordinates
    #match the transformed xys to the right cell inside the full izg...
    icx = cellFromXY(izg@grid, xy)
    #make a new data frame to hold the interior points...
    attrNames = colnames(izg@data)
    attrNames = attrNames[!attrNames==.StemEnv$ppEstimates$depth]
    data = as.data.frame(matrix(NA, nrow=length(cells), ncol=length(attrNames)))
    colnames(data) = attrNames
    for(attr in attrNames) 
      data[cells, attr] = izg@data[icx, attr]

#
#   the depth is one, even in the folded regions, because it denotes the number of trees
#   that lie in the zone, so is one for each tree regardless of the number of mirage
#   reflections that might fold back into certain cells; also note that phantom cells will
#   get a depth of 1 if present...
#
    data$depth = ifelse(data[, attrNames[1]] > 0, 1, 0)
    izg@data$depth = ifelse(izg@data[, attrNames[1]] > 0, 1, 0)


#    
#   and finally build the new mirage izg object...
#
    griz = new('mirageInclusionZoneGrid',
               description = izg@description,
               iz = izg@iz,
               grid = grid,
               data = data,
               bbox = bbox(izg),                       #should still be the same with iz included???? wholeIZ????
               #mirage slots...
               slopOver = gI,
               north.polygon = north.gi,
               north.grid = north.grid,
               south.polygon = south.gi,
               south.grid = south.grid,
               east.polygon = east.gi,
               east.grid = east.grid,
               west.polygon = west.gi,
               west.grid = west.grid,
               izGrid.extended = izg
              )
    
    return(griz)
}   #izGridMirage for'InclusionZone' & 'mirageTract'
)   #setMethod









#================================================================================
#  C.1 plot method for "mirageInclusionZoneGrid" subclass; we could make y a
#      "Tract" subclass object, but then that would require it to be passed, which
#      would be different from the superclass method, so we'll just allow users
#      to pass it as a non-signature argument.
#
setMethod('plot',
          signature(x = 'mirageInclusionZoneGrid', y='missing'),
function(x,
         tract = NULL,
         tractBorder = 'red',
         lwdTract = 1.25,
         showExtended = FALSE,
         showReflectedSlivers = FALSE,
         colReflections = 'blue',
         ...
        )
{
#------------------------------------------------------------------------------
#
#  Arguments...
#    x = a "mirageInclusionZoneGrid" object
#    tract = a "Tract" object or NULL
#    tractBorder = color for the tract border polygon
#    lwdTract = the lwd for the tract border
#    showExtended = TRUE: display the extended mIZG object with external grid
#                         cells w/in the IZ have the estimates;
#                   FALSE: display only the internal (actual) surface
#    showReflectedSlivers = TRUE: display the external sliver polygons for the
#                                 IZs reflected back into the tract if any;
#                           FALSE: don't display them
#                           Note: only applicable when !is.null(tract)
#    colReflections = color for the reflected polygons
#    ... = other arguments passed on to the base "InclusionZoneGrid" plot method
#
#  Returns: invisibly
#
#------------------------------------------------------------------------------
#
#  By default, we use the base "InclusionZoneGrid" plot method to start -- no
#  tract object needed for this...
#
    mizg = x
    if(showExtended && is.null(mizg@izGrid.extended))
      showExtended = FALSE
    if(showExtended)
      #plot(mizg@izGrid.extended, ...)
      callNextMethod(mizg@izGrid.extended, ...)
    else
      #plot(mizg, ...)
      callNextMethod(mizg, ...)

#
#   add the tract perimeter if desired...
#
    if(!is.null(tract)) {
      if(!is(tract, 'Tract'))
        stop('You must pass an object of class "Tract" or subclass in \"Tract\"!')
      p.tr = perimeter(tract)
      suppressWarnings({                            #for 'add'
        plot(p.tr, add=TRUE, border=tractBorder, lwd=lwdTract)
      })

#
#     we need the tract to do slivers...
#
      if(showReflectedSlivers) {
        if(mizg@slopOver['north']) {   #add north reflection...
          xy1.rn = coordinates(mizg@north.polygon@polygons[[1]]@Polygons[[1]])
          ydiff1 = xy1.rn[,2] - ymax(tract)
          xy1.rn[,2] = ymax(tract) - ydiff1
          lines(xy1.rn, col=colReflections, lty='dashed')
        }
        if(mizg@slopOver['east']) {   #add east reflection...
          xy1.re = coordinates(mizg@east.polygon@polygons[[1]]@Polygons[[1]])
          xdiff1 = xy1.re[,1] - xmax(tract)
          xy1.re[,1] = xmax(tract) - xdiff1
          lines(xy1.re, col=colReflections, lty='dashed')
        } 
        if(mizg@slopOver['south']) {   #add south reflection...
          xy1.rs = coordinates(mizg@south.polygon@polygons[[1]]@Polygons[[1]])
          ydiff1 = ymin(tract) - xy1.rs[,2]
          xy1.rs[,2] = ydiff1 + ymin(tract)
          lines(xy1.rs, col=colReflections, lty='dashed')
        }
        if(mizg@slopOver['west']) {   #add east reflection...
          xy1.rw = coordinates(mizg@west.polygon@polygons[[1]]@Polygons[[1]])
          xdiff1 = xmin(tract) - xy1.rw[,1]
          xy1.rw[,1] = xdiff1 + xmin(tract)
          lines(xy1.rw, col=colReflections, lty='dashed')
        } 
      } #reflectedSlivers
    } #tract available
    
    return(invisible())

}    #plot for 'mirageInclusionZoneGrid' base
) #setMethod
