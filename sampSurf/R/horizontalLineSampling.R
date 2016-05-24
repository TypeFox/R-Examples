#---------------------------------------------------------------------------
#
#   This file has the collection of routines necessary to implement horizontal
#   line sampling for standing trees. The routines include...
#
#   1. 'horizontalLineIZ' class structure
#   2. 'horizontalLineIZ' constructor function
#   3. plot method for 'horizontalLineIZ' objects
#   4. summary method for 'horizontalLineIZ' objects
#   5. 'izGrid' method for 'horizontalLineIZ' objects
#   6. perimeter method for 'horizontalLineIZ' objects
#
#   Please note that the "lineSegment" class and associated methods are in
#   the respective collective files (e.g. ArealSamplingClass.R) as it is
#   used in other line-oriented sampling methods.
#
#Author...									Date: 9-Oct-2012
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#



#=================================================================================================
#
# 1. the horizontalLineIZ class is a direct descendant of 'standingTreeIZ'...
#
#
setClass('horizontalLineIZ',
    representation(angleGauge = 'angleGauge',       #angleGauge object
                   lineSegment = 'lineSegment',     #lineSegment object
                   width = 'numeric',               #iz width == 2*R
                   length = 'numeric',              #iz length
                   izPerim = 'matrix',              #holds the iz perimeter in matrix form w/ hc
                   area = 'numeric',                #exact area of the inclusion zone
                   perimeter = 'SpatialPolygons',   #izone perimeter in 'SpatialPolygons' form
                   pgArea = 'numeric'               #polygon izone area approximation                   
                  ),
    contains = 'standingTreeIZ',                    #a subclass of the 'standingTreeIZ' class
    validity = function(object) {

                 if(!identical(object@standingTree@units, object@angleGauge@units))
                   return('measurement units do not match between standingTree and angleGauge objects!')

                 if(!identical(object@standingTree@units, object@lineSegment@units))
                   return('measurement units do not match between standingTree and lineSegment objects!')

                 if(!isTRUE(all.equal(object@area, object@length*object@width))) #close is good enough
                   return('incorrect area for horizontalLineIZ object!')

                 if(!isTRUE(all.equal(object@width, object@standingTree@dbh*object@angleGauge@PRF*2)))
                   return('horizontalLineIZ width not correct according to encapsulated objects!')
                 
                 return(TRUE)
               } #validity check
) #class horizontalLineIZ 









#================================================================================
# 2. generic definition...
#
if(!isGeneric("horizontalLineIZ")) 
  setGeneric('horizontalLineIZ',  
             function(standingTree, angleGauge, lineLength, ...) standardGeneric('horizontalLineIZ'),
             signature = c('standingTree', 'angleGauge', 'lineLength')
            )
 

       
#================================================================================
#  and associated constructor method for class horizontalLineIZ...
#
setMethod('horizontalLineIZ',
          signature(standingTree = 'standingTree', angleGauge = 'angleGauge', lineLength = 'numeric'),
function(standingTree,
         angleGauge,
         lineLength,                           #in ft or m
         orientation = 0,                      #in degrees from North (positive y-axis)
         description = 'inclusion zone for horizontal line sampling method',
         spID = paste('hls',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         ...
        )
{
#------------------------------------------------------------------------------    
#
#   restrict the line segment to be at least the diameter of the tree...
#
    dbh = standingTree@dbh
    if(lineLength <= dbh)
      stop(paste('The line length (',lineLength,'m) must be larger than tree dbh (',
                 dbh,'m)!', sep=''))

    units = standingTree@units
    loc = coordinates(standingTree@location)[1,]  #need a vector from matrix    
    lineSeg = lineSegment(lineLength, orientation, units=units, centerPoint=loc) #this will check for valid arguments

#
#   per unit area estimates...
#
    df = angleGauge@DF                                   #diameter factor yielding ft/ft or m/m
    if(standingTree@units==.StemEnv$msrUnits$metric) {
      baFactor = .StemEnv$baFactor['metric']
      dbhFactor = df/lineLength                          #remember, dbh is stored in m or ft, not cm or inches
      conv = .StemEnv$m2cm
    }
    else {
      baFactor = .StemEnv$baFactor['English']
      dbhFactor = df/lineLength                          #remember, dbh is stored in m or ft, not cm or inches
      conv = .StemEnv$ft2in
    }
    ba = baFactor*dbh*dbh  
    puaBlowup = dbhFactor/dbh 
    puaEstimates = list(standingTree@treeVol*puaBlowup,
                        puaBlowup,
                        ba*puaBlowup, 
                        standingTree@surfaceArea*puaBlowup, 
                        standingTree@biomass*puaBlowup,
                        standingTree@carbon*puaBlowup
                       )
    names(puaEstimates) = .StemEnv$puaEstimates[c('volume', 'Density', 'basalArea',
                                                  'surfaceArea', 'biomass', 'carbon'
                                                )]

    
#
#   create the inclusion zone horizontally, then rotate it...
#
    k = angleGauge@k
    R = dbh/k                                           #remember, dbh is stored in m or ft, not cm or inches
    ls = lineSegment(lineLength, 0, units=units, centerPoint=c(x=0,y=0))  #orient due north
    pts = coordinates(ls@segment)$LinesSeg$LineSeg      #phew...
    izPerim = pts
    izPerim = rbind(izPerim, izPerim[2,], izPerim[1,])
    izPerim[,'x'] = c(R, R, -R, -R)
    izPerim = cbind(izPerim, rep(1,4))                  #homogeneous
    izPerim = rbind(izPerim, izPerim[1,])               #close the polygon

    orientation = .StemEnv$deg2Rad(-orientation)        #lineSegment() handles this already
    trMat = transfMatrix(orientation, loc)    
    izPerim = izPerim %*% trMat
    dimnames(izPerim) = list(NULL,c('x','y','hc'))
    
#
#   and make a SpatialPolygons object...
#
    pg = Polygon(izPerim[,-3])                             #sans hc
    pgs = Polygons(list(pg=pg), ID=spID)
    spObj = SpatialPolygons(list(pgs=pgs), proj4string = spUnits)
    pgArea = pg@area
    

#
#   create the object...
#
    width = 2*R
    hlIZ = new('horizontalLineIZ',
               standingTree = standingTree,
               lineSegment = lineSeg,
               angleGauge = angleGauge,
               width = width,
               length = lineLength,
               izPerim = izPerim,                  #matrix representation of perimeter
               perimeter = spObj,                  #SpatialPolygons perimeter
               bbox = bbox(spObj),
               spUnits = spUnits,
               description = description,
               units = units,
               puaBlowup = puaBlowup,
               puaEstimates = puaEstimates,
               area = lineLength * width,
               pgArea = pgArea
              )

    return(hlIZ)
}   #horizontalLineIZ constructor
)   #setMethod







#================================================================================
#  3. plot method for horizontalLineIZ subclass...
#
setMethod('plot',
          signature(x = 'horizontalLineIZ', y='missing'),
function(x, 
         axes = FALSE,                     #not a par() so can't be passed to callNextMethod, so separate it
         showLineSegment = TRUE,
         ltyLineSegment = 'dashed',
         showTree = TRUE,
         izColor = .StemEnv$izColor,
         izBorder = .StemEnv$izBorderColor,
         add = FALSE,                           #add each IZ to overall plot if TRUE
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
#
#   plots the horizontalLineIZ object...
#
#------------------------------------------------------------------------------
#
    object = x

    if(!add)
      callNextMethod(object, axes=axes, asp=asp, ...)            #setup extents

    suppressWarnings({                                #for object-specific parameters not in par() ...
      plot(object@perimeter, col=izColor, border=izBorder, axes=axes, add=TRUE, ...)
      if(showLineSegment)
        plot(object@lineSegment, add=TRUE, lty=ltyLineSegment, ...)
      if(showTree)
        plot(object@standingTree, add=TRUE, ...)
    })
                     
    return(invisible())
}   #plot for 'horizontalLineIZ'
)   #setMethod






#================================================================================
#  4. summary method for class "horizontalLineIZ"...
#
setMethod('summary',
          signature(object = 'horizontalLineIZ'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to 'InclusionZone' & 'standingTreeIZ' methods for 'horizontalLineIZ'...
#------------------------------------------------------------------------------
    callNextMethod()
    if(object@units == .StemEnv$msrUnits$metric) {
      su = ' square meters'
      lu = ' meters'
    }
    else {
      su = ' square feet'
      lu = ' feet'
    }
    
    cat('\nhorizontalLineIZ...')
    cat('\n  Inclusion zone length = ', object@length, lu, sep='')
    cat('\n  Inclusion zone width = ', object@width, lu, sep='')
    cat('\n  Inclusion zone area = ', object@area, su, sep='')
    cat('\n  Note: use \"summary\" on the angleGauge & lineSegment slots for details')

    cat('\n')
    
    return(invisible())
}   #summary for 'horizontalLineIZ'
) #setMethod





#================================================================================
#  5. izGrid method for 'horizontalLineIZ' and 'Tract' classes...
#
setMethod('izGrid',
          signature(izObject = 'horizontalLineIZ', tract='Tract'),
function(izObject,
         tract,
         description = 'horizontalLineIZ inclusion zone grid object',
         wholeIZ = TRUE,           #TRUE: grid the whole object; FALSE: just grid the IZ
         ...
        )
{
#---------------------------------------------------------------------------
#
#
    griz = izGridConstruct(izObject=izObject, tract=tract, description=description,
                           wholeIZ=wholeIZ, ...)
    return(griz)
}   #izGrid for'horizontalLineIZ'
)   #setMethod




#================================================================================
# 6. perimeter method for horizontalLineIZ object...
#
setMethod('perimeter',
          signature(object = 'horizontalLineIZ'),
function(object, ...)
{
    return(object@perimeter)
}   #horizontalLineIZ
) #setMethod
