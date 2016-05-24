#---------------------------------------------------------------------------
#
#   This file holds the S4 definition for the constructor methods of the
#   InclusionZone class & subclasses...
#
#   The methods include...
#
#        ...downLogIZ subclass constructors...
#
#     1. a constructor for 'standUpIZ'     (Aug 2010)
#     2. for 'chainSawIZ'                  (Fall 2010)
#     3. for 'sausageIZ'                   (Fall 2010)
#     4. for 'fullChainSawIZ'              (Sept 2013)
#     5. for 'pointRelascopeIZ'            (Jan 2011)
#     6. for 'perpendicularDistanceIZ'     (Jan 2011)
#     7. for 'omnibusPDSIZ'                (Feb 2011)
#     8. for 'distanceLimitedPDSIZ'        (Feb-Mar 2011)
#     9. for 'omnibusDLPDSIZ'              (Mar 2011)
#    10. for 'hybridDLPDSIZ'               (July 2011)
#    11. for 'distanceLimitedIZ'           (Mar 2011)
#    12. for 'distanceLimitedMCIZ'         (May 2011)
#
#        ...standingTreeIZ subclass constructors...
#
#     1. for 'circularPlotIZ'              (Dec 2011)
#     2. for 'horizontalPointIZ'           (Dec 2011)
#
#   Note: I adopted a new strategy after HPS was added, new sampling methods
#         have all their respective code for all classes in one file. 
#         Please see those files for their InclusionZone constructors.
#
#   Note that the sp, raster, and rgeos packages must be installed.
#   The gpclib package was originally used for the chainSawIZ method,
#   but it has a restrictive license; rgeos turns out to be equivalent
#   and performs just as well for what is required here.
#
#   In each case, all estimated quatities are calculated for the object;
#   i.e., volume & Density; whatever else might be desired should
#   be added later.
#
#Author...									Date: 23-Aug-2010
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
#if(!isGeneric("standUpIZ")) 
  setGeneric('standUpIZ',  
             function(downLog, plotRadius, ...) standardGeneric('standUpIZ'),
             signature = c('downLog', 'plotRadius')
            )

#if(!isGeneric("chainSawIZ")) 
  setGeneric('chainSawIZ',  
             function(downLog, plotRadius, ...) standardGeneric('chainSawIZ'),
             signature = c('downLog', 'plotRadius')
            )

#if(!isGeneric("sausageIZ")) 
  setGeneric('sausageIZ',  
             function(downLog, plotRadius, ...) standardGeneric('sausageIZ'),
             signature = c('downLog', 'plotRadius')
            )

#if(!isGeneric("fullChainSawIZ")) 
  setGeneric('fullChainSawIZ',  
             function(downLog, plotRadius, ...) standardGeneric('fullChainSawIZ'),
             signature = c('downLog', 'plotRadius')
            )

#if(!isGeneric("pointRelascopeIZ")) 
  setGeneric('pointRelascopeIZ',  
             function(downLog, prs, ...) standardGeneric('pointRelascopeIZ'),
             signature = c('downLog', 'prs')
            )

#if(!isGeneric("perpendicularDistanceIZ")) 
  setGeneric('perpendicularDistanceIZ',  
             function(downLog, pds, ...) standardGeneric('perpendicularDistanceIZ'),
             signature = c('downLog', 'pds')
            )

#if(!isGeneric("omnibusPDSIZ")) 
  setGeneric('omnibusPDSIZ',  
             function(downLog, pds, ...) standardGeneric('omnibusPDSIZ'),
             signature = c('downLog', 'pds')
            )

#if(!isGeneric("distanceLimitedPDSIZ")) 
  setGeneric('distanceLimitedPDSIZ',  
             function(downLog, pds, dls, ...) standardGeneric('distanceLimitedPDSIZ'),
             signature = c('downLog', 'pds', 'dls')
            )

#if(!isGeneric("omnibusDLPDSIZ")) 
  setGeneric('omnibusDLPDSIZ',  
             function(downLog, pds, dls, ...) standardGeneric('omnibusDLPDSIZ'),
             signature = c('downLog', 'pds', 'dls')
            )

#if(!isGeneric("hybridDLPDSIZ")) 
  setGeneric('hybridDLPDSIZ',  
             function(downLog, pds, dls, ...) standardGeneric('hybridDLPDSIZ'),
             signature = c('downLog', 'pds', 'dls')
            )

#if(!isGeneric("distanceLimitedMCIZ")) 
  setGeneric('distanceLimitedMCIZ',  
             function(downLog, dls, ...) standardGeneric('distanceLimitedMCIZ'),
             signature = c('downLog', 'dls')
            )

#if(!isGeneric("distanceLimitedIZ")) 
  setGeneric('distanceLimitedIZ',  
             function(downLog, dls, ...) standardGeneric('distanceLimitedIZ'),
             signature = c('downLog', 'dls')
            )
   



       
#================================================================================
#  1. method for functions and class standUpIZ...
#
setMethod('standUpIZ',
          signature(downLog = 'downLog', plotRadius = 'numeric'),
function(downLog,
         plotRadius,
         description = 'inclusion zone for "standup" method',
         spID = paste('su',.StemEnv$randomID(),sep=':'),
         #spID = unlist(strsplit(tempfile('su:',''),'\\/'))[2],
         spUnits = CRS(projargs=as.character(NA)),
         ...
        )
{
#------------------------------------------------------------------------------
#
#   get bbox from the downLog object...
#
    downLog.bbox = bbox(downLog@spLog)
   
#
#   now also from the circularPlot object...
#   --we must be careful below and not just pass the "..." because we are
#     explicitly passing "units" and if the user also put "units" in the "..."
#     call to the current method, then it will be passed twice throwing an
#     error; so pull out any other arguments from "..." that could go to
#     circularPlot and pass them explicitly or set them to default; i.e.,
#     nptsPerimeter
#
    units = downLog@units
    loc = coordinates(downLog@slNeedleAxis)$naLines$naLine[1,]    #transformed butt center point
    args = list(...)
    if('nptsPerimeter' %in% names(args))          #user can override here
      nptsPerimeter = args$nptsPerimeter          #user's choice
    else
      nptsPerimeter = 100                         #default in circularPlot()
    circularPlot = circularPlot(plotRadius, units=units, centerPoint=loc,
                                spID = spID, spUnits = spUnits,
                                nptsPerimeter = nptsPerimeter)
    cp.bbox = bbox(circularPlot@perimeter)

#
#   combine them for the overall bbox...
#
    min = apply(cbind(downLog.bbox, cp.bbox), 1, min)
    max = apply(cbind(downLog.bbox, cp.bbox), 1, max)
    bbox = matrix(cbind(min,max),ncol=2, dimnames=list(c('x','y'), c('min','max')))

#
#   per unit area estimates...
#
    unitArea = ifelse(downLog@units==.StemEnv$msrUnits$English, .StemEnv$sfpAcre, .StemEnv$smpHectare) 
    puaBlowup = unitArea/circularPlot@area 
    puaEstimates = list(downLog@logVol*puaBlowup, puaBlowup, downLog@logLen*puaBlowup,
                        downLog@surfaceArea*puaBlowup, downLog@coverageArea*puaBlowup,
                        downLog@biomass*puaBlowup, downLog@carbon*puaBlowup
                       )
    names(puaEstimates) = .StemEnv$puaEstimates[c('volume', 'Density', 'Length',
                                                  'surfaceArea', 'coverageArea',
                                                  'biomass', 'carbon'
                                                )]

#
#   create the object...
#
    dlIZ = new('standUpIZ', downLog=downLog, circularPlot=circularPlot,
               bbox=bbox, spUnits=spUnits, description=description,
               units = units, puaBlowup = puaBlowup, puaEstimates = puaEstimates
              )

    return(dlIZ)
}   #standUpIZ constructor
)   #setMethod


   



       
#================================================================================
#  2. method for functions and class chainSawIZ...
#
setMethod('chainSawIZ',
          signature(downLog = 'downLog', plotRadius = 'numeric'),
function(downLog,
         plotRadius,
         plotCenter = c(x=0, y=0),
         description = 'inclusion zone for "chainsaw" method',
         spID = paste('cs',.StemEnv$randomID(),sep=':'),
         #spID = unlist(strsplit(tempfile('cs:',''),'\\/'))[2],
         #spID = paste('cs',format(runif(1,0,10000),digits=8),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         ...
        )
{
#------------------------------------------------------------------------------
#
#   Currently, if there is no intersection between the plot and the log, the
#   routine will stop with an error message. This may not be the best thing
#   in general, but in sampling surface simulations, only plots whose centers
#   fall within the sausage-shaped inclusion zone should be passed here, so
#   unless someone tries to use this for some other purpose, it will not matter.
#------------------------------------------------------------------------------
#
#   rgeos is suggested, make sure the user has it available...
#
    if(!require(rgeos))
      stop('The \"rgeos\" package must be installed to use the chainSaw method!')

  
#
#   get bbox from the downLog object...
#
    downLog.bbox = bbox(downLog@spLog)
    logID = downLog@spLog@polygons$pgsLog@ID
   
#
#   now also from the circularPlot object...
#   --we must be careful below and not just pass the "..." because we are
#     explicitly passing "units" and if the user also put "units" in the "..."
#     call to the current method, then it will be passed twice throwing an
#     error; so pull out any other arguments from "..." that could go to
#     circularPlot and pass them explicitly or set them to default; i.e.,
#     nptsPerimeter
#
    units = downLog@units
    args = list(...)
    if('nptsPerimeter' %in% names(args))          #user can override here
      nptsPerimeter = args$nptsPerimeter          #user's choice
    else
      nptsPerimeter = 100                         #default in circularPlot()
    circularPlot = circularPlot(plotRadius, units=units, centerPoint=plotCenter,
                                spID = spID, spUnits = spUnits,
                                nptsPerimeter = nptsPerimeter)
    cp.bbox = bbox(circularPlot@perimeter)
    cpCoords = coordinates(circularPlot@perimeter@polygons$pgsCircPlot@Polygons$circPlot)


#
#   combine them for the overall bbox...
#
    min = apply(cbind(downLog.bbox, cp.bbox), 1, min)
    max = apply(cbind(downLog.bbox, cp.bbox), 1, max)
    bbox = matrix(cbind(min,max),ncol=2, dimnames=list(c('x','y'), c('min','max')))

#
#   per unit area estimates...
#
    unitArea = ifelse(downLog@units==.StemEnv$msrUnits$English, .StemEnv$sfpAcre, .StemEnv$smpHectare) 
    puaBlowup = unitArea/circularPlot@area 

#
#   this uses the rgeos (or gpclib) routines for intersection of plot & log polygons...
#
    gPlot = as(cpCoords, 'gpc.poly')                     #make a gpc polygon from the plot
    gLog = as(coordinates(downLog@spLog@polygons$pgsLog@Polygons$log), 'gpc.poly') #gpc polygon of the log 
    ov.log = intersect(gPlot, gLog)                      #overlay, result is of class gpc.poly
    if(length(ov.log@pts) == 0)                          #check for no intersection
      stop(paste('**>No plot-log intersection at plotCenter=',plotCenter[1],',',plotCenter[2]))
    mSect = as(ov.log, 'matrix')                         #convert sliver section to matrix for sp
    mSect = rbind(mSect, mSect[1,])                      #close the polygon
    pgSect = Polygon(mSect)                              #sliver matrix to sp classes...
    pgsSect = Polygons(list(sliver=pgSect), ID=paste(logID,'Sliver',sep=':')) 
    spSect = SpatialPolygons(list(pgsSliver=pgsSect))       #takes a list of Polygons objects

    css = chainsawSliver(downLog, sect = mSect, gLog = gLog, ...)

    puaEstimates = list(css$sectVol*puaBlowup, puaBlowup, css$boltLen*puaBlowup,
                        css$sectSA*puaBlowup, css$sectCA*puaBlowup,
                        css$sectBms*puaBlowup, css$sectCarbon*puaBlowup
                       )
    names(puaEstimates) = .StemEnv$puaEstimates[c('volume', 'Density', 'Length',
                                                  'surfaceArea', 'coverageArea',
                                                  'biomass', 'carbon'
                                                )]

#
#   and make a SpatialPolygons object out of the bolt...
#
    pgBolt = Polygon(css$rotBolt[,-3])                             #sans hc
    pgsBolt = Polygons(list(sausage=pgBolt), ID=paste(logID,'BBolt',sep=':'))
    spBolt = SpatialPolygons(list(pgsBolt=pgsBolt),      #takes a list of Polygons objects
                                proj4string = spUnits                       
                               )
    css$spBolt = spBolt
    
#
#   create the object...
#
    csIZ = new('chainSawIZ', downLog=downLog, circularPlot=circularPlot,
               bbox=bbox, spUnits=spUnits, description=description,
               units = units, puaBlowup = puaBlowup, puaEstimates = puaEstimates,
               sliver = spSect,
               bolt = css
              )

    return(csIZ)
        

}   #chainSawIZ constructor
)   #setMethod





    


       
#================================================================================
#  3. method for functions and class sausageIZ...
#
setMethod('sausageIZ',
          signature(downLog = 'downLog', plotRadius = 'numeric'),
function(downLog,
         plotRadius,
         nptsHalfCircle = 50,          #number of points in each end's half circle
         description = 'inclusion zone for downed log "sausage" sampling method',
         spID = paste('saus',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         ...
        )
{
#------------------------------------------------------------------------------
#
#   transformation matrix...
#
    centerOffset = coordinates(downLog@location)
    logAngle = downLog@logAngle
    trMat = transfMatrix(logAngle, centerOffset)

    
#
#   direct calculation of the inclusion zone area and blowup factor...
#
    logLen = downLog@logLen
    izArea = 2.0*logLen*plotRadius + pi*plotRadius*plotRadius 

#
#   per unit area estimates...
#
    unitArea = ifelse(downLog@units==.StemEnv$msrUnits$English, .StemEnv$sfpAcre, .StemEnv$smpHectare) 
    puaBlowup = unitArea/izArea 
    puaEstimates = list(downLog@logVol*puaBlowup, puaBlowup, downLog@logLen*puaBlowup,
                        downLog@surfaceArea*puaBlowup, downLog@coverageArea*puaBlowup,
                        downLog@biomass*puaBlowup, downLog@carbon*puaBlowup
                       )
    names(puaEstimates) = .StemEnv$puaEstimates[c('volume', 'Density', 'Length',
                                                  'surfaceArea', 'coverageArea',
                                                  'biomass', 'carbon'
                                                )]
  
#
#   make sure things are even...
#
    halfLen = logLen/2
    if(nptsHalfCircle > 10)
      npts = nptsHalfCircle
    else
      npts = 50
    
    if(npts%%2 > 0)             #augment if odd
      npts = npts+1
    nptsHalfCircle = npts/2

#
#   left half of the circle, then right...
#
    circ.left = rev( seq(pi/2, 3*pi/2, len=npts) )
    circ.right = rev( c(seq(3*pi/2, 2*pi, len=nptsHalfCircle), seq(0, pi/2, len=nptsHalfCircle)) )
    
#
#   make some sausage...
#
    sausage = matrix(c(-halfLen, plotRadius, 1), nrow=1)                      #left point on top "side"
    
#   right-half circle starts at top-right point, ends at bottom-right point...    
    sausage = rbind(sausage, cbind(halfLen + plotRadius*cos(circ.right),   
                                   plotRadius*sin(circ.right), rep(1,npts)) )
    
#   left half of circle starts at bottom-left, ends at top-left (beginning) point...
    sausage = rbind(sausage, cbind(-halfLen + plotRadius*cos(circ.left),     
                                   plotRadius*sin(circ.left), rep(1,npts)) )

#    
#   any little difference between start & end pts with identical() can mess up the
#   the sp package Polygon routine, so set the end point exactly to start, then transform...
#
    sausage[2*npts+1,] = sausage[1,]
    
    sausage = sausage %*% trMat
    dimnames(sausage) = list(NULL,c('x','y','hc'))

#
#   and make a SpatialPolygons object...
#
    pgSausage = Polygon(sausage[,-3])                             #sans hc
    pgsSausage = Polygons(list(sausage=pgSausage), ID=spID)
    spSausage = SpatialPolygons(list(pgsSausage=pgsSausage),      #takes a list of Polygons objects
                                proj4string = spUnits                       
                               )
    pgSausageArea = pgSausage@area

#
#   create the object...
#
    sausIZ = new('sausageIZ', downLog=downLog,
                 sausage = sausage,                  #matrix representation of perimeter
                 perimeter = spSausage,              #SpatialPolygons perimeter
                 pgSausageArea = pgSausageArea,      #area of SpatialPolygons sausage: approximate
                 spUnits = spUnits,                  #CRS units
                 description = description,          #a comment
                 units = downLog@units,              #units of measure
                 area = izArea,                      #exact inclusion zone area
                 puaBlowup = puaBlowup,              #sausage per unit area blowup factor
                 puaEstimates = puaEstimates,        #per unit area estimates
                 radius = plotRadius,                #plot radius for sausage
                 bbox = bbox(spSausage)              #overall bounding box--redundant here
                )

    return(sausIZ)
}   #sausageIZ constructor
)   #setMethod
    





    


       
#================================================================================
#  4. method for functions and class fullChainSawIZ -- it just calls its parent
#     for a suasageIZ...
#
setMethod('fullChainSawIZ',
          signature(downLog = 'downLog', plotRadius = 'numeric'),
function(downLog,
         plotRadius,
         nptsHalfCircle = 50,          #number of points in each end's half circle
         description = 'Inclusion zone for "fullChainSaw" sampling method',
         spID = paste('fcs', .StemEnv$randomID(), sep=':'),
         spUnits = CRS(projargs = as.character(NA)),
         ...
        )
{
#------------------------------------------------------------------------------
#
    sausIZ = sausageIZ(downLog = downLog,
                       plotRadius = plotRadius,
                       nptsHalfCircle = nptsHalfCircle,
                       description = description,
                       spID = spID,
                       spUnits = spUnits,
                       ...)

    fcsIZ = as(sausIZ, 'fullChainSawIZ')

#
#   volume depends on position, so assign NA here, same for carbon & biomass, also Length,
#   surface and coverage are bolt-based, these are computed in izGrid; Density is okay,
#   see chainSawSliver for details...
#
    fcsIZ@puaEstimates[c('volume', 'Length', 'surfaceArea', 'coverageArea', 'biomass', 'carbon')] = NA_real_

    return(fcsIZ)
}   #fullChainSawIZ constructor
)   #setMethod
                   


    


       
#================================================================================
#  5. method for functions and class pointRelascopeIZ...
#
setMethod('pointRelascopeIZ',
          signature(downLog = 'downLog', prs = 'pointRelascope'),
function(downLog,
         prs,
         nptsCircle = 100,          #number of points in each dual circle
         description = 'inclusion zone for down log point relascope sampling method',
         spID = paste('prs',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         ...
        )
{
#------------------------------------------------------------------------------
#
#   Note that we make the blob as if it were centered at (0,0), with the top
#   half in positive y, and bottom in negative y--this allows us to then translate
#   and rotate to the same specs as the downLog quite easily
#
#   transformation matrix...
#
    centerOffset = coordinates(downLog@location)
    logAngle = downLog@logAngle
    trMat = transfMatrix(logAngle, centerOffset)
    
#
#   direct calculation of the inclusion zone area and blowup factor...
#
    logLen = downLog@logLen
    izArea = prs@phi*logLen*logLen 
    halfLen = logLen/2

#
#   per unit area estimates...
#
    unitArea = ifelse(downLog@units==.StemEnv$msrUnits$English, .StemEnv$sfpAcre, .StemEnv$smpHectare) 
    puaBlowup = unitArea/izArea 
    puaEstimates = list(downLog@logVol*puaBlowup, puaBlowup, downLog@logLen*puaBlowup,
                        downLog@surfaceArea*puaBlowup, downLog@coverageArea*puaBlowup,
                        downLog@biomass*puaBlowup, downLog@carbon*puaBlowup
                       )
    names(puaEstimates) = .StemEnv$puaEstimates[c('volume', 'Density', 'Length',
                                                  'surfaceArea', 'coverageArea',
                                                  'biomass', 'carbon'
                                                )]
  
#
#   make sure things are even...
#
    if(nptsCircle > 10)
      npts = nptsCircle
    else
      npts = 100
    
    if(npts%%2 > 0)             #augment if odd
      npts = npts+1

#
#   common radius,
#
    nu = prs@angleRadians
    R = halfLen/sin(nu)

#
#   dual circle centers...
#
    a = R*cos(nu)
    dualCenters = matrix(c(0,a,1, 0,-a,1), nrow=2, byrow=TRUE)
    dualCenters = (dualCenters %*% trMat)[,-3]
    dimnames(dualCenters) = list(NULL, c('x','y'))
    
#
#   get the circular points for the top half of the blob...
#
    halfBlob = seq(0-(pi/2-nu), 3*pi/2-nu, len=npts)     #begin/end points are easy to show
    blob = matrix(NA, nrow=2*npts+1, ncol=3)             #matrix of the blob/dual circle coordinates
    
#   top part of blob...    
    x = R*cos(halfBlob)
    y = a + R*sin(halfBlob)
    hc = rep(1,2*npts)
#   duplicate it for the bottom half...
    xx = c(x, rev(x))
    yy = c(y, rev(-y))
    blob[1:(2*npts),] = cbind(xx, yy, hc)


#    
#   any little difference between start & end pts with identical() can mess up the
#   the sp package Polygon routine, so set the end point exactly to start, then transform...
#
    blob[2*npts+1,] = blob[1,]
    
    blob = blob %*% trMat
    dimnames(blob) = list(NULL,c('x','y','hc'))


#
#   and make a SpatialPolygons object...
#
    pgPRS = Polygon(blob[,-3])                             #sans hc
    pgsPRS = Polygons(list(pgPRS=pgPRS), ID=spID)
    spPRS = SpatialPolygons(list(pgsPRS=pgsPRS),        #takes a list of Polygons objects
                             proj4string = spUnits                       
                            )
    pgPRSArea = pgPRS@area

#
#   create the object...
#
    prsIZ = new('pointRelascopeIZ', downLog=downLog,
                 prs = prs,                          #point realscope sampling object
                 dualCircle = blob,                  #matrix representation of perimeter
                 perimeter = spPRS,                  #SpatialPolygons perimeter
                 pgDualArea = pgPRSArea,             #area of SpatialPolygons blob: approximate
                 dualCenters = dualCenters,          #centers of each dual circle
                 spUnits = spUnits,                  #CRS units
                 description = description,          #a comment
                 units = downLog@units,              #units of measure
                 area = izArea,                      #exact inclusion zone area
                 puaBlowup = puaBlowup,              #sausage per unit area blowup factor
                 puaEstimates = puaEstimates,        #per unit area estimates
                 radius = R,                         #plot radius for each half blob
                 bbox = bbox(spPRS)                  #overall bounding box--redundant here
                )

    return(prsIZ)
}   #pointRelascopeIZ constructor
)   #setMethod

    




    


       
#================================================================================
#  6. method for functions and class perpendicularDistanceIZ...
#
setMethod('perpendicularDistanceIZ',
          signature(downLog = 'downLog', pds = 'perpendicularDistance'),
function(downLog,
         pds,
         description = 'inclusion zone for down log perpendicular distance sampling',
         spID = paste('pds',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         pdsType = .StemEnv$pdsTypes,
         ...
        )
{
#------------------------------------------------------------------------------
#
#   a quick check...
#
    if(downLog@units != pds@units)
      stop('units are not the same for downLog and pds!')
  
#
#   which kind of object do we build? And what PPS type in pds?
#
    pdsType = match.arg(pdsType)
  
#
#   put in a warning for number of points in taper data frame too few <20 for now<<<<<******
#
#   transformation matrix...
#
    centerOffset = coordinates(downLog@location)
    logAngle = downLog@logAngle
    trMat = transfMatrix(logAngle, centerOffset)
    trMatInv = solve(trMat)
    
#
#   direct calculation of the inclusion zone area and blowup factor...
#
    izArea = switch(pdsType,
                    volume = 2*pds@kpds*downLog@logVol,
                    surfaceArea = 2*pds@kpds*downLog@surfaceArea,
                    coverageArea = 2*pds@kpds*downLog@coverageArea,
                    stop('Illegal pdsType in perpendicularDistanceIZ!')
                   )    

#
#   per unit area estimates...
#
    unitArea = ifelse(downLog@units==.StemEnv$msrUnits$English, .StemEnv$sfpAcre, .StemEnv$smpHectare) 
    puaBlowup = unitArea/izArea 
    puaEstimates = list(downLog@logVol*puaBlowup, puaBlowup, downLog@logLen*puaBlowup,
                        downLog@surfaceArea*puaBlowup, downLog@coverageArea*puaBlowup,
                        downLog@biomass*puaBlowup, downLog@carbon*puaBlowup
                       )
    names(puaEstimates) = .StemEnv$puaEstimates[c('volume', 'Density', 'Length',
                                                  'surfaceArea', 'coverageArea',
                                                  'biomass', 'carbon'
                                                )]

#    
#   back transform the log to canonical postion, then determine the inclusion zone perimeter...
#
    izPerim = downLog@rotLog
    dn = dimnames(izPerim)
    izPerim = izPerim  %*% trMatInv                                    #to canonical log position
    dimnames(izPerim) = dn
    izPerim[,'y'] = switch(pdsType,                                     #y==radius, already in feet or meters
                           volume = {
                                     sense = sign(izPerim[,'y'])        #preserve the sign in squaring...
                                     sense*pds@kpds * pi*izPerim[,'y']*izPerim[,'y']
                                    },
                           surfaceArea = pds@kpds * 2*pi*izPerim[,'y'], #remember y=radius, not diameter
                           coverageArea = pds@kpds * 2*izPerim[,'y']    #remember y=radius, not diameter
                          )
    izPerim = izPerim %*% trMat
    dimnames(izPerim) = dn

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
    pdsIZ = new('perpendicularDistanceIZ', downLog=downLog,
                 pds = pds,                          #pds sampling object
                 pdsType = pdsType,                  #PPS version of PDS
                 izPerim = izPerim,                  #matrix representation of perimeter
                 perimeter = spObj,                  #SpatialPolygons perimeter
                 pgArea = pgArea,                   #area of SpatialPolygons izone: approximate
                 spUnits = spUnits,                  #CRS units
                 description = description,          #a comment
                 units = downLog@units,              #units of measure
                 area = izArea,                      #exact inclusion zone area
                 puaBlowup = puaBlowup,              #sausage per unit area blowup factor
                 puaEstimates = puaEstimates,        #per unit area estimates
                 bbox = bbox(spObj)                  #overall bounding box--redundant here
                )

    return(pdsIZ)
}   #perpendicularDistanceIZ constructor
)   #setMethod


    




    


       
#================================================================================
#  7. method for class omnibusPDSIZ construction--just call parent method...
#
setMethod('omnibusPDSIZ',
          signature(downLog = 'downLog', pds = 'perpendicularDistance'),
function(downLog,
         pds,
         description = 'inclusion zone for down log omnibus perpendicular distance sampling',
         spID = paste('opds',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         pdsType = .StemEnv$pdsTypes,
         ...
        )
{
#------------------------------------------------------------------------------
#
#   create an object of the correct class...
#
    opdsIZ = perpendicularDistanceIZ(downLog, pds, description = description,
                                     spID = spID, spUnits = spUnits,
                                     #pdsClass = 'omnibusPDSIZ',
                                     pdsType = pdsType,
                                     ...
                                    )
    opdsIZ = as(opdsIZ, 'omnibusPDSIZ')   #cast

#
#   the omnibus pua estimates all depend on perpendicular diameter, so there are no
#   overall estimates to assign...
#
    puaEstimates = as.list(rep(NA, 7))
    names(puaEstimates) = .StemEnv$puaEstimates[c('volume', 'Density', 'Length',
                                                     'surfaceArea', 'coverageArea',
                                                     'biomass', 'carbon'
                                                   )]
    opdsIZ@puaEstimates = puaEstimates

    return(opdsIZ)
}   #omnibusPDSIZ constructor
)   #setMethod


    




    

    



       
#================================================================================
#
#  8. DLPDS...
#
#  method for class distanceLimitedPDSIZ construction using a "distanceLimited"
#  object...
#
#================================================================================
#
setClassUnion('dlsNumeric', c('distanceLimited','numeric'))       #let this constructor handle both
setMethod('distanceLimitedPDSIZ',
          signature(downLog = 'downLog', pds = 'perpendicularDistance', dls='dlsNumeric'), 
function(downLog,
         pds,
         dls,
         description = 'inclusion zone for down log distance limited PDS',
         spID = paste('dlpds',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         pdsType = .StemEnv$pdsTypes,
         ...
        )
{
#------------------------------------------------------------------------------
#
#   The flow of this routine is as follows...
#     1. Check to see if the DL affects the log, if not, just return a regular
#        PDS object since that is what it reduces to
#
#        Note: If the DL is such that the entire log diameter is too small to
#              have any distances beyond this, then we stop here as the entire
#              log is regular PDS; the other slots are NULLs
#        Otherwise...
#
#     2. Determine the part of the log that is affected by the distance limit
#        and create a DLMCIZ object with this component
#
#        Note: the whole log can be under DL if both end diameters are large 
#              enough, in that case, we stop here and fill other slots with
#              NULLs
#        Otherwise...
#
#     3. Determine the part of the log that is treated as a regular PDS object
#        and create this object
#     4. Merge the two components in steps 2 & 3 into an overall distance
#        limited PDS object
#
#   8&10-Mar-2011
#
#------------------------------------------------------------------------------
#  
#   convert distance limit if necessary...
#
    if(is.numeric(dls))
      dls = distanceLimited(dls, units=downLog@units, description=description, ...)
    distanceLimit = dls@distanceLimit
    pdsType = match.arg(pdsType)

    
#-------------------------------------------------------------------------------
#   Step 1:   See if DLPDS is required...
#
#   first just see if the entire log is already in, if so, no further work is
#   necessary, create a regular pds object...
#
    full.pdsIZ = perpendicularDistanceIZ(downLog, pds, description = description,
                                         spID = spID, spUnits = spUnits,
                                         pdsClass = 'perpendicularDistanceIZ',
                                         pdsType = pdsType,
                                         ...
                                        )
    kpds = pds@kpds
    dl.diam = switch(pdsType,                                         #DL cutoff diameter
                     volume = sqrt(4*distanceLimit/(kpds*pi)),   
                     surfaceArea = distanceLimit/(kpds*pi),
                     coverageArea = distanceLimit/kpds
                    )
#
#   done if the whole log is pds, just set up the object...
#
    if(dl.diam > max(downLog@taper$diameter)) {        #off the large end--whole log is pds
      dlpds = as(full.pdsIZ, 'distanceLimitedPDSIZ')   #cast object
      dlpds@dls = dls
      dlpds@dlsDiameter = dl.diam                      #limiting diameter: can be off the log
      dlpds@pdsPart = NULL
      dlpds@dlsPart = NULL
      dlpds@pdsFull = full.pdsIZ                       #just for consistency, does take more room
      return(dlpds)
    }

    full.centerOffset = coordinates(downLog@location)[1,]    #for the full log
    halfLen = downLog@logLen/2
    logAngle = downLog@logAngle
    trMat = transfMatrix(logAngle)



#-------------------------------------------------------------------------------
#   Step 2:    DL component...
#
#   now do the butt end of the log (or whole log), which has a rectangular
#   inclusion zone...
#
    if(dl.diam <= min(downLog@taper$diameter)) {                    #off small end--whole log==DL inclusion zone
      dl.taper = downLog@taper
      dl.centerOffset = full.centerOffset
      dl.length = downLog@logLen
    }
    else {                                                          #someplace on the log
      dl.length = taperInterpolate(downLog,'length', dl.diam)       #length to DL cutoff
      dl.taper = downLog@taper[downLog@taper$diameter >= dl.diam,]  #taper below DL cutoff
      dl.taper = rbind(dl.taper, c(dl.diam, dl.length))             #with new "top"
      os.len = halfLen - dl.length/2                                #offset length to new center point from old
      dl.centerOffset = full.centerOffset - (matrix(c(os.len,0,1),nr=1) %*% trMat)[,-3]
      names(dl.centerOffset) = c('x','y')
    }

    dl.dlog = downLog(dl.taper,
                      solidType = NULL,                             #NULL for taper data
                      logAngle = logAngle,
                      vol2wgt = downLog@conversions[['volumeToWeight']],    #[[]]to avoid concatenating names
                      wgt2carbon = downLog@conversions[['weightToCarbon']], #[[]]to avoid concatenating names
                      centerOffset = dl.centerOffset,
                      species = downLog@species,
                      logID = downLog@spLog@polygons$pgsLog@ID,
                      description = downLog@description,
                      units = downLog@units,
                      spUnits = downLog@spUnits,
                      runQuiet=TRUE
                     )

      
#
#   now create just the distance limited component of the inclusion zone from this truncated
#   smaller down log butt portion; this will be a distanceLimitedMCIZ object...
#
    dls = distanceLimited(distanceLimit, units=downLog@units)
    dlIZ = distanceLimitedIZ(dl.dlog, dls, description = description,
                             spID = spID, spUnits = spUnits, ...)

#
#   if the whole log is DL, then we are done, make the full object here and exit...
#
    if(dl.diam <= min(downLog@taper$diameter)) {        #off small end--whole log==DL inclusion zone
      dlpdsIZ = new('distanceLimitedPDSIZ', downLog=downLog,
                    pds = pds,                          #pds sampling object
                    pdsType = pdsType,                  #PPS version of PDS
                    izPerim = dlIZ@izPerim,             #matrix representation of perimeter
                    perimeter = dlIZ@perimeter,         #SpatialPolygons perimeter
                    pgArea = dlIZ@pgArea,               #area of SpatialPolygons izone: approximate
                    spUnits = spUnits,                  #CRS units
                    description = description,          #a comment
                    units = downLog@units,              #units of measure
                    area = dlIZ@area,                   #exact inclusion zone area
                    puaBlowup = dlIZ@puaBlowup,         #NA per unit area blowup factor
                    puaEstimates = dlIZ@puaEstimates,   #per unit area estimates
                    bbox = bbox(dlIZ),                  #overall bounding box--redundant here
                    dls = dls,                          #distanceLimited object
                    dlsDiameter = dl.diam,              #the limiting diameter
                    pdsPart = NULL,                     #pdsIZ component object
                    dlsPart = dlIZ,                     #dlmcsIZ component object
                    pdsFull = full.pdsIZ                #as if it were plain pdsIZ              
                 ) 
      return(dlpdsIZ)
    }

    
#-------------------------------------------------------------------------------
#   Step 3:   PDS component...
#
#   okay, if we got to here, then both components exist in the log; do the top (pds)
#   portion here: create a downLog object from the top part...
#
    pds.taper = downLog@taper[downLog@taper$diameter < dl.diam,]    #taper to DL cutoff
    pds.taper = rbind(c(dl.diam, dl.length), pds.taper)             #with new "butt"
    pds.taper$length = pds.taper$length - dl.length                 #translate to butt at zero

    os.len = dl.length + diff(range(pds.taper$length))/2 - halfLen #offset length to new center point from old
    pds.centerOffset = full.centerOffset + (matrix(c(os.len,0,1),nr=1) %*% trMat)[,-3]
    names(pds.centerOffset) = c('x','y')
    pds.dlog = downLog(pds.taper,
                       solidType = NULL, #downLog@solidType,
                       logAngle = logAngle,
                       vol2wgt = downLog@conversions[['volumeToWeight']],    #[[]]to avoid concatenating names
                       wgt2carbon = downLog@conversions[['weightToCarbon']], #[[]]to avoid concatenating names
                       centerOffset = pds.centerOffset,
                       species = downLog@species,
                       logID = downLog@spLog@polygons$pgsLog@ID,
                       description = downLog@description,
                       units = downLog@units,
                       spUnits = downLog@spUnits,
                       runQuiet=TRUE
                      )
 
#
#   now create just the pds component of the inclusion zone from this truncated
#   smaller down log...
#
    pdsIZ = perpendicularDistanceIZ(pds.dlog, pds, description = description,
                                     spID = spID, spUnits = spUnits,
                                     pdsClass = 'perpendicularDistanceIZ',
                                     pdsType = pdsType,
                                     ...
                                    )
    



#-------------------------------------------------------------------------------
#   Step 4:   Complete object...
#
#   Combine the two pieces into one IZ object...
#
#   direct calculation of the full inclusion zone area and blowup factor...
#
    izArea = pdsIZ@area + dlIZ@area                #for any pdsType
    puaBlowup = NA_real_                           #no true, single blowup for dlpds
    
#
#   per unit area estimates for the full log...
#
    puaEstimates = as.list(unlist(pdsIZ@puaEstimates) + unlist(dlIZ@puaEstimates))  #not in general
    
#
#   full perimeter is in two parts that must be joined together...
#
    nn = nrow(pdsIZ@izPerim)
    pds.izPerim = pdsIZ@izPerim                     #already in correct orientation, etc.

    izPerim = rbind(dlIZ@izPerim[1:2,], pds.izPerim[-nn,], dlIZ@izPerim[3:5,])  #without repeat in pds
    dimnames(izPerim) = list(NULL,c('x','y','hc'))   

#
#   and make a SpatialPolygons object...
#
    pg = Polygon(izPerim[,-3])                             #sans hc
    pgs = Polygons(list(pg=pg), ID=spID)
    spObj = SpatialPolygons(list(pgs=pgs), proj4string = spUnits )
    pgArea = pg@area

#
#   create the distance limited IZ object...
#
    dlpdsIZ = new('distanceLimitedPDSIZ', downLog=downLog,
                  pds = pds,                          #pds sampling object
                  pdsType = pdsType,                  #PPS version of PDS
                  izPerim = izPerim,                  #matrix representation of perimeter
                  perimeter = spObj,                  #SpatialPolygons perimeter
                  pgArea = pgArea,                    #area of SpatialPolygons izone: approximate
                  spUnits = spUnits,                  #CRS units
                  description = description,          #a comment
                  units = downLog@units,              #units of measure
                  area = izArea,                      #exact inclusion zone area
                  puaBlowup = puaBlowup,              #NA per unit area blowup factor
                  puaEstimates = puaEstimates,        #per unit area estimates
                  bbox = bbox(spObj),                 #overall bounding box--redundant here
                  dls = dls,                          #distanceLimited object
                  dlsDiameter = dl.diam,              #the limiting diameter
                  pdsPart = pdsIZ,                    #pdsIZ component object
                  dlsPart = dlIZ,                     #dlmcsIZ object
                  pdsFull = full.pdsIZ                #as if it were plain pdsIZ              
                 )
    

    
    return(dlpdsIZ)
}   #distanceLimitedPDSIZ constructor
)   #setMethod

    




    


       
#================================================================================
#  9. method for class omnibusDLPDSIZ construction--just call parent method...
#
setMethod('omnibusDLPDSIZ',
          signature(downLog = 'downLog', pds = 'perpendicularDistance', dls = 'dlsNumeric'),
function(downLog,
         pds,
         dls,
         description = 'inclusion zone for down log omnibus distance limited PDS',
         spID = paste('odlpds',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         pdsType = .StemEnv$pdsTypes,
         ...
        )
{
#------------------------------------------------------------------------------
#
#   create an object of the correct class; note that we cannot use as() to
#   directly coerce below because it chokes on the classUnions with NULL, so
#   we can just create an object with new() instead...
#
    pdsType = match.arg(pdsType)
    dlpdsIZ = distanceLimitedPDSIZ(downLog, pds, dls, description = description,
                                     spID = spID, spUnits = spUnits,
                                     pdsType = pdsType,
                                     ...
                                    )

#
#   this just mirrors what we have in omnibusPDSIZ and distanceLimitedMCIZ
#   using coercion without having to recreate everything from scratch...
#
    npua = length(dlpdsIZ@puaEstimates)

    #pdsPart...
    if(!is.null(dlpdsIZ@pdsPart)) {
      pdsPart = as(dlpdsIZ@pdsPart, 'omnibusPDSIZ')   #cast/coerce
      pdsPart@puaEstimates[1:npua] = rep(NA, npua)    #will keep names okay
    }
    else
      pdsPart = NULL

    #dlsPart...
    if(!is.null(dlpdsIZ@dlsPart)) {
      dlsPart = as(dlpdsIZ@dlsPart, 'distanceLimitedMCIZ')   #cast/coerce
      unitArea = ifelse(dlsPart@downLog@units==.StemEnv$msrUnits$English, .StemEnv$sfpAcre, .StemEnv$smpHectare)
      izArea = dlsPart@area
      logLen = dlsPart@downLog@logLen
      dlsPart@puaBlowup = unitArea/izArea * logLen     #<<<*****Note
      dlsPart@puaEstimates[c('volume','surfaceArea','coverageArea','biomass', 'carbon')] = NA
    }
    else
      dlsPart = NULL

    #pdsFull...
    if(!is.null(dlpdsIZ@pdsFull)) {
      pdsFull = as(dlpdsIZ@pdsFull, 'omnibusPDSIZ')   #cast/coerce
      pdsFull@puaEstimates[1:npua] = rep(NA, npua)    #will keep names okay
    }
    else
      pdsFull = NULL

#
#   per unit area estimates for the full log are unknown at this point as they
#   depend on some function of perpendicular diameter...
#
    puaEstimates = as.list(rep(NA, 7))             #all to NA to cover omnibusDLPDSIZ
    names(puaEstimates) = .StemEnv$puaEstimates[c('volume', 'Density', 'Length',
                                                  'surfaceArea', 'coverageArea',
                                                  'biomass', 'carbon'
                                                )]
    
    
#
#   create the new object afresh...
#
    odlpdsIZ = new('omnibusDLPDSIZ', downLog=downLog,
                   pds = pds,                            #pds sampling object
                   pdsType = pdsType,                    #PPS version of PDS
                   izPerim = dlpdsIZ@izPerim,            #matrix representation of perimeter
                   perimeter = dlpdsIZ@perimeter,        #SpatialPolygons perimeter
                   pgArea = dlpdsIZ@pgArea,              #area of SpatialPolygons izone: approximate
                   spUnits = spUnits,                    #CRS units
                   description = description,            #a comment
                   units = dlpdsIZ@units,                #units of measure
                   area = dlpdsIZ@area,                  #exact inclusion zone area
                   puaBlowup = dlpdsIZ@puaBlowup,        #NA per unit area blowup factor
                   puaEstimates = puaEstimates,          #per unit area estimates
                   bbox = dlpdsIZ@bbox,                  #overall bounding box--redundant here
                   dls = dlpdsIZ@dls,                    #distanceLimited object
                   dlsDiameter = dlpdsIZ@dlsDiameter,    #the limiting diameter
                   pdsPart = pdsPart,                    #pdsIZ component object
                   dlsPart = dlsPart,                    #DL component pdsIZ object
                   pdsFull = pdsFull                     #as if it were plain pdsIZ
                  )


    return(odlpdsIZ)
}   #omnibusDLPDSIZ constructor
)   #setMethod


    


       
#================================================================================
#  10. method for class hybridDLPDSIZ construction--just call parent method...
#
#     this code is very similar to omnibusDLPDSIZ, which was written first, here
#     we have distanceLimitedMCIZ + canonicalPDSIZ iz components...
#
setMethod('hybridDLPDSIZ',
          signature(downLog = 'downLog', pds = 'perpendicularDistance', dls = 'dlsNumeric'),
function(downLog,
         pds,
         dls,
         description = 'inclusion zone for down log hybrid distance limited PDS',
         spID = paste('hdlpds',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         pdsType = .StemEnv$pdsTypes,
         ...
        )
{
#------------------------------------------------------------------------------
#
#   create an object of the correct class; note that we cannot use as() to
#   directly coerce below because it chokes on the classUnions with NULL, so
#   we can just create an object with new() instead...
#
    pdsType = match.arg(pdsType)
    dlpdsIZ = distanceLimitedPDSIZ(downLog, pds, dls, description = description,
                                     spID = spID, spUnits = spUnits,
                                     pdsType = pdsType,
                                     ...
                                    )

#
#   this just mirrors what we have in omnibusPDSIZ and distanceLimitedMCIZ
#   using coercion without having to recreate everything from scratch;
#   note that the pds portion is canonical, so we leave it alone...
#
    npua = length(dlpdsIZ@puaEstimates)

    #dlsPart...
    if(!is.null(dlpdsIZ@dlsPart)) {
      dlsPart = as(dlpdsIZ@dlsPart, 'distanceLimitedMCIZ')   #cast/coerce
      unitArea = ifelse(dlsPart@downLog@units==.StemEnv$msrUnits$English, .StemEnv$sfpAcre, .StemEnv$smpHectare)
      izArea = dlsPart@area
      logLen = dlsPart@downLog@logLen
      dlsPart@puaBlowup = unitArea/izArea * logLen     #<<<*****Note
      dlsPart@puaEstimates[c('volume','surfaceArea','coverageArea','biomass', 'carbon')] = NA
    }
    else
      dlsPart = NULL

#
#   per unit area estimates for the full log are unknown at this point as they
#   depend on some function of perpendicular diameter...
#
    puaEstimates = as.list(rep(NA, 7))             #all to NA to cover hybridDLPDSIZ
    names(puaEstimates) = .StemEnv$puaEstimates[c('volume', 'Density', 'Length',
                                                  'surfaceArea', 'coverageArea',
                                                  'biomass', 'carbon'
                                                )]
    
    
#
#   create the new object afresh...
#
    hdlpdsIZ = new('hybridDLPDSIZ', downLog=downLog,
                   pds = pds,                            #pds sampling object
                   pdsType = pdsType,                    #PPS version of PDS
                   izPerim = dlpdsIZ@izPerim,            #matrix representation of perimeter
                   perimeter = dlpdsIZ@perimeter,        #SpatialPolygons perimeter
                   pgArea = dlpdsIZ@pgArea,              #area of SpatialPolygons izone: approximate
                   spUnits = spUnits,                    #CRS units
                   description = description,            #a comment
                   units = dlpdsIZ@units,                #units of measure
                   area = dlpdsIZ@area,                  #exact inclusion zone area
                   puaBlowup = dlpdsIZ@puaBlowup,        #NA per unit area blowup factor
                   puaEstimates = puaEstimates,          #per unit area estimates
                   bbox = dlpdsIZ@bbox,                  #overall bounding box--redundant here
                   dls = dlpdsIZ@dls,                    #distanceLimited object
                   dlsDiameter = dlpdsIZ@dlsDiameter,    #the limiting diameter
                   pdsPart = dlpdsIZ@pdsPart,            #pdsIZ component object
                   dlsPart = dlsPart,                    #DL component pdsIZ object
                   pdsFull = dlpdsIZ@pdsFull             #as if it were plain pdsIZ
                  )


    return(hdlpdsIZ)
}   #hybridDLPDSIZ constructor
)   #setMethod











       
#================================================================================
# 11. method for functions and class distanceLimitedIZ...
#
setMethod('distanceLimitedIZ',
          signature(downLog = 'downLog', dls = 'distanceLimited'), #change second argument!!
function(downLog,
         dls,
         description = 'inclusion zone for down log DL sampling',
         spID = paste('dl',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         ...
        )
{
#------------------------------------------------------------------------------
#
#   a quick check...
#
    if(downLog@units != dls@units)
      stop('units are not the same for downLog and dls!')

    
#
#   put in a warning for number of points in taper data frame too few <20 for now<<<<<******
#
#   transformation matrix...
#
    centerOffset = coordinates(downLog@location)
    logAngle = downLog@logAngle
    trMat = transfMatrix(logAngle, centerOffset)
    trMatInv = solve(trMat)
    
#
#   direct calculation of the inclusion zone area and blowup factor...
#
    distanceLimit = dls@distanceLimit
    logLen = downLog@logLen
    unitArea = ifelse(downLog@units==.StemEnv$msrUnits$English, .StemEnv$sfpAcre, .StemEnv$smpHectare) 
    izArea = 2*logLen*distanceLimit          #DL component
    puaBlowup = unitArea/izArea 
    puaEstimates = list(downLog@logVol*puaBlowup, puaBlowup, logLen*puaBlowup,
                        downLog@surfaceArea*puaBlowup, downLog@coverageArea*puaBlowup,
                        downLog@biomass*puaBlowup, downLog@carbon*puaBlowup
                       )
    names(puaEstimates) = .StemEnv$puaEstimates[c('volume', 'Density', 'Length',
                                                  'surfaceArea', 'coverageArea',
                                                  'biomass', 'carbon'
                                                )]


#    
#   back transform the log to canonical postion, then determine the inclusion zone perimeter...
#
    halfLen = logLen/2
    izPerim = matrix(c(-halfLen, distanceLimit, 1), nrow=1)          #left-most top point
    izPerim = rbind(izPerim, c(halfLen, distanceLimit,1))            #right-most top point
    izPerim = rbind(izPerim, c(halfLen, -distanceLimit,1))           #right-most bottom point
    izPerim = rbind(izPerim, c(-halfLen, -distanceLimit,1))          #left-most bottom point
    izPerim = rbind(izPerim, izPerim[1,])                            #closed polygon
    trMat = transfMatrix(logAngle, centerOffset)
    izPerim = izPerim %*% trMat
    dimnames(izPerim) = dimnames(downLog@rotLog)


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
    dlsIZ = new('distanceLimitedIZ',
                 downLog=downLog,
                 dls = dls,                          #distanceLimited sampling object
                 izPerim = izPerim,                  #matrix representation of perimeter
                 perimeter = spObj,                  #SpatialPolygons perimeter
                 pgArea = pgArea,                    #area of SpatialPolygons izone: approximate
                 spUnits = spUnits,                  #CRS units
                 description = description,          #a comment
                 units = downLog@units,              #units of measure
                 area = izArea,                      #exact inclusion zone area
                 puaBlowup = puaBlowup,              #sausage per unit area blowup factor
                 puaEstimates = puaEstimates,        #per unit area estimates
                 bbox = bbox(spObj)                  #overall bounding box--redundant here
                )

    return(dlsIZ)
}   #distanceLimitedIZ constructor
)   #setMethod








       
#================================================================================
#  12. method for functions and class distanceLimitedMCIZ; everything is the same
#      as under dls, except the per unit area estimates for all but Density
#      and Length...
#
setMethod('distanceLimitedMCIZ',
          signature(downLog = 'downLog', dls = 'distanceLimited'), #change second argument!!
function(downLog,
         dls,
         description = 'inclusion zone for down log DLMC sampling',
         spID = paste('dlmc',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         ...
        )
{
#------------------------------------------------------------------------------
#
    dlsIZ = distanceLimitedIZ(downLog=downLog, dls=dls, description=description,
                              spID=spID, spUnits=spUnits, ...)
    dlmcIZ = as(dlsIZ, 'distanceLimitedMCIZ')   #cast

#
#   we must reset the blowup for use in the InclusionZoneGrid routine for dlmc...
#
    unitArea = ifelse(downLog@units==.StemEnv$msrUnits$English, .StemEnv$sfpAcre, .StemEnv$smpHectare)
    izArea = dlsIZ@area
    logLen = dlsIZ@downLog@logLen
    dlmcIZ@puaBlowup = unitArea/izArea * logLen     #<<<*****Note
    
#
#   per unit area estimates...
#     most of the per unit area estimates for this component depend on the perpendicular
#     diameter in some form, but a couple (Density and Length) are constant...
#
    dlmcIZ@puaEstimates[c('volume','surfaceArea','coverageArea','biomass', 'carbon')] = NA
    
    return(dlmcIZ)
}   #distanceLimitedMCIZ constructor
)   #setMethod
    




#---------------------------------------------------------------------------
#
#   This section defines the constructors for standTreeIZ objects...
#
#Author...									Date: 1-Dec-2011
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#

#if(!isGeneric("circularPlotIZ")) 
  setGeneric('circularPlotIZ',  
             function(standingTree, plotRadius, ...) standardGeneric('circularPlotIZ'),
             signature = c('standingTree', 'plotRadius')
            )

#if(!isGeneric("horizontalPointIZ")) 
  setGeneric('horizontalPointIZ',  
             function(standingTree, angleGauge, ...) standardGeneric('horizontalPointIZ'),
             signature = c('standingTree', 'angleGauge')
            )
   



       
#================================================================================
#  1. method for functions and class circularPlotIZg...
#
setMethod('circularPlotIZ',
          signature(standingTree = 'standingTree', plotRadius = 'numeric'),
function(standingTree,
         plotRadius,
         description = 'inclusion zone for circular plot method',
         spID = paste('cp',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         ...
        )
{
#------------------------------------------------------------------------------
#
#   get bbox from the standingTree object...
#
    standingTree.bbox = bbox(standingTree@spDBH)
    
#
#   get bbox from the circularPlot object---tied to the tree center...
#   --we must be careful below and not just pass the "..." because we are
#     explicitly passing "units" and if the user also put "units" in the "..."
#     call to the current method, then it will be passed twice throwing an
#     error; so pull out any other arguments from "..." that could go to
#     circularPlot and pass them explicitly or set them to default; i.e.,
#     nptsPerimeter
#
    units = standingTree@units
    loc = coordinates(standingTree@location)[1,]  #need a vector from matrix
    args = list(...)
    if('nptsPerimeter' %in% names(args))          #user can override here
      nptsPerimeter = args$nptsPerimeter          #user's choice
    else
      nptsPerimeter = 100                         #default in circularPlot()
    circularPlot = circularPlot(plotRadius, units=units, centerPoint=loc,
                                spID = spID, spUnits = spUnits,
                                nptsPerimeter = nptsPerimeter)
    cp.bbox = bbox(circularPlot@perimeter)

#
#   combine them for the overall bbox (we can't assume that the circular plot
#   will be larger than the tree)...
#
    min = apply(cbind(standingTree.bbox, cp.bbox), 1, min)
    max = apply(cbind(standingTree.bbox, cp.bbox), 1, max)
    bbox = matrix(cbind(min,max),ncol=2, dimnames=list(c('x','y'), c('min','max')))

#
#   per unit area estimates...
#
    baFactor =  ifelse(standingTree@units==.StemEnv$msrUnits$English, .StemEnv$baFactor['English'],
                       .StemEnv$baFactor['metric'])
    unitArea = ifelse(standingTree@units==.StemEnv$msrUnits$English, .StemEnv$sfpAcre, .StemEnv$smpHectare) 
    puaBlowup = unitArea/circularPlot@area 
    puaEstimates = list(standingTree@treeVol*puaBlowup,
                        puaBlowup,
                        standingTree@dbh^2*baFactor*puaBlowup,
                        standingTree@surfaceArea*puaBlowup, 
                        standingTree@biomass*puaBlowup,
                        standingTree@carbon*puaBlowup
                       )
    names(puaEstimates) = .StemEnv$puaEstimates[c('volume', 'Density', 'basalArea',
                                                  'surfaceArea', 'biomass', 'carbon'
                                                )]

#
#   create the object...
#
    cpIZ = new('circularPlotIZ', standingTree=standingTree, circularPlot=circularPlot,
               bbox=bbox, spUnits=spUnits, description=description,
               units = units, puaBlowup = puaBlowup, puaEstimates = puaEstimates
              )

    return(cpIZ)
}   #circularPlotIZ constructor
)   #setMethod





       
#================================================================================
#  2. method for functions and class horizontalPointIZ...
#
setMethod('horizontalPointIZ',
          signature(standingTree = 'standingTree', angleGauge = 'angleGauge'),
function(standingTree,
         angleGauge,
         description = 'inclusion zone for horizontal point sampling method',
         spID = paste('hps',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         ...
        )
{
#------------------------------------------------------------------------------
#
#   get bbox from the standingTree object...
#
    standingTree.bbox = bbox(standingTree@spDBH)
    
#
#   get bbox from the circularPlot object---tied to the tree center...
#   --we must be careful below and not just pass the "..." because we are
#     explicitly passing "units" and if the user also put "units" in the "..."
#     call to the current method, then it will be passed twice throwing an
#     error; so pull out any other arguments from "..." that could go to
#     circularPlot and pass them explicitly or set them to default; i.e.,
#     nptsPerimeter
#
    units = standingTree@units                    #user has no choice here, must match
    dbh = standingTree@dbh
    plotRadius = angleGauge@PRF * dbh             #remember, dbh is stored in m or ft, not cm or inches
    loc = coordinates(standingTree@location)[1,]  #need a vector from matrix
    args = list(...)
    if('nptsPerimeter' %in% names(args))          #user can override here
      nptsPerimeter = args$nptsPerimeter          #user's choice
    else
      nptsPerimeter = 100                         #default in circularPlot()
    circularPlot = circularPlot(plotRadius, units=units, centerPoint=loc,
                                spID = spID, spUnits = spUnits,
                                nptsPerimeter = nptsPerimeter)
    cp.bbox = bbox(circularPlot@perimeter)

#
#   combine them for the overall bbox (we can't assume that the circular plot
#   will be larger than the tree)...
#
    min = apply(cbind(standingTree.bbox, cp.bbox), 1, min)
    max = apply(cbind(standingTree.bbox, cp.bbox), 1, max)
    bbox = matrix(cbind(min,max),ncol=2, dimnames=list(c('x','y'), c('min','max')))

#
#   per unit area estimates...
#
    baFactor =  ifelse(standingTree@units==.StemEnv$msrUnits$English, .StemEnv$baFactor['English'],
                       .StemEnv$baFactor['metric'])
    #unitArea = ifelse(standingTree@units==.StemEnv$msrUnits$English, .StemEnv$sfpAcre, .StemEnv$smpHectare)
    ba = baFactor*dbh*dbh                               #remember, dbh is stored in m or ft, not cm or inches
    baf = angleGauge@baf
    puaBlowup = baf/ba                                  #unitArea constant is in the baf already
    puaEstimates = list(standingTree@treeVol*puaBlowup,
                        puaBlowup,
                        baf, 
                        standingTree@surfaceArea*puaBlowup, 
                        standingTree@biomass*puaBlowup,
                        standingTree@carbon*puaBlowup
                       )
    names(puaEstimates) = .StemEnv$puaEstimates[c('volume', 'Density', 'basalArea',
                                                  'surfaceArea', 'biomass', 'carbon'
                                                )]

#
#   create the object...
#
    hpIZ = new('horizontalPointIZ',
               standingTree = standingTree,
               circularPlot = circularPlot,
               angleGauge = angleGauge,
               bbox = bbox,
               spUnits = spUnits,
               description = description,
               units = units,
               puaBlowup = puaBlowup,
               puaEstimates = puaEstimates
              )

    return(hpIZ)
}   #horizontalPointIZ constructor
)   #setMethod

