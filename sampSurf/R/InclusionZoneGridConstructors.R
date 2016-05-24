#---------------------------------------------------------------------------
#
#   This file holds the S4 definition for the constructor methods of the
#   InclusionZoneGrid class...
#
#   The methods include...
#     0. matrix constructor
#
#        ...downLogIZ subclass constructors...
#
#     1. a constructor for 'standUpIZ'
#     2. for 'sausageIZ'
#     3. for 'chainSawIZ'
#     4. for 'fullChainSawIZ'
#     5. for 'pointRelascopeIZ'
#     6. for 'perpendicularDistanceIZ'
#     7. for 'omnibusPDSIZ'
#     8. for 'distanceLimitedPDSIZ'
#     9. for 'distanceLimitedMCIZ'
#    10. for 'distanceLimtedIZ'
#
#        ...standingTreeIZ subclass constructors...
#
#     1. for 'circularPlotIZ'
#        for 'horizontalPointIZ', nothing is required, it uses 'circularPlotIZ'
#
#   Note: I adopted a new strategy after HPS was added, new sampling methods
#         have all their respective code for all classes in one file. 
#         Please see those files for their InclusionZoneGrid constructors.
#
#   Of course, if chainSaw was not so strange, we would only require one
#   method for every different technique! (Well except for omnibus, and
#   other methods with varying height IZs...)
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
#   generic definition...
#
setGeneric('izGrid',  
           function(izObject, tract, ...) standardGeneric('izGrid'),
             signature = c('izObject', 'tract')
            )




          
#================================================================================
#  This is a helper routine, not normally called from outside the construction
#  of an InclusionZoneGrid object...
#
#  0. method for 'matrix' and 'Tract' classes; this is used within the following
#     constructors and should not normally be called by itself as it simply
#     established the minimal bounding grid that encompasses the inclusion zone
#
setMethod('izGrid',
          signature(izObject = 'matrix', tract='Tract'),
function(izObject,        #a bbox object
         tract,
         data = 0,        #background value for all cells in the minimal bounding grid
         useCrop = TRUE,  #T: use crop(); F: use extent method
         ...
        )
{
#---------------------------------------------------------------------------
#
#   need to get the inclusion zone bbox aligned with the tract grid...
#
    if(!useCrop) {
      bbex = extent(izObject)
      j = intersectExtent(bbex, extent(tract) )   #get extent of overlap    !!!>this has been replaced by intersect()
      jj = alignExtent(j, tract)                  #and align to tract grid
    }
    else
      jj = extent(izObject)

#
#   pad one extra cell on each side for good measure...
#
    cr = xres(tract)            #square cells always
    jj@xmin = jj@xmin - cr
    jj@xmax = jj@xmax + cr
    jj@ymin = jj@ymin - cr
    jj@ymax = jj@ymax + cr

#
#   now create a raster object large enough to surround the inclusion zone bbox
#   out of the above extent...
#
#   to use crop() instead of the above, the entire IZ must lie within the tract...
#
    if(!useCrop) {
      nx = (jj@xmax - jj@xmin)/cr                   #remember x=columns
      ny = (jj@ymax - jj@ymin)/cr                   #and y=rows
      sr = raster(jj, nrows=as.integer(ny), ncols=as.integer(nx), crs=NA)
    }
    else
      sr = crop(tract, jj)
    sr[] = rep(data, ncell(sr))

    return(sr)
}   #izGrid base method
)   #setMethod



  

  


          
#================================================================================
#  1. method for 'standUpIZ' and 'Tract' classes...
#
setMethod('izGrid',
          signature(izObject = 'standUpIZ', tract='Tract'),
function(izObject,
         tract,
         description = 'standUpIZ inclusion zone grid object',
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
}   #izGrid for'standUpIZ'
)   #setMethod


  

  


          
#================================================================================
#  2. method for 'sausageIZ' and 'Tract' classes...
#
setMethod('izGrid',
          signature(izObject = 'sausageIZ', tract='Tract'),
function(izObject,
         tract,
         description = 'sausageIZ inclusion zone grid object',
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
}   #izGrid for'sausageIZ'
)   #setMethod




  
 


          
#================================================================================
#  3. method for 'chainSawIZ' and 'Tract' classes--of course chainsaw is weird
#     since it is really a point inclusion zone, so it must have a more involved
#     constructor method...
#
setMethod('izGrid',
          signature(izObject = 'chainSawIZ', tract='Tract'),
function(izObject,
         tract,
         description = 'chainSaw grid point inclusion zone grid object',
         wholeIZ = TRUE,           #TRUE: grid the whole object; FALSE: just grid the IZ
         ...
        )
{
#---------------------------------------------------------------------------
#
#   find the wholeIZ bounding grid to begin with...
#
    iz.bbox = bbox(izObject)                           #grid point plus log (always internal)    
    izg = izGrid(iz.bbox, tract)                       #get the minimal bounding grid
    cpt = perimeter(izObject, whatSense='point')       #circularplot centerpoint only
    grid = rasterize(cpt,                              #plot center location
                     izg,                              #bounding grid
                     fun = function(x,na.rm) 0,        #value at x,y--anything other than NA
                     background = NA                   #all other values are zero
                    )

#
#   now, if we are just interested in the single grid cell inclusion zone for the plot
#   centerpoint, we must do a little more work to narrow it down within the properties
#   of the raster grid...
#
    if(!wholeIZ) {
      j = cellFromXY(grid, cpt)        #cell number for center point of the circular plot
      xy = xyFromCell(grid, j)         #cell center point, may not be same as cpt
      res = xres(grid)/2               #half resolution of the grid cells
      ex = extent(c(xy[1,'x']-res, xy[1,'x']+res,  #extent for one grid cell which contains
                  xy[1,'y']-res, xy[1,'y']+res))   #the circular plot center point
      grid = crop(grid, ex)            #ensures the result is on the same grid as parent
      }

#
#   a data frame with each pua estimate...
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
}   #izGrid for'chainSawIZ'
)   #setMethod

    


 


          
#================================================================================
#  4. method for 'fullChainSawIZ' within the encompassing "sausage" inclusion
#     zone and 'Tract' classes...
#
#     Recall that 'fullChainSawIZ' is a direct subclass to 'sausageIZ' so this is
#     really using the sausage inclusion zone...
#
setMethod('izGrid',
          signature(izObject = 'fullChainSawIZ', tract='Tract'),
function(izObject,
         tract,
         description = 'Full chainSaw (sausage) inclusion zone grid object',
         wholeIZ = TRUE,           #TRUE: grid the whole object; FALSE: just grid the IZ
         runQuiet = FALSE,
         ...
        )
{
#---------------------------------------------------------------------------
#
#   First, create an InclusionZoneGrid object that encompasses the entire
#   sausageIZ (fullChainSawIZ) inclusion zone and log. Then simply apply
#   the chain saw (chainSawIZ) method to each grid cell within the sausage
#   inclusion zone...
#
    izGridSausage = izGridConstruct(izObject=izObject, tract=tract, description=description,
                                    wholeIZ=wholeIZ, ...)
  
    plotRadius = izObject@radius
    downLog = izObject@downLog

    if(!runQuiet) 
      cat('\nThis can take some time...\ngrid cell: ')

#
#   now we need to assign all internal grid cells the correct value based
#   on a applying a chainSawIZ to each cell...
#
    grid = izGridSausage@grid
    numCells = ncell(grid)
    chiz = vector('list', numCells)
    names(chiz) = paste('izgCS',1:numCells,sep='.')
    mask = getValues(grid)              #vector valued
    #npua = length(izGridSausage@iz@puaEstimates)    
    npua = ncol(izGridSausage@data)    
    df = data.frame(matrix(0, nrow = numCells, ncol = npua)) #background grid defaults to zero
    #colnames(df) = names(izGridSausage@iz@puaEstimates)
    colnames(df) = names(izGridSausage@data)
    for(i in seq_len(numCells)) {
      if(!runQuiet)
        cat(i,', ',sep='')
      if(!is.na(mask[i])) {
        xy = xyFromCell(grid, i)[1,]                                       #make it a vector
        izCS = chainSawIZ(downLog, plotRadius=plotRadius, plotCenter = xy) #get point estimates
        izgCS = izGrid(izCS, tract, wholeIZ = FALSE)                       #one grid cell/point only!!!
        df[i, ] = izgCS@data                                               #nrow==1
        chiz[[i]] = izgCS
      }
      else
        chiz[[i]] = NA               #the rest of the background surface objects are missing 
    }
    if(!runQuiet)
      cat('\n')

#
#   create the object...
#
    griz = new('csFullInclusionZoneGrid',
               description = description,
               iz = izGridSausage@iz,
               chiz = chiz,
               grid = grid,
               data = df,
               bbox = izGridSausage@bbox
              )
    return(griz)
}   #izGridCSFull for the full chainSaw-sausage inclusion zone
)   #setMethod

               
  

  


          
#================================================================================
#  5. method for 'pointRelascopeIZ' and 'Tract' classes...
#
setMethod('izGrid',
          signature(izObject = 'pointRelascopeIZ', tract='Tract'),
function(izObject,
         tract,
         description = 'pointRelascopeIZ inclusion zone grid object',
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
}   #izGrid for'pointRelascopeIZ'
)   #setMethod

               
  

  


          
#================================================================================
#  6. method for 'perpendicularDistanceIZ' and 'Tract' classes...
#
setMethod('izGrid',
          signature(izObject = 'perpendicularDistanceIZ', tract='Tract'),
function(izObject,
         tract,
         description = 'perpendicularDistanceIZ inclusion zone grid object',
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
}   #izGrid for'perpendicularDistanceIZ'
)   #setMethod
               
  

  


          
#================================================================================
#  7. method for 'omnibusPDSIZ' and 'Tract' classes...  
#
setMethod('izGrid',
          signature(izObject = 'omnibusPDSIZ', tract='Tract'),
function(izObject,
         tract,
         description = 'omnibusPDSIZ inclusion zone grid object',
         wholeIZ = TRUE,           #TRUE: grid the whole object; FALSE: just grid the IZ
         runQuiet = TRUE,
         ...
        )
{
#---------------------------------------------------------------------------
#
#   first get the overall PDS inclusion zone grid object to work with...
#
    griz = izGridConstruct(izObject=izObject, tract=tract, description=description,
                           wholeIZ=wholeIZ, ...)
    
#
#   get the transformation matrices...
#
    dlog = izObject@downLog
    halfLen = dlog@logLen/2                                    #half length of log
    centerOffset = coordinates(dlog@location)
    logAngle = dlog@logAngle
    trMat = transfMatrix(logAngle, centerOffset)
    trMatInv = solve(trMat)

    factor = izObject@pds@factor
   
#
#   now we need to assign all internal grid cells the correct value based
#   on a applying omnibus estimates to each cell...
#
    grid = griz@grid
    numCells = ncell(grid)
    mask = getValues(grid)                                   #vector valued (either NA or zero)
    df = griz@data                                           #data frame of pua estimates
    for(i in seq_len(numCells)) {
      if(!runQuiet && identical(i%%10,0))
        cat(i,', ',sep='')
      if(!is.na(mask[i])) {
        xy = cbind(xyFromCell(grid, i), 1)                          #1x3 matrix in hc
        xy = xy %*% trMatInv                                        #transform back to canonical (names stripped)
        xLength = xy[1,1] + halfLen                                 #shift so butt is at zero
        #get log diameter at this grid point
        if(!is.null(dlog@solidType))                                #use default taper equation
          diam = .StemEnv$wbTaper(dlog@buttDiam, dlog@topDiam, dlog@logLen, nSegs=1,
                                  solidType=dlog@solidType, hgt=xLength)$diameter
        else                                                        #spline
          diam = spline(dlog@taper$length, dlog@taper$diameter, xout=xLength)$y
        diam = ifelse(identical(diam,0), 1e-04, diam)               #check for singularity point at tip
        #omnibus estimates...
        xsecArea = pi*diam^2/4
        denom = switch(izObject@pdsType,
                       volume = xsecArea, 
                       surfaceArea = pi*diam,
                       coverageArea = diam,
                       stop('Illegal pdsType in izGrid!')
                      )
        df[i, 'volume'] = factor*xsecArea/denom
        df[i, 'Density'] = factor/(denom * dlog@logLen)               #number of logs
        df[i, 'Length'] = factor/denom                                #length of logs
        df[i, 'surfaceArea'] = factor*pi*diam/denom                   #surface area of logs
        df[i, 'coverageArea'] = factor*diam/denom                     #coverage area of logs
        df[i, 'biomass'] = factor*dlog@conversions['volumeToWeight']*
                                  xsecArea/denom                      #woody biomass
        df[i, 'carbon'] = factor*dlog@conversions['volumeToWeight'] *
                                 dlog@conversions['weightToCarbon'] *
                                 xsecArea/denom                       #carbon content
      }
    }  #for
    if(!runQuiet)
      cat('\n')

    griz@data = df
    
    return(griz)
}   #izGrid for 'omnibusPDSIZ'
)   #setMethod
  






          
#================================================================================
#  8. method for 'distanceLimitedPDSIZ', 'omnibusDLPDSIZ' or 'hybridDLPDSIZ'
#     and 'Tract' classes...
#
#  Note that this will work for both of the above PDS classes "as is" since
#  the first is a superclass of the second/third...
#
#  The routine is set up in 3 main steps...
#  1. the entire stem is under PDS, i.e., not DL portion, so we just treat
#     it as normal PDS and use pdsFull to get the grid
#
#     otherwise construct the entire dlpds izGrid object so that we have a
#     bounding grid for the whole polygon--do not use the estimates here as 
#     they are not correct, just do this to get the grid object
#
#  2. if it exists, compute and expand the dlsPart of the izGrid and substitute
#     its values for the dummy values in step 1
#
#  3. if it exists, compute and expand the pdsPart of the izGrid and substitute
#     its values for the dummy values in step 1
#
#================================================================================
#
setMethod('izGrid',
          signature(izObject = 'distanceLimitedPDSIZ', tract='Tract'),
function(izObject,
         tract,
         description = 'a distance limited PDSIZ inclusion zone grid object',
         wholeIZ = TRUE,           #TRUE: grid the whole object; FALSE: just grid the IZ
         runQuiet = TRUE,
         ...
        )
{
#---------------------------------------------------------------------------
#   1. first see if the inclusion zone is completely PDS, if so, we just
#      dispatch to the correct constructor and we are done...
#
    if(is.null(izObject@dlsPart) && is.null(izObject@pdsPart)) {
      griz = izGrid(izObject=izObject@pdsFull, tract=tract, description=description,
                       wholeIZ=wholeIZ, runQuiet=runQuiet, ...)
      return(griz)
    }    

    densityAdjust = 0  #see end of routine for an explanation of this
  
#---------------------------------------------------------------------------
#
#   ---Now we have a situation where we are either all DLPDS, or a hybrid
#      of the two methods...
#
#      first get the overall izGrid object to work with--DO NOT use these 
#      puaEstimates as they are combined (if both DL and PDS in the same log)...
#
    griz = izGridConstruct(izObject=izObject, tract=tract, description=description,
                           wholeIZ=wholeIZ, ...)


#---------------------------------------------------------------------------
#
#   2. now, if there is a DL component (as a distanceLimitedMCIZ object),
#      then we must find its portion of the above raster object via overlay, and
#      finally assign the correct estimate values to the DL component grid cells...
#
    if(!is.null(izObject@dlsPart)) {
      if(!runQuiet)
        cat(' (dlsPart)')
#     get the variable inclusion zone, hence the use of izGrid here...      
      dlsgriz = izGrid(izObject=izObject@dlsPart, tract=tract, description=description,
                       wholeIZ=wholeIZ, ...)
      dlsgriz.vals = getValues(dlsgriz@grid)
      dlsFrom.idx = which(!is.na(dlsgriz.vals))            #index of DL values to replace with
      df = dlsgriz@data

#     overlay dls part with above full izGrid to get the right cells with same extent as griz...      
      izg = griz@grid
      grid = rasterize(perimeter(izObject@dlsPart), izg, mask=TRUE, silent=TRUE) #same extent as griz

      dls.vals = getValues(grid)   
      dlsTo.idx = which(!is.na(dls.vals))                  #index of the DL values to be replaced
      
#
#     replace the DLS portion of the grid for each [est]imate variable in the data frame;
#     note that the dlsFrom.idx and dlsTo.idx index sets will differ unless the entire log
#     is DL, but they will always have the same number of grid cells (length) to swap...
#
      if(length(dlsTo.idx > 0)) {
        data = griz@data
        #for(est in colnames(data)) 
          est = colnames(data)
          data[dlsTo.idx, est] = df[dlsFrom.idx, est]          

        griz@data = data

        densityAdjust = densityAdjust + ifelse(any(data$Density>0),1,0)
      }
    }  #DL component

    

#---------------------------------------------------------------------------
#
#   3. similarly, if there is a PDS component, then we must find its portion of the above
#      raster object via overlay, then assign the correct estimate values to the
#      PDS component grid cells; note that we are call izGrid below for omnibus, so this
#      will work for the variable surface in that method...
#
    if(!is.null(izObject@pdsPart)) {
      if(!runQuiet)
        cat(' (pdsPart)')
      if(is(izObject@pdsPart, 'omnibusPDSIZ')) {
#       get the variable omnibus inclusion zone from izGrid...      
        pdsgriz = izGrid(izObject=izObject@pdsPart, tract=tract, description=description,
                         wholeIZ=wholeIZ, ...)
        pdsgriz.vals = getValues(pdsgriz@grid)
        pdsFrom.idx = which(!is.na(pdsgriz.vals))            #index of PDS values to replace with
        df = pdsgriz@data
      } #omnibus only
      
#     overlay pds part with above full izGrid to get the right cells with same extent as griz...      
      izg = griz@grid
      grid = rasterize(perimeter(izObject@pdsPart), izg, mask=TRUE, silent=TRUE) #same extent as griz

      pds.vals = getValues(grid)   
      pdsTo.idx = which(!is.na(pds.vals))                  #index of the PDS values to be replaced
      
#
#     replace the PDS portion of the grid for each [est]imate variable in the data frame;
#     note that the pdsFrom.idx and pdsTo.idx index sets will differ unless the entire log
#     is PDS, but they will always have the same number of grid cells (length) to swap...
#
      if(length(pdsTo.idx > 0)) {
        data = griz@data
#       this test must be first as it is also a 'perpendicularDistanceIZ' object...
        if(is(izObject@pdsPart, 'omnibusPDSIZ')) {  
            est = colnames(data)
            data[pdsTo.idx, est] = df[pdsFrom.idx, est]                   #variable within each est
          }
        else {
          cnames = colnames(data)
          cnames = cnames[cnames!=.StemEnv$ppEstimates$depth]           #no depth in puaEstimates below
          for(est in cnames) 
            data[pdsTo.idx, est] = izObject@pdsPart@puaEstimates[[est]] #constant for each est
        }

        griz@data = data

        densityAdjust = densityAdjust + ifelse(any(data$Density>0),1,0)
      }
      
    } #PDS component

    
#   ------------------------
#   The estimates of density are both==1 if both methods are present (i.e., both components
#   are in the overall inclusion zone), so we have to factor the estimate of "2 logs" down by
#   half in this case...
#
#   However, it can happen that the inclusion zone for one or the other is just a sliver, and
#   even though present, will not contain any grid cell centers, and thus no estimate; this
#   case is taken care of by adjusting by densityAdjust, which will be either 1 or 2,
#   depending on whether any sample points (grid cell centers) actually fall within
#   the respective IZ component; thus, we get...
#      a. Both components present and contain sample points: densityAdjust=2
#      b. Both components present but one does not contain sample points: densityAdjust=1
#      c. Only DLS present (PDS-only is handled in section 1. above): densityAdjust=1
#      d. Both or one present, and the overall inclusion zone is so small it contains no grid
#         cell centers: densityAdjust=0 (gosh, let's hope this never happens!) 
#
    if(!(densityAdjust %in% c(1,2)))
      stop('densityAdjust = ',densityAdjust,' out of range in distanceLimitedPDSIZ (izGrid)!')
    griz@data[,'Density'] = griz@data[,'Density']/densityAdjust
    
    return(griz)
}   #izGrid for'distanceLimitedPDSIZ'
)   #setMethod
  
  











          
#================================================================================
#  9. method for 'distanceLimitedMCIZ' and 'Tract' classes...
#================================================================================
#
setMethod('izGrid',
          signature(izObject = 'distanceLimitedMCIZ', tract='Tract'),
function(izObject,
         tract,
         description = 'distanceLimitedMCIZ inclusion zone grid object',
         wholeIZ = TRUE,           #TRUE: grid the whole object; FALSE: just grid the IZ
         runQuiet = TRUE,
         ...
        )
{
#---------------------------------------------------------------------------
#
#   first get the overall izGrid object to work with...
#
    griz = izGridConstruct(izObject=izObject, tract=tract, description=description,
                           wholeIZ=wholeIZ, ...)
  
#
#   get the transformation matrices...
#
    dlog = izObject@downLog
    halfLen = dlog@logLen/2                                    #half length of log
    centerOffset = coordinates(dlog@location)
    logAngle = dlog@logAngle
    trMat = transfMatrix(logAngle, centerOffset)
    trMatInv = solve(trMat)

    puaBlowup = izObject@puaBlowup

#
#   now we need to assign all internal grid cells the correct value based
#   on a applying omnibus estimates to each cell...
#
    grid = griz@grid
    numCells = ncell(grid)
    mask = getValues(grid)                                          #vector valued (either NA or zero)
    df = griz@data                                                  #data frame of pua estimates
    for(i in seq_len(numCells)) {
      if(!runQuiet && identical(i%%10,0))
        cat(i,', ',sep='')
      if(!is.na(mask[i])) {
        xy = cbind(xyFromCell(grid, i), 1)                          #1x3 matrix in hc
        xy = xy %*% trMatInv                                        #transform back to canonical (names stripped)
        xLength = xy[1,1] + halfLen                                 #shift so butt is at zero
        #get log diameter at this grid point
        if(!is.null(dlog@solidType))                                #use default taper equation
          diam = .StemEnv$wbTaper(dlog@buttDiam, dlog@topDiam, dlog@logLen, nSegs=1,
                                  solidType=dlog@solidType, hgt=xLength)$diameter
        else                                                        #spline
          diam = spline(dlog@taper$length, dlog@taper$diameter, xout=xLength)$y
        diam = ifelse(identical(diam,0), 1e-04, diam)               #check for singularity point at tip
        xsecArea = pi*diam^2/4
        df[i, 'volume'] = xsecArea*puaBlowup 
        df[i, 'Density'] = puaBlowup / dlog@logLen                                  #number of logs
        df[i, 'Length'] = puaBlowup                                                 #length of logs
        df[i, 'surfaceArea'] = puaBlowup*pi*diam                                    #surface area of logs
        df[i, 'coverageArea'] = puaBlowup*diam                                      #coverage area of logs
        df[i, 'biomass'] = df[i, 'volume']*dlog@conversions['volumeToWeight']       #woody biomass
        df[i, 'carbon'] = df[i, 'biomass']*dlog@conversions['weightToCarbon']       #carbon content
      }
    }  #for loop

    if(!runQuiet)
      cat('\n')

    griz@data = df
    
    return(griz)
}   #izGrid for'distanceLimitedMCIZ'
)   #setMethod
  
               
  

  


          
#================================================================================
#  10. method for 'distanceLimitedIZ' and 'Tract' classes...
#================================================================================
#
setMethod('izGrid',
          signature(izObject = 'distanceLimitedIZ', tract='Tract'),
function(izObject,
         tract,
         description = 'distanceLimitedIZ inclusion zone grid object',
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
}   #izGrid for'distanceLimitedIZ'
)   #setMethod











#---------------------------------------------------------------------------
#
#   For standingTreeIZ subclass objects...
#
#   1. "circularPlotIZ"
#      -- "horizontalPointIZ" is a subclass of "circularPlotIZ" so that no
#         new method is required for it.
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
  


          
#================================================================================
#  1. method for 'circularPlotIZ' and 'Tract' classes...
#
setMethod('izGrid',
          signature(izObject = 'circularPlotIZ', tract='Tract'),
function(izObject,
         tract,
         description = 'circularPlotIZ inclusion zone grid object',
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
}   #izGrid for'circularPlotIZ'
)   #setMethod
