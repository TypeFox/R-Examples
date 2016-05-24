#---------------------------------------------------------------------------
#
#   This file has the collection of routines necessary to implement critical
#   height sampling for standing trees. The routines include...
#
#   1. 'criticalHeightIZ' class structure
#   2. 'criticalHeightIZ' constructor
#   3. izGrid constructor for this class & 'Tract'
#   4. summary for 'criticalHeightIZ' class
#
#Author...									Date: 29-Jan-2013
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
# 1. the criticalHeightIZ class is a direct descendant of 'horizontalPointIZ'...
#
#
setClass('criticalHeightIZ',
    representation(referenceHeight = 'character'        #the reference height for the inclusion zone
                  ),
    contains = 'horizontalPointIZ',                     #a subclass of the 'horizontalPointIZ' class
    validity = function(object) {
      
                 if(!(object@referenceHeight %in% .StemEnv$referenceCHSIZ))
                   return(paste('referenceHeight must be one of:',
                                paste(.StemEnv$referenceCHSIZ,collapse=', ')))
                
                 return(TRUE)
               } #validity check
) #class criticalHeightIZ 
 






#================================================================================
# 2. generic definition...
#
if(!isGeneric("criticalHeightIZ")) 
  setGeneric('criticalHeightIZ',  
             function(standingTree, angleGauge, ...) standardGeneric('criticalHeightIZ'),
             signature = c('standingTree', 'angleGauge')
            )
 

       
#================================================================================
#  and associated constructor method for class criticalHeightIZ...
#
setMethod('criticalHeightIZ',
          signature(standingTree = 'standingTree', angleGauge = 'angleGauge'),
function(standingTree,
         angleGauge,
         referenceHeight = .StemEnv$referenceCHSIZ,
         description = 'inclusion zone for critical height sampling method',
         spID = paste('chs',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         ...
        )
{
#------------------------------------------------------------------------------
#
#   We just want to generate the inclusion zone using the inherited code from
#   the HPS class. But we want to be able to do this either at the base (completely
#   unbiased) or at dbh (which imparts a little bias as the "slabs" outside the
#   BHxDBH cylinder will not be accounted for, producing a small percent bias).
#   The simplest way to do this is to swap buttDiam into dbh in the standingTree
#   object, thereby fooling the HPS constructor into making an inclusion zone based
#   on it. Since we don't use anything else from the horizontalPointIZ object,
#   like the estimates, we can just do this and set everything else to NA below...
#
    referenceHeight = match.arg(referenceHeight)
    if(referenceHeight == 'butt') {  
      buttDiam = standingTree@buttDiam
      dbh = standingTree@dbh
      standingTree@dbh = buttDiam                       #fool it
    }
    hpsIZ = horizontalPointIZ(standingTree=standingTree, angleGauge=angleGauge,
                              description=description,
                              spID=spID, spUnits=spUnits, ...)
    chsIZ = as(hpsIZ, 'criticalHeightIZ')               #cast
    if(referenceHeight == 'butt')                       #reset dbh  
      chsIZ@standingTree@dbh = dbh

#
#   per unit area estimates: set all to NA...
#     the only applicable attribute for chs is volume, and it gets set based on
#     the critical height in the corresponding inclusionZoneGrid routine...
#
    chsIZ@puaBlowup = chsIZ@angleGauge@baf
    chsIZ@puaEstimates[c('volume','Density', 'basalArea','surfaceArea','biomass', 'carbon')] = NA_real_

    chsIZ@referenceHeight = referenceHeight
    
    return(chsIZ)
}   #criticalHeightIZ constructor
)   #setMethod







#================================================================================
#  3. izGrid method for 'criticalHeightIZ' and 'Tract' classes...
#
setMethod('izGrid',
          signature(izObject = 'criticalHeightIZ', tract='Tract'),
function(izObject,
         tract,
         description = 'criticalHeightIZ inclusion zone grid object',
         wholeIZ = TRUE,           #TRUE: grid the whole object; FALSE: just grid the IZ
         runQuiet = TRUE,
         ...
        )
{
#---------------------------------------------------------------------------
#
#   set up the grid object...
#
    griz = izGridConstruct(izObject=izObject, tract=tract, description=description,
                           wholeIZ=wholeIZ, ...)

    #just get everything handy from the object for ease of reading...
    standingTree = izObject@standingTree
    treeCenter = coordinates(standingTree@location)   #a 1x2 matrix of x,y
    taper = standingTree@taper
    height = standingTree@height
    buttDiam = standingTree@buttDiam
    topDiam = standingTree@topDiam
    solidType = standingTree@solidType
    
    angleGauge = izObject@angleGauge
    baf = angleGauge@baf
    PRF = angleGauge@PRF                              #remember PRF is in feet/foot (or m/m) of diameter

#
#   now we need to assign all internal grid cells the correct value based
#   on the critical height in each cell center/sample point...
#
    grid = griz@grid
    numCells = ncell(grid)
    mask = getValues(grid)                                          #vector valued (either NA or zero)
    df = griz@data                                                  #data frame of pua estimates
    
    #df$nHeap = rep(1,numCells) #**** could be used to identify joint inclusion zones easily
    for(i in seq_len(numCells)) {
      if(!runQuiet && identical(i%%10,0))
        cat(i,', ',sep='')
      if(is.na(mask[i]))                                          #skip background cells == NA
        next

      xy = xyFromCell(grid, i)                                    #1x2 grid cell center for sample point
      critDist = sqrt(sum((xy-treeCenter)^2))                     #critical distance to tree
      critDiameter = critDist/PRF                                 #critical diameter
      if(critDiameter < topDiam)                                  #allow for trees with broken tops
        critHeight = height 
      else
        critHeight = taperInterpolate(standingTree, 'height', critDiameter)
      
      df[i, 'volume'] = critHeight*baf                            #puaBlowup==baf here
    }  #for loop

    if(!runQuiet)
      cat('\n')

    griz@data = df
     
    return(griz)
}   #izGrid for'criticalHeightIZ'
)   #setMethod





#================================================================================
#  4. summary method for class "criticalHeightIZ"...
#
setMethod('summary',
          signature(object = 'criticalHeightIZ'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   add a little to 'horizontalPointIZ' method for 'criticalHeightIZ'...
#------------------------------------------------------------------------------
    callNextMethod()
   
    cat('\ncriticalHeightIZ...')
    cat('\n  Inclusion zone reference height is:', object@referenceHeight)

    cat('\n')
    
    return(invisible())
}   #summary for 'criticalHeightIZ'
) #setMethod
