#---------------------------------------------------------------------------
#
#   This file has the collection of routines necessary to implement critical
#   height sampling variants of Lynch based on importance sampling and antithetic
#   variates for standing trees. The routines include...
#
#   1a. 'importanceCHSIZ' class structure
#   1a. 'antitheticICHSIZ' class structure
#   1b. 'pairedAICHSIZ' class structure
#
#   2.  'importanceCHSIZ' class generic 
#       'antitheticICHSIZ' class generic 
#       'pairedAICHSIZ' class generic 
#   2a. 'importanceCHSIZ' class constructor
#   2b. 'antitheticICHSIZ' class constructor
#   2c. 'pairedAICHSIZ' class constructor
#
#   3a. izGrid constructor for 'importanceCHSIZ' & 'Tract' classes
#   3b. izGrid constructor for 'antitheticICHSIZ' & 'Tract' classes
#   3c. izGrid constructor for 'pairedAICHSIZ' & 'Tract' classes
#
#   Note that 3b and 3c are just wrapper constructors that call the main
#   constructor (3a) through inheritance. It does all the work and decides what
#   kind of surface to calculate based on the class of the izObject passed;
#   see the routine for details.
#
#Author...									Date: 31-Jan-2013
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
# 1a. the importanceCHSIZ class is a direct descendant of 'criticalHeightIZ'...
#
#
setClass('importanceCHSIZ',
    representation(
                  ),
    contains = 'criticalHeightIZ',                     #a subclass of the 'criticalHeightIZ' class
    validity = function(object) {
                 return(TRUE)
               } #validity check
) #class importanceCHSIZ 


#=================================================================================================
#
# 1a. the antitheticICHSIZ class is a direct descendant of 'importanceCHSIZ'...
#
#
setClass('antitheticICHSIZ',
    representation(
                  ),
    contains = 'importanceCHSIZ',                     #a subclass of the 'importanceCHSIZ' class
    validity = function(object) {
                 return(TRUE)
               } #validity check
) #class antitheticICHSIZ 


#=================================================================================================
#
# 1b. the pairedAICHSIZ class is a direct descendant of 'importanceCHSIZ'...
#
#
setClass('pairedAICHSIZ',
    representation(
                  ),
    contains = 'importanceCHSIZ',                     #a subclass of the 'importanceCHSIZ' class
    validity = function(object) {
                 return(TRUE)
               } #validity check
) #class pairedAICHSIZ 










#================================================================================
# 2. generic definitions...
#================================================================================
#
if(!isGeneric("importanceCHSIZ")) 
  setGeneric('importanceCHSIZ',  
             function(standingTree, angleGauge, ...) standardGeneric('importanceCHSIZ'),
             signature = c('standingTree', 'angleGauge')
            )
 
if(!isGeneric("antitheticICHSIZ")) 
  setGeneric('antitheticICHSIZ',  
             function(standingTree, angleGauge, ...) standardGeneric('antitheticICHSIZ'),
             signature = c('standingTree', 'angleGauge')
            )
 
if(!isGeneric("pairedAICHSIZ")) 
  setGeneric('pairedAICHSIZ',  
             function(standingTree, angleGauge, ...) standardGeneric('pairedAICHSIZ'),
             signature = c('standingTree', 'angleGauge')
            )
 


       
#================================================================================
# 2a. and associated constructor method for class importanceCHSIZ...
#
setMethod('importanceCHSIZ',
          signature(standingTree = 'standingTree', angleGauge = 'angleGauge'),
function(standingTree,
         angleGauge,
         referenceHeight = .StemEnv$referenceCHSIZ,
         description = 'inclusion zone for importance CH sampling method',
         spID = paste('ichs',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         ...
        )
{
#------------------------------------------------------------------------------
#   it's the same as criticalHeightIZ, so just cast it...
#
    chsIZ = criticalHeightIZ(standingTree, angleGauge, referenceHeight, description,
                             spID, spUnits, ...)
    ichsIZ = as(chsIZ, 'importanceCHSIZ')               #cast

    return(ichsIZ)
}   #importanceCHSIZ constructor
)   #setMethod




#================================================================================
# 2b. and associated constructor method for class antitheticICHSIZ...
#
setMethod('antitheticICHSIZ',
          signature(standingTree = 'standingTree', angleGauge = 'angleGauge'),
function(standingTree,
         angleGauge,
         referenceHeight = .StemEnv$referenceCHSIZ,
         description = 'inclusion zone for antithetic ICH sampling method',
         spID = paste('aichs',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         ...
        )
{
#------------------------------------------------------------------------------
#   it's the same as importanceCHSIZ, so just cast it...
#
    ichsIZ = importanceCHSIZ(standingTree, angleGauge, referenceHeight, description,
                             spID, spUnits, ...)
    aichsIZ = as(ichsIZ, 'antitheticICHSIZ')               #cast

    return(aichsIZ)
}   #antitheticICHSIZ constructor
)   #setMethod



#================================================================================
# 2c. and associated constructor method for class pairedAICHSIZ...
#
setMethod('pairedAICHSIZ',
          signature(standingTree = 'standingTree', angleGauge = 'angleGauge'),
function(standingTree,
         angleGauge,
         referenceHeight = .StemEnv$referenceCHSIZ,
         description = 'inclusion zone for paired antithetic ICH sampling method',
         spID = paste('paichs',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         ...
        )
{
#------------------------------------------------------------------------------
#   it's the same as importanceCHSIZ, so just cast it...
#
    ichsIZ = importanceCHSIZ(standingTree, angleGauge, referenceHeight, description,
                             spID, spUnits, ...)
    paichsIZ = as(ichsIZ, 'pairedAICHSIZ')               #cast

    return(paichsIZ)
}   #pairedAICHSIZ constructor
)   #setMethod








#================================================================================
#================================================================================
#  3a. izGrid method for 'importanceCHSIZ' and 'Tract' classes...
#
#  Note in the individual point/grid cell loop below that it might appear that
#  we could have b>B when reference height=BH; but recall that the largest critical
#  distance will be that for dbh, so the smallest CH==BH, and the largest b==B; thus,
#  there is no need to check for a b>B in this case.
#
#  Note: uncomment the "##" lines below to play with the type of correction for the
#        butt (below dbh) section. The "uCrit" correction is what we finally decided
#        worked the best: 'cmc' will add randomness to the surface, but is unbiased
#        (stands for crude Monte Carlo) and 'cyl' just adds the volume of a clyinder,
#        which has a little bias because it is missing the volume in the flare area.
#
setMethod('izGrid',
          signature(izObject = 'importanceCHSIZ', tract='Tract'),
function(izObject,
         tract,
         description = 'importanceCHSIZ-based inclusion zone grid object',
         wholeIZ = TRUE,           #TRUE: grid the whole object; FALSE: just grid the IZ
         runQuiet = TRUE,
         ##buttType = c('uCrit','cmc','cyl'),    #butt vol correction when referenceHeight='dbh'
         ...
        )
{
#---------------------------------------------------------------------------
#
#   get the base constructor setup for this tree...
#
    griz = izGridConstruct(izObject=izObject, tract=tract, description=description,
                           wholeIZ=wholeIZ, ...)

#
#   determine the variant based on the type of object passed...
#
    if(is(izObject, 'pairedAICHSIZ')) 
      ichsType = 'paired'
    else if(is(izObject, 'antitheticICHSIZ'))
      ichsType = 'antithetic'
    else if((is(izObject, 'importanceCHSIZ')))
      ichsType = 'importance'
    else                                       #should never get here because of dispatching
      stop('izGrid for ICHS: You have selected an incompatible inclusion zone!')

    ##buttType = match.arg(buttType)

#    
#   just get everything handy from the object for ease of reading...
#
    standingTree = izObject@standingTree
    treeCenter = coordinates(standingTree@location)   #a 1x2 matrix of x,y
    taper = standingTree@taper
    height = standingTree@height
    buttDiam = standingTree@buttDiam
    dbh = standingTree@dbh
    topDiam = standingTree@topDiam
    solidType = standingTree@solidType

    #cross-sectional area (ba) conversion factor and height to dbh...
    baFactor =  ifelse(standingTree@units==.StemEnv$msrUnits$English, .StemEnv$baFactor['English'],
                       .StemEnv$baFactor['metric'])
    dbhHgt = ifelse(standingTree@units==.StemEnv$msrUnits$English, .StemEnv$dbhHgt['English'],
                       .StemEnv$dbhHgt['metric'])
    
    #reference cross-sectional area depends on reference height...
    referenceHeight = izObject@referenceHeight
    if(referenceHeight == 'butt')
      B = baFactor*buttDiam*buttDiam
    else                                     #dbh
      B = standingTree@ba
    
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
      if(is.na(mask[i]))                                            #skip background cells == NA
        next
      
      xy = xyFromCell(grid, i)                                      #1x2 grid cell center for sample point
      critDist = sqrt(sum((xy-treeCenter)^2))                       #critical distance to tree
      critDiameter = critDist/PRF                                   #critical diameter
      critCSA = critDiameter*critDiameter*baFactor                  #critical cross-sectional area b.c
      uCrit = critCSA/B                                             #u* == b.c/B

      #get a random diameter in base section for bias correction if sampling at dbh...
      if(referenceHeight == 'dbh') {
        ##if(buttType == 'cmc')
        ##  hFlare = runif(1, taper[1,'height'], dbhHgt)                #butt height should always be zero
        ##else if(buttType == 'uCrit') {
          hFlare = dbhHgt * uCrit
        ##}
        ##else if(buttType == 'cyl')
        ##  hFlare = dbhHgt
        dFlare = taperInterpolate(standingTree, 'diameter', hFlare)
        csaFlare = baFactor*dFlare*dFlare
      }        

      #now for the estimators: paired or ichs & aichs...
      if(ichsType == 'paired') {
        u.a = rep(0,2)
        b.a = u.a
        if(uCrit < 0.5)
          u.a[1] = uCrit + 0.5
        else
          u.a[1] = uCrit - 0.5
        u.a[2] = 1 - u.a[1]
        b.a[1] = B*(1-sqrt(1-u.a[1]))
        b.a[2] = B*(1-sqrt(u.a[1]))                                 #==B*(1-sqrt(1-u.a[2]))
        d.a = sqrt(b.a/baFactor)                                    #and their associated diameters
        # don't use the following ifelse, it evaluates the whole d.a vector in taperInterpolate and can
        # fail if one of the diameters meets the TRUE condition!...
        #bHeight = ifelse(d.a <= topDiam, height,                      #allow for trees with broken tops 
        #              taperInterpolate(standingTree, 'height', d.a))  #get height at diameter corresponding to 'b.a'
        bHeight = rep(0,2)
        for(j in 1:2) {                                             #this is safer, if slower
          if(d.a[j] <= topDiam)                                     #allow for trees with broken tops
            bHeight[j] = height  
          else                                                      #get height at diameter corresponding to 'b.a'
            bHeight[j] = taperInterpolate(standingTree, 'height', d.a[j])
        }
        if(referenceHeight == 'butt')
          volume = 0.25*baf*(bHeight[1]/sqrt(1-u.a[1]) + bHeight[2]/sqrt(u.a[1]) )
        if(referenceHeight == 'dbh')                                #correction for breast-height reference
          volume = 0.25*baf*((bHeight[1]-dbhHgt)/sqrt(1-u.a[1]) +
                             (bHeight[2]-dbhHgt)/sqrt(u.a[1]) + 4*csaFlare*dbhHgt/B)
      } #paired
      else {                                                        #ichs and antithetic ichs
        if(ichsType == 'importance')
          b = B*(1-sqrt(1-uCrit))                                   #upper stem CSA for sampled height (11)
        else if(ichsType == 'antithetic')
          b = B*(1-sqrt(uCrit))                                     #upper stem CSA for sampled height (14)
        d = sqrt(b/baFactor)                                        #and its associated diameter
        if(d < topDiam)                                             #allow for trees with broken tops
          bHeight = height 
        else                                                        #get height at diameter corresponding to 'b'
          bHeight = taperInterpolate(standingTree, 'height', d)

        if(ichsType ==  'importance') {
          if(referenceHeight == 'butt')
            volume = baf*bHeight/(2*sqrt(1-uCrit))         #estimator (13) 
          if(referenceHeight == 'dbh')                              
            volume = baf*((bHeight-dbhHgt)/(2*sqrt(1-uCrit)) + csaFlare*dbhHgt/B)
        }
        else if(ichsType ==  'antithetic') {
          if(referenceHeight == 'butt')
            volume = baf*bHeight/(2*sqrt(uCrit))           #estimator (16)
          if(referenceHeight == 'dbh')                              #correction for breast-height reference
            volume = baf*((bHeight-dbhHgt)/(2*sqrt(uCrit))  + csaFlare*dbhHgt/B)
        }
      } #importance & antithetic

      df[i, 'volume'] = volume
    }  #for loop

    if(!runQuiet)
      cat('\n')
    
    griz@data = df
     
    return(griz)
}   #izGrid for'importanceCHSIZ'
)   #setMethod




#================================================================================
#  3b. izGrid method for 'antitheticICHSIZ' and 'Tract' classes--just use
#      inherited method...
#
#
setMethod('izGrid',
          signature(izObject = 'antitheticICHSIZ', tract='Tract'),
function(izObject,
         tract,
         description = 'antitheticICHSIZ inclusion zone grid object',
         wholeIZ = TRUE,           #TRUE: grid the whole object; FALSE: just grid the IZ
         runQuiet = TRUE,
         ...
        )
{
#---------------------------------------------------------------------------
#
    griz = callNextMethod()
     
    return(griz)
}   #izGrid for'antitheticICHSIZ'
)   #setMethod
    



#================================================================================
#  3c. izGrid method for 'pairedAICHSIZ' and 'Tract' classes--just use
#      inherited method...
#
#
setMethod('izGrid',
          signature(izObject = 'pairedAICHSIZ', tract='Tract'),
function(izObject,
         tract,
         description = 'pairedAICHSIZ inclusion zone grid object',
         wholeIZ = TRUE,           #TRUE: grid the whole object; FALSE: just grid the IZ
         runQuiet = TRUE,
         ...
        )
{
#---------------------------------------------------------------------------
#
    griz = callNextMethod()
     
    return(griz)
}   #izGrid for'pairedAICHSIZ'
)   #setMethod
    
