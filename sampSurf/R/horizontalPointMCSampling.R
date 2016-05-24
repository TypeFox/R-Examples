#---------------------------------------------------------------------------
#
#   This file has the collection of routines for implementing various forms
#   of Monte Carlo subsampling within a horizontal point sampling areal frame.
#
#   The routines include...
#
#   A. Class definitions for inclusion zones...
#
#     1. 'horizontalPointCMCIZ' class structure for crude Monte Carlo sampling
#     2. 'horizontalPointISIZ' class structure for importance sampling
#     3. 'horizontalPointCVIZ' class structure for control variate sampling
#     4. 'horizontalPointMonteCarloSamplingIZ' class union of izGrid
#
#   B. Constructor generics and methods...
#
#     1. 'horizontalPointCMCIZ' generic and constructor method
#     2. 'horizontalPointISIZ' generic and constructor method
#     3. 'horizontalPointCVIZ' generic and constructor method
#
#   C. "InclusionZoneGrid" constructor method...
#
#     1. izGrid constructor which works on any of the A1-A3 classes via the A4
#        class union, plus a 'Tract' subclass object
#
#   Please note that we are using all the default plot, etc. methods for
#   a usual horizontalPointIZ object.
#
#***Important...
#
#   Each of the classes in part A contain both the "horizontalPointIZ" class
#   (which gives the hps functionality) and the "MonteCarloSamplingIZ" virtual
#   class (which gives slots for MC sampling methods). T
#
#   In the above, the idea in each of the IZ class definitions is to create
#   and store a "dummy" object of the MCS sampling class desired for the
#   inclusion zone. This makes all of the information like the number of
#   samples (n.s), etc. that is passed in the "..." list to the stem-based MCS
#   constructor available for future use in the corresponding izGrid routine,
#   which is where the actual MC application happens in each grid cell. The
#   call to the respective MC constructor within izGrid is done using do.call
#   so that we can stop any duplicate passing of arguments in the "..." list
#   (say from the sampSurf constructor where number of trees is passed, which
#   will require MC arguments like n.s to be in the "..." list) being passed
#   to more than the inclusion zone constructor; otherwise, R will throw an
#   error if duplicate arguments are passed.
#
#   So, the MCS object within each IZ class defined below, is just to hold the
#   desired MC argument options for future use in izGrid. And, it also gives
#   the class some meaning associated with the MC aspect of the inclusion zone.
#
#   The consequence of the above is that all IZ classes defined here are
#   direct descendents of "horizontalPointIZ" & "MonteCarloSamplingIZ", so there
#   is no MC hierarchy within the new HPS-MC IZ classes.
#
#Author...									Date: 14-May-2013
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   recall that the "MonteCarloSamplingIZ" class is defined in MonteCarloClass.R 
#
 





###########################################################################################
#
#   A. Class definitions...
#
###########################################################################################


#=================================================================================================
#
#  A.1 the horizontalPointCMCIZ class is a direct descendant of "MonteCarloSamplingIZ" and
#      'horizontalPointIZ'...
#
#  This is for crude Monte Carlo within an hps inclusion zone...
#
setClass('horizontalPointCMCIZ',
    contains = c('MonteCarloSamplingIZ', 'horizontalPointIZ')
) #class horizontalPointCMCIZ 


#=================================================================================================
#
#  A.2. the horizontalPointISIZ class is a direct descendant of "MonteCarloSamplingIZ" and
#       'horizontalPointIZ'...
#
#  This is for importance sampling within an hps inclusion zone...
#
setClass('horizontalPointISIZ',
    contains = c('MonteCarloSamplingIZ', 'horizontalPointIZ')
) #class horizontalPointISIZ 


#=================================================================================================
#
#  A.3. the horizontalPointCVIZ class is a direct descendant of "MonteCarloSamplingIZ" and
#       'horizontalPointIZ'...
#
#  This is for importance sampling within an hps inclusion zone...
#
setClass('horizontalPointCVIZ',
    contains = c('MonteCarloSamplingIZ', 'horizontalPointIZ')
) #class horizontalPointCVIZ 


#=================================================================================================
#
#  A.4. This class union allows us to use just on izGrid constructor for all of the methods
#       and makes things a little simpler, without repeated code...
#
setClassUnion('horizontalPointMonteCarloSamplingIZ',
              c('horizontalPointCMCIZ', 'horizontalPointISIZ', 'horizontalPointCVIZ')
             )












###########################################################################################
#
#   B. "InclusionZone" generics and constructor methods...
#
###########################################################################################

 

#================================================================================
#  B.1. crude Monte Carlo IZ generic definition...
#
if(!isGeneric("horizontalPointCMCIZ")) 
  setGeneric('horizontalPointCMCIZ',  
             function(standingTree, angleGauge, ...) standardGeneric('horizontalPointCMCIZ'),
             signature = c('standingTree', 'angleGauge')
            )
 
       
#================================================================================
#  and associated constructor method for class horizontalPointCMCIZ...
#
setMethod('horizontalPointCMCIZ',
          signature(standingTree = 'standingTree', angleGauge = 'angleGauge'),
function(standingTree,
         angleGauge,
         antithetic = FALSE,
         description = 'Inclusion zone for horizontal point CMC sampling method',
         spID = paste('hps.cmc',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         ...   #any arguments that should be passed on the the CMC constructor
        )
{
#------------------------------------------------------------------------------
#
#   We just want to generate the inclusion zone using the inherited code from
#   the HPS class.
       
    hpsIZ = horizontalPointIZ(standingTree=standingTree, angleGauge=angleGauge,
                              description=description,
                              spID=spID, spUnits=spUnits, ...)
    hpsCMCIZ = as(hpsIZ, 'horizontalPointCMCIZ')               #cast

#    
#   call constructor to do checks and to get userArgs if any--not used for
#   anything else at this point, the estimates are never used...
#
    mcsObj = crudeMonteCarlo(standingTree, ...)
    hpsCMCIZ@mcsObj = mcsObj
    hpsCMCIZ@antithetic = antithetic
    hpsCMCIZ@proxy = mcsObj@proxy
    
#
#   per unit area estimates: set all to NA...
#     the only applicable attribute for hpsCMC is volume, and it gets set 
#     in the corresponding inclusionZoneGrid routine...
#
    hpsCMCIZ@puaEstimates[c('volume','Density', 'basalArea','surfaceArea','biomass', 'carbon')] = NA_real_
    
    return(hpsCMCIZ)
}   #horizontalPointCMCIZ constructor
)   #setMethod





#================================================================================
#   B.2. importance sampling IZ generic definition...
#
if(!isGeneric("horizontalPointISIZ")) 
  setGeneric('horizontalPointISIZ',  
             function(standingTree, angleGauge, ...) standardGeneric('horizontalPointISIZ'),
             signature = c('standingTree', 'angleGauge')
            )
 

       
#================================================================================
#  and associated constructor method for class horizontalPointISIZ...
#
setMethod('horizontalPointISIZ',
          signature(standingTree = 'standingTree', angleGauge = 'angleGauge'),
function(standingTree,
         angleGauge,
         antithetic = FALSE,
         description = 'Inclusion zone for horizontal point importance sampling method',
         spID = paste('hps.is',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         ...
        )
{
#------------------------------------------------------------------------------
#
#   We just want to generate the inclusion zone using the inherited code from
#   the HPS class.
       
    hpsIZ = horizontalPointIZ(standingTree=standingTree, angleGauge=angleGauge,
                              description=description,
                              spID=spID, spUnits=spUnits, ...)
    hpsISIZ = as(hpsIZ, 'horizontalPointISIZ')               #cast

#    
#   call constructor to do checks and to get userArgs if any--not used for
#   anything else at this point, the estimates are never used...
#
    mcsObj = importanceSampling(standingTree, ...)           #proxy can be in "..."s
    hpsISIZ@mcsObj = mcsObj
    hpsISIZ@antithetic = antithetic
    hpsISIZ@proxy = mcsObj@proxy
    
#
#   per unit area estimates: set all to NA...
#     the only applicable attribute for hpsIS is volume, and it gets set 
#     in the corresponding inclusionZoneGrid routine...
#
    hpsISIZ@puaEstimates[c('volume','Density', 'basalArea','surfaceArea','biomass', 'carbon')] = NA_real_
    
    return(hpsISIZ)
}   #horizontalPointISIZ constructor
)   #setMethod





#================================================================================
#   B.3. control variate sampling IZ generic definition...
#
if(!isGeneric("horizontalPointCVIZ")) 
  setGeneric('horizontalPointCVIZ',  
             function(standingTree, angleGauge, ...) standardGeneric('horizontalPointCVIZ'),
             signature = c('standingTree', 'angleGauge')
            )
 

       
#================================================================================
#  and associated constructor method for class horizontalPointCVIZ...
#
setMethod('horizontalPointCVIZ',
          signature(standingTree = 'standingTree', angleGauge = 'angleGauge'),
function(standingTree,
         angleGauge,
         antithetic = FALSE,
         description = 'Inclusion zone: horizontal point with control variate sampling',
         spID = paste('hps.cv',.StemEnv$randomID(),sep=':'),
         spUnits = CRS(projargs=as.character(NA)),
         ...
        )
{
#------------------------------------------------------------------------------
#
#   We just want to generate the inclusion zone using the inherited code from
#   the HPS class.
       
    hpsIZ = horizontalPointIZ(standingTree=standingTree, angleGauge=angleGauge,
                              description=description,
                              spID=spID, spUnits=spUnits, ...)
    hpsCVIZ = as(hpsIZ, 'horizontalPointCVIZ')               #cast

#    
#   call constructor to do checks and to get userArgs if any--not used for
#   anything else at this point, the estimates are never used...
#
    mcsObj = controlVariate(standingTree, ...)           #proxy can be in "..."s
    hpsCVIZ@mcsObj = mcsObj
    hpsCVIZ@antithetic = antithetic
    hpsCVIZ@proxy = mcsObj@proxy
    
#
#   per unit area estimates: set all to NA...
#     the only applicable attribute for hpsIS is volume, and it gets set 
#     in the corresponding inclusionZoneGrid routine...
#
    hpsCVIZ@puaEstimates[c('volume','Density', 'basalArea','surfaceArea','biomass', 'carbon')] = NA_real_
    
    return(hpsCVIZ)
}   #horizontalPointCVIZ constructor
)   #setMethod













###########################################################################################
#
#   C. "InclusionZoneGrid" constructor for all MCS methods...
#
###########################################################################################


#================================================================================
#  2. izGrid method for 'horizontalPoint???IZ' and 'Tract' classes...
#
setMethod('izGrid',
          signature(izObject = 'horizontalPointMonteCarloSamplingIZ', tract='Tract'),
function(izObject,
         tract,
         description = 'horizontalPoint w/ MC subsampling inclusion zone grid object',
         wholeIZ = TRUE,           #TRUE: grid the whole object; FALSE: just grid the IZ
         ...
        )
{
#---------------------------------------------------------------------------
#
#   Note that we must visit each grid cell individually and get an estimate
#   based on the n.s subsamples; it would be nice if we could draw all the
#   the subsamples in one shot from the tree, but then we'd still need to
#   get the individual grid cell estimates based on these, so what is done
#   below is reasonable, even though a bit slow with the class structure...
#
#   Please see "horizontalPointCMCIZ" class definition, all of the MC arguments
#   are defined in the CMC object stored there, so if any are passed to
#   izGrid in the "..." call here, they are ignored (see also comments at
#   the intro to this file).
#
#---------------------------------------------------------------------------
#
#   set up the grid object...
#
    griz = izGridConstruct(izObject=izObject, tract=tract, description=description,
                           wholeIZ=wholeIZ, ...)

    #just get everything handy from the object for ease of reading...
    standingTree = izObject@standingTree    
    puaBlowup = izObject@puaBlowup         #the normal HPS blowup factor
    
    antithetic = izObject@antithetic

#
#   initialize the random number stream once, for all trees; note that is should not be
#   passed in args below, or the stream will be re-initialized (if !NA) for each tree,
#   which results in zero variance...
#
    initRandomSeed(izObject@mcsObj@startSeed)

#    
#   set up the argument list for the constructor in the do.call below...
#
    args = list(object = standingTree,
                n.s = izObject@mcsObj@n.s,
                u.s = NA,
                segBnds = izObject@mcsObj@segBnds,
                startSeed = NA,                              #all trees go from current stream above
                alphaLevel = izObject@mcsObj@alphaLevel
               )
    userArgs = izObject@mcsObj@userArgs                      #not for cmc normally, but it can't hurt
    if(length(userArgs)>0)                                   #any extras to be passed?
      args = c(args, userArgs)                               #include them if so

#
#   set up the function for the call, and augment the proxy argument if need be...
#
    if(is(izObject, 'horizontalPointCMCIZ'))
      mcs.f = 'crudeMonteCarlo'
    else {
      if(is(izObject, 'horizontalPointISIZ')) 
        mcs.f = 'importanceSampling'
      else if(is(izObject, 'horizontalPointCVIZ')) 
        mcs.f = 'controlVariate'
      else
        stop('No known MCS method for izGrid: this should not happen!')
      args = c(args, list(proxy = izObject@proxy) )
    }

#
#   now we need to assign all internal grid cells the correct value based
#   on a random volume estimate in each cell center/sample point...
#
    grid = griz@grid                                       #raster grid object
    mask = getValues(grid)                                 #vector (values either NA or zero)
    df = griz@data                                         #data frame of pua estimates
    df$variance = rep(0, nrow(df))                         #we will save with-tree variance too

    idx = which(!is.na(mask))                              #internal IZ grid cell value indices vector
    if(length(idx) > 0) {
      for(i in idx) {
        mcsObj = do.call(mcs.f, args)                      #call the correct method
        if(antithetic)                                     #do we want an antithetic estimate?
          mcsObj = antitheticSampling(mcsObj)
        df[i, 'volume'] = mcsObj@volEst * puaBlowup
        df[i, 'variance'] = mcsObj@volVar                        #need to be weighted???? <<<*****?????????
      }
      griz@data = df
    }
         
    return(griz)
}   #izGrid for'horizontalPoint???IZ'
)   #setMethod
