checkStemDimensions = function(stem,
                               tolerancePercent = 1, #percentage of given attribute
                               ...
                              )
{
#---------------------------------------------------------------------------
#
#   One place where there may be too much freedom in sampSurf is the ability
#   to have volumes, etc., that are different from those that are given by
#   the taper data; these can be specified in the constructor. The problem
#   comes with Monte Carlo methods that use the taper data for intermediate
#   dimensions that are then used to estimate the attribute of interest. So
#   if, e.g., volume is different from the integrated taper volume, a bias
#   could result, depending on the magnitude of the difference.
#
#   This routine is used to check for such inconsistencies and warn the user
#   if they are found. It will be called from within the constructors for
#   "Stem" objects.
#
#   Arguments...
#     stem = A "Stem" subclass object
#     tolerancePercent = the tolerance limit in percent difference; the default
#                        is set small because it really is not known how little
#                        variation can be tolerated in Monte Carlo methods that
#                        depend on reconciling the taper and actual attribute
#                        values.
#
#   Returns...
#     A list invisibly with...
#
#      stemID = the stem ID
#      tolerancePercent = the tolerance as described above
#      volume = 1. volume from the slot, 2. volume from either the taper or
#               Smalian's, depending on solidType, 3. splined volume from taper,
#               4. the percent difference
#      surfaceArea = 1. surface area slot, 2. splined or taper SA depending on
#                    solidType, 3. percent difference
#      coverageArea = as above but for CA
#      biomass = 1. biomass slot, 2. conversion from (2) volume, 3. percent diff
#      carbon = as for biomass
#
#   Note that we could also add a list to the return above with the actual
#   warning messages in it.
#
#Author...									Date: 16-Jan-2014
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   a quick check...
#
    if(!is(stem, 'Stem'))
       stop('Please pass a "Stem" class object!')

#
#   attribute & slot names & messages...
#
    volName = 'volume'
    saName = 'surface area'
    saSlot = 'surfaceArea'
    caName = 'coverage area'
    caSlot = 'coverageArea'
    bmsName = 'biomass'
    bmsSlot = bmsName
    cName = 'carbon'
    cSlot = cName

    stemID = getID(stem)
    prose = c('A significant difference in stem',
              'was found when comparing the contents of the',
              '\nslot to the',
              paste('calculated from the taper slot data frame in ',stemID,'.',sep='')
             )
    problemMessage = paste('\nThis can cause a bias in the sampling surface simulation, which could potentially be',
                           'serious when using Monte Carlo methods.')

#
#   stem attributes...
#
    solidType = stem@solidType
    taper = stem@taper
    topDiam = stem@topDiam
    buttDiam = stem@buttDiam
    surfaceArea = stem@surfaceArea
    biomass = stem@biomass
    carbon = stem@carbon
    volumeToWeight = as.numeric(stem@conversions['volumeToWeight']) #strip off the name
    weightToCarbon = as.numeric(stem@conversions['weightToCarbon']) #strip off the name
    
    if(is(stem, 'downLog')) {
      isLog = TRUE
      coverageArea = stem@coverageArea
      height = stem@logLen
      treeVol = stem@logVol
      volSlot = 'logVol'
    }
    else {
      isLog = FALSE
      coverageArea = NA_real_
      height = stem@height
      treeVol = stem@treeVol
      volSlot = 'treeVol'
    }     

#
#   volume check...
#
    if(is.null(solidType)) {                                     #user-defined taper
      treeVol2 = .StemEnv$SmalianVolume(taper, isLog)$logVol     #logVol is treeVol here
      splineVol = .StemEnv$splineVolume(taper, 0, height, isLog) #via spline vs Smalian's
    }
    else {                                                   #from default taper equation
      treeVol2 = .StemEnv$wbVolume(buttDiam, topDiam, height, solidType)
      splineVol = NA_real_
    }

    volDiff = abs(treeVol-treeVol2)/treeVol2*100     #relative to the taper version
    if(volDiff > tolerancePercent) {
      diffMess = paste('(',format(volDiff,dig=4),'%)',sep='')
      mess = paste(prose[1], volName, diffMess, prose[2], volSlot, prose[3], volName,
                  prose[4], problemMessage)
      warning(mess)
    }
    volumes = c(stemVol=treeVol, taperVol=treeVol2, spline=splineVol, volDiffPct=volDiff)


#
#   surface area check...
#
    if(is.null(solidType))                        #user-defined taper
      surfaceArea2 = .StemEnv$splineSurfaceArea(taper, lenBot=0, lenTop=height, isLog=isLog)   #spline function
    else                                          #default taper equation
      surfaceArea2 = .StemEnv$wbSurfaceArea(buttDiam, topDiam, height, solidType, lenBot=0, lenTop=height)

    saDiff = abs(surfaceArea-surfaceArea2)/surfaceArea2*100     #relative to the taper version
    if(saDiff > tolerancePercent) {
      diffMess = paste('(',format(saDiff,dig=4),'%)',sep='')
      mess = paste(prose[1], saName, diffMess, prose[2], saSlot, prose[3], saName,
                   prose[4], problemMessage)
      warning(mess)
    }
    surfaceAreas = c(surfaceArea=surfaceArea, taperSA=surfaceArea2, saDiffPct=saDiff)

#
#   coverage area only if a donwLog...
#
    if(isLog) {
      if(is.null(solidType))                        #user-defined taper
        coverageArea2 = .StemEnv$splineCoverageArea(taper, lenBot=0, lenTop=height)   #spline function
      else                                          #default taper equation
        coverageArea2 = .StemEnv$wbCoverageArea(buttDiam, topDiam, height, solidType, lenBot=0, lenTop=height)

      caDiff = abs(coverageArea-coverageArea2)/coverageArea2*100     #relative to the taper version
      if(caDiff > tolerancePercent) {
        diffMess = paste('(',format(caDiff,dig=4),'%)',sep='')
        mess = paste(prose[1], caName, diffMess, prose[2], caSlot, prose[3], caName,
                     prose[4], problemMessage)
        warning(mess)
      }
      coverageAreas = c(coverageArea=coverageArea, taperCA=coverageArea2, caDiffPct=caDiff)
    }
    else
      coverageAreas = c(coverageArea=NA_real_, taperCA=NA_real_, caDiffPct=NA_real_)


#
#   biomass--note that since it is calculated from treeVol/logVol, we can again
#            compare it to the taper data...
#
    if(!is.na(biomass)) {
      biomass2 = treeVol2*volumeToWeight
      bmsDiff = abs(biomass-biomass2)/biomass2*100     #relative to the taper version
      if(bmsDiff > tolerancePercent) {
        diffMess = paste('(',format(bmsDiff,dig=4),'%)',sep='')
        mess = paste(prose[1], bmsName, diffMess, prose[2], bmsSlot, prose[3], bmsName,
                    prose[4], problemMessage)
        warning(mess)
      }
      Biomass = c(biomass=biomass, taperBms=biomass2, bmsDiffPct=bmsDiff)
    }
    else
      Biomass = c(biomass=NA_real_, taperBms=NA_real_, bmsDiffPct=NA_real_)


#
#   carbon--note that since it is calculated from biomass via treeVol/logVol, we can again
#            compare it to the taper data...
#
    if(!is.na(carbon)) {
      carbon2 = biomass2*weightToCarbon
      cDiff = abs(carbon-carbon2)/carbon2*100     #relative to the taper version
      if(cDiff > tolerancePercent) {
        diffMess = paste('(',format(cDiff,dig=4),'%)',sep='')
        mess = paste(prose[1], cName, diffMess, prose[2], cSlot, prose[3], cName,
                    prose[4], problemMessage)
        warning(mess)
      }
      Carbon = c(carbon=carbon, taperC=carbon2, cDiffPct=cDiff)
    }
    else
      Carbon = c(carbon=NA_real_, taperC=NA_real_, cDiffPct=NA_real_)
    
    
#
#   send back a list with the relevent data...
#    
    return(invisible(list(stemID = stemID,
                          tolerancePercent = tolerancePercent,
                          volume = volumes,
                          surfaceArea = surfaceAreas,
                          coverageArea = coverageAreas,
                          biomass = Biomass,
                          carbon = Carbon
                         )
                    )
          )
}   #checkStemDimensions
