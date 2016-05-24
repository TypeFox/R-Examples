#---------------------------------------------------------------------------
#
#   Hidden global environment for class 'Stem' stuff.
#
#   Note that this now holds constants/parameters, etc. for other classes
#   within sampSurf as well. JHG 16-Dec-2010.
#
#   Note that this environment and its bindings are locked so that no one
#   can change it when it is attached. If it is not locked, then it is quite
#   simple to mess with it by changing any of the quantities inside. However,
#   locking can not be done at the end of this file, it must be done inside
#   the .onLoad() function as the package is loaded. 23-Mar-2011, JHG.
#
#Author...									Date: 6-Aug-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------





#---------------------------------------------------------------------------
#
#   Specifically for 'Stem' class and subclass objects and pertinent methods.
#
#   Creates a hidden environment for any global constants that might be
#   shared among different routines. Just run this to re-create. Add any
#   new assignments below.
#   
#---------------------------------------------------------------------------
#
#   create hidden environment to store things in within the package...
#
.StemEnv = new.env()



#
#   units of measure types, conversions and constants...
#
.StemEnv$msrUnits = list(English='English', metric='metric')
.StemEnv$cm2m = 1/100
.StemEnv$m2cm = 100
.StemEnv$in2ft = 1/12
.StemEnv$ft2in = 12

.StemEnv$dbhHgt = c(English = 4.5, metric = 1.3716)  #in English and metric

.StemEnv$smpHectare = 10000
.StemEnv$sfpAcre = 43560

#used in, e.g., mirage, walkthrough...
.StemEnv$cardinal = c('north', 'south', 'east', 'west')
.StemEnv$Cardinal = c('North', 'South', 'East', 'West')




#the following is ba in ft^2 or m^2 conversion, NOT the BAF for prism cruising!...
.StemEnv$baFactor = c( English = pi/4, metric = pi/4)     #dbh in feet or meters
#.StemEnv$baFactor = c( English = pi/(4*144), metric = pi/(4*10000) ) #dbh in inches or cm

.StemEnv$angleGaugeMaxDegrees = 6.5   #maximum angle for angleGauge objects in ArealSampling Class

#horizontal line sampling (hls) constants...
.StemEnv$HLSSegment = c(English = 66, metric = 20) #hls segment base length in ft or m






#
# per unit area names for the list object slot in the InclusionZone class or subclasses;
# I have perhaps made this too difficult, but one can assign the slots based on the names
# of the list, while the values could someday change; e.g.
# eval(parse(text=paste('list(',.StemEnv$puaEstimates$volume,'=3)')))
# will do it...
#
.StemEnv$puaEstimates = list(volume = 'volume',                    #cubic volume in correct units
                             #bfVolume = 'bfVolume',                #board foot volume--for trees!
                             Density = 'Density',                  #number of logs or trees
                             Length = 'Length',                    #total length of logs
                             surfaceArea = 'surfaceArea',          #surface area
                             coverageArea = 'coverageArea',        #projected coverage area--down logs
                             biomass = 'biomass',                  #could be green or dry
                             carbon = 'carbon',                    #carbon content
                             basalArea = 'basalArea',              #basal area for standing trees
                             variance = 'variance'                 #required for MC methods 
                            )
#not pua or estimates, other sampling related attributes: pp stands for per point...
.StemEnv$ppEstimates = list(depth = 'depth')                       #sample or overlap depth per point
                        
#
#  these are the valid estimates from the above list for each main subclass of Stem object...
#
.StemEnv$validEstimates = list(downLogs = c('volume','Density','Length','surfaceArea','coverageArea',
                                            'biomass','carbon'),
                               standingTrees = c('volume','Density','surfaceArea','basalArea',
                                                 'biomass','carbon')
                              )
  
#.StemEnv$puaNames = list(English =
#                          list(volume = 'volume in cubic feet per acre',
#                               bfVolume = 'board foot volume per acre',
#                               logDensity = 'number of logs per acre'
#                              ),
#                          metric =
#                           list(volume = 'volume in cubic meters per hectare',
#                                bfVolume = NULL,
#                                logDensity = 'number of logs per hectare'
#                               )
#                         ) #puaNames


#
#   stuff for PDS...
#
.StemEnv$pdsTypes = c('volume', 'surfaceArea', 'coverageArea')  #main PPS variable


#
#   some plausible species codes/names--note that they can be any character string...
#
.StemEnv$species = c('wp','rm','sm','hemlock','Picea glauca','shagbark hickory','BW')


#
#   possible log lie angles & other angular constants...
#
.StemEnv$logAngles = c(0, 2*pi)
.StemEnv$deg2rad = pi/180
.StemEnv$rad2deg = 1/.StemEnv$deg2rad


#
#   some useful defaults for plotting consistently in different classes...
#
.StemEnv$alphaTrans = 0.5


#log colors, etc.
.StemEnv$logBorderColor = 'brown4'  #perimeter color of downLog objects
.StemEnv$logColor = transparentColorBase('sandybrown', .StemEnv$alphaTrans)  #internal shade for down logs
.StemEnv$logAttributeColor = transparentColorBase('tan3',              #needle & center color
                                    ifelse(1.5*.StemEnv$alphaTrans>1, 1, 1.5*.StemEnv$alphaTrans))


#standing tree colors, etc. (could be for dbh circle or standing tree as needed)...
.StemEnv$treeBorderColor = 'brown4'  #perimeter color for dbh circle 
.StemEnv$treeColor = transparentColorBase('sandybrown', .StemEnv$alphaTrans)  #internal shade for dbh circle
.StemEnv$treeAttributeColor = transparentColorBase('tan3',                    #center color
                                    ifelse(1.5*.StemEnv$alphaTrans>1, 1, 1.5*.StemEnv$alphaTrans))

#inclusion zones or plot-related...
.StemEnv$izBorderColor = transparentColorBase('slategray',                      #plot perimeter color
                                    ifelse(1.5*.StemEnv$alphaTrans>1, 1, 1.5*.StemEnv$alphaTrans)) 
.StemEnv$izColor = transparentColorBase('gray95', .StemEnv$alphaTrans)  #interior color for sample plots
.StemEnv$izCenterColor = transparentColorBase('slategray', .StemEnv$alphaTrans)  #center point color


#zero=lightest
.StemEnv$blue.colors = colorRampPalette(c('lightsteelblue1','steelblue','steelblue4'))
.StemEnv$gray.colors = colorRampPalette(c('grey90','grey50'))
#blue.colors = colorRampPalette(c('lightsteelblue1','steelblue','steelblue4'))
#zero=darkest...
#.StemEnv$blue.colors = colorRampPalette(c('steelblue4','steelblue','lightsteelblue'))
#blue.colors = colorRampPalette(c('steelblue4','steelblue','lightsteelblue'))

.StemEnv$gridLineColor = transparentColorBase('slategray',.StemEnv$alphaTrans)
.StemEnv$gridCenterColor = transparentColorBase('firebrick4',.StemEnv$alphaTrans)


#
#   critical height and importance variants of CHS...
#
.StemEnv$referenceCHSIZ = c('butt', 'dbh') #reference height for critical height inclusion zone

#
#   Monte Carlo Sampling (importance, etc.)...
#
.StemEnv$isHeightColor = transparentColorBase('tan4',                    #height locations
                                    ifelse(1.5*.StemEnv$alphaTrans>1, 1, 1.5*.StemEnv$alphaTrans))


#
#   sampleLogs & sampleTrees stuff...
#
.StemEnv$sampleLogsNames = c('species', 'logLen', 'buttDiam', 'topDiam', 'solidType',
                             'x', 'y', 'logAngle', 'logAngle.D')
.StemEnv$sampleTreesNames = c('species', 'height', 'dbh', 'topDiam', 'solidType',
                             'x', 'y')


#
#   taper/volume stuff...
#
.StemEnv$solidTypes = c(1,10)  #range for valid solidType in log taper/volume




#================================================================================
#
#  default taper function...
#
#  Note that the taper function must have diameters in the same units as length,
#  so if length is in meters, diameters must be as well...
#
#  Arguments...
#    botDiam = the diameter at the large end in height/length units
#    topDiam = the diameter at the small end in height/length units
#    logLen = the length of the log
#    nSegs = the number of segments desired
#    solidType = between 1 and 10 is legal
#    hgt = NULL to calculate sections here; otherwise a vector of hgt/length
#          section information
#    isLog = TRUE: a down log, so use "length" in taper; FALSE: a standing tree
#            so use "hgt" in taper
#
#   Note: little error checking has been added yet!!!!
#
wbTaper = function(botDiam, topDiam, logLen, nSegs=20, solidType, hgt=NULL, isLog=TRUE) {
    if(nSegs < 1)
      stop('Must have positive number of log segments for taper!')
    nSegs = nSegs + 1               #becomes the number of diameters required for taper
    if(is.null(solidType) || solidType < .StemEnv$solidTypes[1] || solidType > .StemEnv$solidTypes[2])
      stop('solidType=',solidType,' out of range, must be in: [',solidTypes[1],',',solidTypes[2],']')
    r = solidType
    
    if(is.null(hgt))
      hgt = seq(0, logLen, length.out=nSegs)
    diameter = topDiam + (botDiam - topDiam) * ((logLen - hgt)/logLen)^(2/r)
    taper = data.frame(diameter=diameter, hgt=hgt)
    if(isLog)
      colnames(taper) = c('diameter','length')
    else
      colnames(taper) =  c('diameter','height')
    return(taper)
}   #wbTaper
assign('wbTaper', wbTaper, envir=.StemEnv)            #move to .StemEnv
environment(.StemEnv$wbTaper) = .StemEnv              #assign its environment
rm(wbTaper)                                           #and remove from .GlobalEnv



#================================================================================
#
#  default volume function based on default taper function...
#
#  k is the conversion factor that takes diameter to radius and puts it into the
#  same units as length. But diameters should be in length units for downLogs,
#  so k just represents taking squared diameter to squared radius...
#
#  To get actual bolt volume of a segment somewhere on the stem, call this twice,
#  first with the shorter length, then the longer length (both defining the bolt)
#  and get the volume by subtraction--see segmentVolumes() and chainsawSliver
#  for examples--but now one should simply use segmentVolumes 25-Apr-2013.
#
wbVolume = function(botDiam, topDiam, logLen, solidType, boltLen=NULL) {
    if(is.null(solidType) || solidType < .StemEnv$solidTypes[1] || solidType > .StemEnv$solidTypes[2])
      stop('solidType=',solidType,' out of range, must be in: (',
           .StemEnv$solidTypes[1],',',.StemEnv$solidTypes[2],')')
    
    r = solidType
    k = 1/4                              #diameter to radius; diam to length units conversion==1
    if(is.null(boltLen))
      h = logLen                         #any height/length is possible, default is for total volume
    else
      h = boltLen                        #some intermediate volume
                                         
    logVol = pi*k*topDiam^2*h +
             pi*k*logLen*(botDiam - topDiam)^2*r/(4+r)*(1-(1-h/logLen)^((4+r)/r)) +
             2*pi*k*logLen*topDiam*(botDiam - topDiam)*r/(2+r)*(1-(1-h/logLen)^((2+r)/r))
    #logVol = pi*k*logLen*((buttD-topD)^2*(r/(r+4)) + topD*(topD + 2*(buttD-topD)*(r/(r+2))))
    return(logVol)
}   #wbVolume
assign('wbVolume', wbVolume, envir=.StemEnv)               #move to .StemEnv
environment(.StemEnv$wbVolume) = .StemEnv                  #assign its environment
rm(wbVolume)                                               #and remove from .GlobalEnv




#================================================================================
#
#  Smalian's volume function for passed taper...
#
#  k is the conversion factor that takes diameter to radius and puts it into the
#  same units as length. But diameters should be in length units for downLogs,
#  so k just represents taking squared diameter to squared radius...
#
SmalianVolume = function(taper, isLog=TRUE) {
    k = 1/4                              #diameter to radius; diam to length units conversion==1
    nSegs = nrow(taper) - 1
    if(nSegs < 1)
      stop("Must have positive number of log segments for Smalian's!")
    if(isLog)
      hgtName = 'length'
    else
      hgtName = 'height'
    vol = matrix(NA, nrow=nSegs, ncol=1)
    diam = taper[,'diameter']
    csArea = diam^2
    length = taper[,hgtName]
    for(i in 1:nSegs) {
      sectLen = length[i+1] - length[i]
      if(isTRUE(all.equal(diam[i+1],0.0)))
        vol[i,1] = pi*k*csArea[i+1]*sectLen/3                 #cone for tip
      else
        vol[i,1] = pi*k*(csArea[i] + csArea[i+1])*sectLen/2   #Smalian's
    }
    sumVol = colSums(vol)
    return(list(boltVol = vol, logVol = sumVol))
}   #SmalianVolume
assign('SmalianVolume', SmalianVolume, envir=.StemEnv)     #move to .StemEnv
environment(.StemEnv$SmalianVolume) = .StemEnv             #assign its environment
rm(SmalianVolume)                                          #and remove from .GlobalEnv



#================================================================================
#
#  spline volume function based on taper points...
#
#  taper = the taper data frame for the stem object
#  lenBot = length at the beginning of the section (default whole stem)
#  lenTop = length at the top of the section (default whole stem)
#  isLog = TRUE: a "downLog" object; FALSE: a "standingTree" object
#  units = units for the stem (doesn't really matter as csa conversion is the
#          same for both--see .StemEnv$baFactor above)
#
#  Please note: Always pass the entire log's taper data frame as the first argument,
#               even if you only want some intermediate bolt volume, as
#               the spline function is defined on the entirety of what is in
#               taper, and will not reflect the entire stem if only that portion
#               is passed that contains the bolt, for example.
#
#  Also: The reason we do not simply pass a "Stem" subclass object instead of
#        taper, units, etc. here is that we want this routine to possibly be of
#        use in the constructor of a "Stem" object. We could, of course, make
#        this a generic with "data.frame" and "Stem" signatures--see commented-out
#        section below...
#
#  Additionally, lenBot is the length to the bottom of the bolt to be integrate,
#  and lenTop is to the top of the bolt, both are relative to the butt of the
#  stem, which is always zero in a "Stem" subclass object. 
#
splineVolume = function(taper, lenBot, lenTop, isLog=TRUE, units='metric') {
    if(isLog)
      length = taper$length
    else
      length = taper$height
    if(lenBot > lenTop || lenBot < 0 || lenTop < 0 || lenTop > max(length))
      stop('Nonsensical lengths in splineVolume!')
    diameter = taper$diameter

#splineVolume = function(stemObject, lenBot, lenTop) {
  #  if(!is(stemObject, 'Stem'))
  #    stop('You must pass a \"Stem"\" subclass object to splineVolume!')
  #  taper = stemObject@taper
  #  if(is(stemObject, 'downLog'))
  #    length = taper$length
  #  else
  #    length = taper$height      
  #  diameter = taper$diameter
    
  #  if(lenBot > lenTop || lenBot < 0 || lenTop < 0 || lenTop > max(length))
  #    stop('Nonsensical lengths in splineVolume!')
  #  units = stemObject@units

    #cross-sectional area (ba) conversion factor and height to dbh...
    csaFactor =  ifelse(units==.StemEnv$msrUnits$English, .StemEnv$baFactor['English'],
                       .StemEnv$baFactor['metric'])
    
    csa.spline = splinefun(length, csaFactor*diameter^2)    #spline the cross-sectional area
    vol = integrate(csa.spline, lenBot, lenTop)$value       #integrates "hgt" from lenBot to lenTop
    return(vol)
} #splineVolume
assign('splineVolume', splineVolume, envir=.StemEnv)       #move to .StemEnv
environment(.StemEnv$splineVolume) = .StemEnv                   #assign its environment
rm(splineVolume)                                                #and remove from .GlobalEnv
  



#================================================================================
#
#  default surface area function based on default taper function...
#
#  botDiam = buttDiam for the log
#  topDiam = topTiam for the log
#  logLen = logLen for the log
#  solidType = solidType for the log
#  lenBot = length at the beginning of the section (default whole log)
#  lenTop = length at the top of the section (default whole log)
#
wbSurfaceArea = function(botDiam, topDiam, logLen, solidType, lenBot=0, lenTop=logLen) {
    if(lenBot > lenTop || lenBot < 0 || lenTop < 0)
      stop('Nonsensical lengths in wbSurfaceArea!')
    
    sa.taper = function(hgt, botDiam, topDiam, logLen, solidType) {
               diam = .StemEnv$wbTaper(botDiam, topDiam, logLen, 1, solidType, hgt)$diameter
               deriv = -2*(logLen-hgt)^(2/solidType-1) * (botDiam-topDiam)/(solidType*logLen^(2/solidType))
               tsa = pi*diam* sqrt(1 + deriv^2/4)              
               return(tsa)
    } #sa.taper
    sa = integrate(sa.taper, lenBot, lenTop, botDiam=botDiam, topDiam=topDiam, logLen=logLen,
                   solidType=solidType)$value
    return(sa)
} #wbSurfaceArea
assign('wbSurfaceArea', wbSurfaceArea, envir=.StemEnv)               #move to .StemEnv
environment(.StemEnv$wbSurfaceArea) = .StemEnv                       #assign its environment
rm(wbSurfaceArea)                                                    #and remove from .GlobalEnv
  



#================================================================================
#
#  spline surface area function based on taper points...
#
#  taper = the taper data frame for the log
#  lenBot = length at the beginning of the section (default whole log)
#  lenTop = length at the top of the section (default whole log)
#
#  Please note: Always pass the entire log's taper data frame as the first argument,
#               even if you only want some intermediate bolt surface area, as
#               the spline function is defined on the entirety of what is in
#               taper, and will not reflect the entire log if only that portion
#               is passed that contains the bolt, for example.
#
#  Additionally, lenBot is the length to the bottom of the bolt to be integrate,
#  and lenTop is to the top of the bolt, both are relative to the butt of the
#  log, which is always zero in a downLog object. 
#
splineSurfaceArea = function(taper, lenBot, lenTop, isLog=TRUE) {
    if(lenBot > lenTop || lenBot < 0 || lenTop < 0)
      stop('Nonsensical lengths in splineSurfaceArea!')
    if(isLog)
      length = taper$length
    else
      length = taper$height
    diameter = taper$diameter
    
    taper.spline = splinefun(length, diameter)
    sa.spline = function(hgt) {                           #"hgt" is length dummy argument
      diam = taper.spline(hgt)
      deriv = taper.spline(hgt, 1)
      ssa = pi*diam * sqrt(1 + deriv^2/4)
      return(ssa)
    } #sa.spline
    sa = integrate(sa.spline, lenBot, lenTop)$value       #integrates "hgt" from lenBot to lenTop
    return(sa)
} #splineSurfaceArea
assign('splineSurfaceArea', splineSurfaceArea, envir=.StemEnv)       #move to .StemEnv
environment(.StemEnv$splineSurfaceArea) = .StemEnv                   #assign its environment
rm(splineSurfaceArea)                                                #and remove from .GlobalEnv

  



#================================================================================
#
#  default coverage area function based on default taper function...
#
#  botDiam = buttDiam for the log
#  topDiam = topTiam for the log
#  logLen = logLen for the log
#  solidType = solidType for the log
#  lenBot = length at the beginning of the section (default whole log)
#  lenTop = length at the top of the section (default whole log)
#
wbCoverageArea = function(botDiam, topDiam, logLen, solidType, lenBot=0, lenTop=logLen) {
    if(lenBot > lenTop || lenBot < 0 || lenTop < 0)
      stop('Nonsensical lengths in wbCoverageArea!')
    
      if(identical(logLen, (lenBot - lenTop)))
        ca = (botDiam*solidType + 2*topDiam)*logLen/(solidType+2)
      else {
        r = solidType
        a = lenBot
        b = lenTop
        Du = topDiam
        Dl = botDiam
        L = logLen
        p1 = ((b - a)*Du*r + (2*b - 2*a)*Du)*L^(2/r)
        p2 = (L - b)^(2/r) * ((Du - Dl)*r*L + (b*Dl - b*Du)*r)
        p3 = (L - a)^(2/r) * ((Dl - Du)*r*L + (a*Du - a*Dl)*r)
        ca = (p1 + p2 + p3)/((r+2)*L^(2/r))
      }
    
    return(ca)
} #wbCoverageArea
assign('wbCoverageArea', wbCoverageArea, envir=.StemEnv)              #move to .StemEnv
environment(.StemEnv$wbCoverageArea) = .StemEnv                       #assign its environment
rm(wbCoverageArea)                                                    #and remove from .GlobalEnv



#================================================================================
#
#  spline coverage area function based on taper points...
#
#  taper = the taper data frame for the log
#  lenBot = length at the beginning of the section (default whole log)
#  lenTop = length at the top of the section (default whole log)
#
#  Please note: Always pass the entire log's taper data frame as the first argument,
#               even if you only want some intermediate bolt coverage area, as
#               the spline function is defined on the entirety of what is in
#               taper, and will not reflect the entire log if only that portion
#               is passed that contains the bolt, for example.
#
#  Additionally, lenBot is the length to the bottom of the bolt to be integrate,
#  and lenTop is to the top of the bolt, both are relative to the butt of the
#  log, which is always zero in a downLog object. 
#
splineCoverageArea = function(taper, lenBot, lenTop) {
    if(lenBot > lenTop || lenBot < 0 || lenTop < 0)
      stop('Nonsensical lengths in splineCoverageArea!')
    length = taper$length
    diameter = taper$diameter
    
    taper.spline = splinefun(length, diameter)
    ca = integrate(taper.spline, lenBot, lenTop)$value       #integrates "hgt" from lenBot to lenTop
    return(ca)
} #splineCoverageArea
assign('splineCoverageArea', splineCoverageArea, envir=.StemEnv)       #move to .StemEnv
environment(.StemEnv$splineCoverageArea) = .StemEnv                   #assign its environment
rm(splineCoverageArea)                                                #and remove from .GlobalEnv
 






#================================================================================
#
#  converts degrees to radians on [0,2pi]...
#
#  or use simply angle%%360 ??...
#
.deg2Rad = function(angle) {
    twoPi = 2.0*pi
    if(angle > 360) {
      fact = floor(angle/360)    
      angle = angle - fact*360
    }
    return(angle*twoPi/360)
} #.deg2Rad
assign('deg2Rad', .deg2Rad, envir=.StemEnv)               #move to .StemEnv
environment(.StemEnv$deg2Rad) = .StemEnv                  #assign its environment
rm(.deg2Rad)                                              #and remove from .GlobalEnv

#================================================================================
#
#  converts radians to degrees on [0,360]...
#

.rad2Deg = function(angle) {
    twoPi = 2.0*pi
    if(angle > twoPi) {
      fact = floor(angle/twoPi)    
      angle = angle - fact*twoPi
    }
    return(angle*360/twoPi)
} #.rad2Deg
assign('rad2Deg', .rad2Deg, envir=.StemEnv)               #move to .StemEnv
environment(.StemEnv$rad2Deg) = .StemEnv                  #assign its environment
rm(.rad2Deg)                                              #and remove from .GlobalEnv
  
  





#================================================================================
#
.underLine = function(lineLength = 0,   #length of line
                      lineChar = '-',   #character for line
                      prologue = '\n',  #some character (vector) to be output first
                      postfix = '\n'    #character (vector) to be output last
                     )
{
#------------------------------------------------------------------------------
#   just cat()'s an underline for nicer output...
#------------------------------------------------------------------------------
    if(lineLength>0 && lineLength<200)
      cat(prologue,rep(lineChar, lineLength),postfix,sep='')
    return(invisible())
}   #.underLine
assign('underLine', .underLine, envir=.StemEnv)           #move to .StemEnv
environment(.StemEnv$underLine) = .StemEnv                #assign its environment
rm(.underLine)                                            #and remove from .GlobalEnv






#================================================================================
#
#   just generates a random ID in character string that can be combined with
#   any prefix for spatial IDs in the package...
#
randomID = function(lenDigits=4,                        #number of digits to use
                    lenAlpha=4,                         #number of alpha characters
                    ...)
{
    if(lenDigits<1 || lenAlpha<1)
      stop('random IDs must have at least one digit and one letter!')
    dig = sample(0:9, lenDigits)
    alpha = sample(letters, lenAlpha)
    comb = c(dig, alpha)
    nn = length(comb)
    idx = sample(1:nn, nn)
    ranid = paste(comb[idx], collapse='')
    return(ranid)
}   #randomID
assign('randomID', randomID, envir=.StemEnv)              #move to .StemEnv
environment(.StemEnv$randomID) = .StemEnv                 #assign its environment
rm(randomID)                                              #and remove from .GlobalEnv


