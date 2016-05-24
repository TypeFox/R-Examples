chainsawSliver = function(downLog,
                 sect,
                 gLog,
                 nSegs = 25,             #bolt can be whole stem, so make enough
                 runQuiet = TRUE,
                 ...
                )
{
#---------------------------------------------------------------------------
#
#   This routine will estimate the volume of the intersection of the sample
#   plot with the log under protocol 1 for the chain saw method; i.e., where
#   the locus of all sample plot centers falls within the sausage inclusion
#   zone. The slivers will at least be tangent to the central longitudinal
#   axis/needle of the log at a minimum under this protocol.
#
#   However, the above is controlled by the calling routine, so actually this
#   will cut a slice out of the log for any intersection with the plot. Since
#   the intersection is passed as an argument, it is up to the calling routine
#   again to make sure that in fact, there is an intersection to begin with!
#
#   The main idea here is summarized in the following steps...
#     1. transform back so that the log and section are situated with the
#        buttD at x=0, facing east to make calculations simpler
#     2. determine the minimal bounding bolt on the log that encompases
#        the sliced section resulting from the intersection of the log
#        and the sample plot, as passed in gLog
#     3. determine the area of the section and of the bolt to get the proportion
#        of area the section is of the bolt
#     4. determine the bolt volume
#     5. assume that for the bolt section area and volume will share the same
#        proportion and apply the proportion from step 3 to the bolt volume
#        in step 4
#
#   The above is not going to be exact, there is some error associated with it
#   because the area is a projection and does not take the log shape into
#   consideration. But there seems to be no other way to get the sliver/section
#   volume. Update (17-June): It looks like this method works fine as it
#   is unbiased when all plots within the inclusion zone are considered. Indeed,
#   it is remarkably good.
#
#   24-Feb-2011: Added 'Length', 'surfaceArea', 'coverageArea' to the variables
#                that can be estimated with CS method. 
#
#   Arguments...
#     downLog = an object of class "downLog"
#     sect = a matrix of (x,y) points delineating the sliver/section that
#            the chain saw method P1 cuts out of the log in the form of
#            a closed polygon
#     gLog = a "gpc.polygon" object of the log itself
#     nSegs = number of segments for the taper curve (on one side) of the bounding
#             bolt segment (see above, step 2)
#     runQuiet = TRUE: no summary printin; FALSE: print some info
#     ... = passes/gobbles other args
#
#     Please note that if too few segments are specified in nSegs, the apportionment
#     will be off; for example, if the entire log is enclosed within the plot,
#     propArea could be greater than one. The more segments, the less the error.
#
#   Returns...
#     a list with...
#     --gBolt = a "gpc.polygon" object of the minimal bounding bolt
#     --rotBolt = a matrix representing the minimal bonding bolt
#     --boltVol = the bounding bolt volume
#     --sectVol = the sliver/section volume
#     --area = a vector of the proportional area, bolt ad section areas
#
#   This code was converted from the fixed-radius plot study to be used in
#   the sampSurf package, 30-Aug-2010, JHG.
#
#   14-April-2015: added explicit reference to rgeos in calling area.poly as
#   it was causing problems with example checks on CRAN for some reason.
#
#Author...									Date: 15-June-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#
#---------------------------------------------------------------------------
# 
    if(!(is(downLog,'downLog') || !validObject(downLog)) )
       stop('***>Invalid "downLog" object passed!')
    
#
#   tansformation matrix and its inverse...
#
    trMat = transfMatrix(downLog@logAngle, coordinates(downLog@location))
    trMatInv = solve(trMat)

    logLen = downLog@logLen
    buttDiam = downLog@buttDiam
    topDiam = downLog@topDiam
    solidType = downLog@solidType
    #if(is.null(solidType))
      #stop('***>chainSaw method requires use of the built-in taper function for log: solidType must be non-NULL!')

#    
#   transform log back so it is laying east with butt at x=0...
#
    halfLen = logLen/2
    log = as(gLog, 'matrix')
    log = cbind(log, rep(1, nrow(log)) )
    dn = dimnames(log)
    log = log %*% trMatInv
    dimnames(log) = dn
    log[,'x'] = log[,'x'] + halfLen    #shift so base is at x=0, still a closed polygon

#    
#   convert the slice-section the same way and make a gpc.poly from it...
#
    sect = cbind(sect, rep(1, nrow(sect)) )   #add hc's
    sect = sect %*% trMatInv
    dimnames(sect) = dn
    sect[,'x'] = sect[,'x'] + halfLen
    gSect = as(sect[,-3],'gpc.poly')
    sectArea = rgeos::area.poly(gSect)   #make this explicit on rgeos

#    
#   determine the minimal bounding bolt that will enclose the sliced section;
#   note that in very few cases (odd log angles) a slight rounding error can
#   be introduced by the transformations, so we need to catch that or it could
#   introduce NaNs...
#    
    sectLens = range(sect[,'x'])                        #min/max lengths for bounding bolt
    if(sectLens[1] < 0)                                 #catch little rounding error
      sectLens[1] = 0
    if(sectLens[2] > logLen)                            #again
      sectLens[2] = logLen
    boltLen = sectLens[2] - sectLens[1]                 #bounding bolt length
    length = seq(sectLens[1], sectLens[2], len=nSegs)   #lengths for the bolt section taper
    if(!is.null(solidType)) {
      taper = .StemEnv$wbTaper(buttDiam, topDiam, logLen, nSegs=nSegs, solidType, length)
      diameter = taper$diameter
    }
    else {
      diameter = taperInterpolate(downLog, 'diameter', length)
      taper = data.frame(diameter, length)
    }
   
#   from makeLog; get the entire bolt profile outline...    
    rad = diameter/2                                         #radius
    profile = data.frame( rad = c(-rad, rev(rad)) )      #center at diameter=x=0 
    profile$length = c(length, rev(length))              #duplicate lengths too
    profile = rbind(profile, profile[1,])                #close the polygon
    np = nrow(profile)
    
#    
#   lie the bolt on the ground, aligned with tip towards positive x-axis...
#
    rotBolt = profile
    rotBolt$hc = rep(1, np)
    rotBolt = as.matrix(rotBolt)
    trMat = transfMatrix(offset = c(0,-halfLen))   #to center at (0,0)
    rotBolt = rotBolt %*% trMat                    #translate
    trMat = transfMatrix(angle = -pi/2)            #lay the section down centered at (0,0)
    rotBolt = rotBolt %*% trMat
    
#
#   now rotate and translate back to the desired position of the original log...
#
    trMat = transfMatrix(downLog@logAngle, offset = coordinates(downLog@location)) 
    rotBolt = rotBolt %*% trMat
    dimnames(rotBolt) = list(NULL, c('x','y','hc'))
    gBolt = as(rotBolt[,-3], 'gpc.poly')

#
#   get the areas required for volume apportionment...
#
    boltArea = rgeos::area.poly(gBolt)   #make this explicit on rgeos
    propArea = sectArea/boltArea         #proportion of bolt area in the plot inter-section

#    
#   calculate bolt and sliver section volume, the latter assuming volume gets apportioned
#   the same as the sliver area w/r to the bolt...
#
    if(!is.null(solidType)) { 
      boltVol = .StemEnv$wbVolume(buttDiam, topDiam, logLen, solidType, sectLens[2]) -
                .StemEnv$wbVolume(buttDiam, topDiam, logLen, solidType, sectLens[1])
      boltSA = .StemEnv$wbSurfaceArea(buttDiam, topDiam, logLen, solidType,
                                      sectLens[1], sectLens[2])
      boltCA = .StemEnv$wbCoverageArea(buttDiam, topDiam, logLen, solidType, 
                                       sectLens[1], sectLens[2])
    }
    else {
      boltVol = .StemEnv$SmalianVolume(taper)$logVol           #note, logVol==whole bolt
      boltSA = .StemEnv$splineSurfaceArea(taper, sectLens[1], sectLens[2])
      boltCA = .StemEnv$splineCoverageArea(taper, sectLens[1], sectLens[2])
    }
    sectVol = propArea*boltVol
    sectSA = propArea*boltSA
    sectCA = propArea*boltCA
    boltBms = boltVol*downLog@conversions['volumeToWeight']
    sectBms = propArea*boltBms
    boltCarbon = boltBms*downLog@conversions['weightToCarbon']
    sectCarbon = propArea*boltCarbon
    

    if(!runQuiet) {
      cat('\nPercentage sliver is of bolt area =', propArea*100)
      cat('\nBolt volume (not expanded) =', boltVol)
      cat('\nSection/sliver volume (not expanded) =', sectVol)
      cat('\n')
    }  
    
    z = list(#gBolt = gBolt,                #"gpc.polygon" object of the minimal bounding bolt
             rotBolt = rotBolt,             #matrix outline of the minimal bounding bolt
             boltVol = boltVol,             #bounding bolt volume
             sectVol = sectVol,             #sliver section volume
             area = c(propArea=propArea,    #proportion sliver is to bolt w/r to area
                      boltArea=boltArea,    #bolt polygon area
                      sectArea=sectArea),   #sliver polygon area
             boltLen = boltLen,             #bolt length == section length
             boltSA = boltSA,               #bolt surface area
             sectSA = sectSA,               #sliver surface area
             boltCA = boltCA,               #bolt coverage area
             sectCA = sectCA,               #sliver coverage area
             boltBms = boltBms,             #bolt woody biomass  (can be NA)
             sectBms = sectBms,             #sliver woody biomass  (can be NA)
             boltCarbon = boltCarbon,       #bolt carbon mass  (can be NA)
             sectCarbon = sectCarbon        #sliver carbon  (can be NA)
            ) 
    return(z)
        
}   #chainsawSliver


    
