sampleLogs = function(nLogs = 2,
                      buttDiams = c(8, 40),   #cm or inches in default
                      topDiams = c(0, 0.9),   #proportion multiplier
                      logLens = c(1,10),
                      logAngles = c(0, 2*pi),
                      solidTypes = c(1,10),
                      species = .StemEnv$species,
                      sampleRect = NULL,
                      startSeed = NA,
                      runQuiet = FALSE,
                      ...
                     )
{
#---------------------------------------------------------------------------
#
#   This routine will draw random logs from uniform distributions with
#   the specs given by the arguments. This routine should work equally
#   for English and metric units, just specify the ranges desired that make
#   sense in either system. In addition, you can assume any units for the
#   dimensions based on what you specify for logical ranges.
#
#   There is little error checking done here, please follow the guidelines
#   in the arguments or things could be a bit unexpected.
#
#   Note that log volumes will be calculated in downLog() or similar
#   constructor.
#
#   Arguments...
#     nLogs = the number of sample logs to draw in the population
#     buttDiams = the large-end diameter range as a vector of (min,max)
#     topDiams = a range of proportional multipliers to the butt diameters
#                which will have been previously drawn: (min,max) in [0,1]
#     logLens = the range in log lengths as (min,max)
#     logAngles = range of log angles usually (0, 2*pi), but for methods
#                 like sausage where the orientation only matters to 180
#                 it could be (0, pi)
#     solidTypes = a range of (min,max) solid types for the taper equation
#                  see .StemEnv$wbTaper() for details
#     sampleRect = a matrix with row names c('x','y') and column names
#                  c('min', 'max') giving the bounding rectangle from which
#                  the log centers will be drawn for possible placement
#                  into a larger grid; if NULL or NA, the log centers in
#                  x,y will be dimensionless from (0,1)
#     startSeed = see initializeRandSeed for details; if NA, use the
#                 current seed; if a number, then this will be used and
#                 is the way to make duplicate log populations from the
#                 same rn stream for different simulation runs when
#                 comparing, e.g., sausage and standUp methods
#     ... = ignored arguments right now--just gobbled into the gloom
#
#   Returns...
#     A data frame with the following columns invisibly...
#       --species
#       --logLen = log lengths
#       --buttDiam = large-end/butt diameters
#       --topDiam = small-end/top diameters
#       --solidType = solid of revolution type
#       --cx = x center point
#       --cy = y center point
#       --logAngle = log angles in radians
#       --logAngle.D = log angles in degrees (for dummies like me)
#
#
#Author...									Date: 1-June-2010
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
#   some checks...
#
    if(!is.numeric(nLogs) || length(nLogs)>1 || nLogs<=0 || nLogs>10000 )
      stop('Number of logs is specified incorrectly: ',nLogs)
    nLogs = trunc(nLogs)
    if(length(buttDiams)!=2 || diff(buttDiams)<0 || any(buttDiams<0) )
      stop('Butt diameters specified incorrectly: ', buttDiams[1],',',buttDiams[2])
    if(length(topDiams)!=2 || diff(topDiams)<0 || any(topDiams<0) || any(topDiams>1) )
      stop('Top diameter proportions specified incorrectly: ', topDiams[1],',',topDiams[2])
    if(length(logAngles)!=2 || diff(logAngles)<0 || any(logAngles<min(.StemEnv$logAngles)) ||
       any(logAngles>max(.StemEnv$logAngles)) )
      stop('Log angles specified incorrectly: ', logAngles,'must be within:',.StemEnv$logAngles)
    if(length(solidTypes)!=2 || diff(solidTypes)<0 || any(solidTypes < 1) || any(solidTypes > 10) )
      stop('Solid types specified incorrectly: ', solidTypes[1],',',solidTypes[2])
              
  
#
#   initialize the random seed...
#
    ranSeed = initRandomSeed(startSeed)
  
#
#   generate the logs...
#
    logLen = round( runif(nLogs, logLens[1], logLens[2]), 2 )
    logAngle = runif(nLogs, logAngles[1], logAngles[2])
    solidType = round( runif(nLogs, solidTypes[1], solidTypes[2]), 1 )
    buttDiam = round( runif(nLogs, buttDiams[1], buttDiams[2]), 2)
    topDiam = round( runif(nLogs, topDiams[1], topDiams[2])*buttDiam, 2)

#
#   species codes, names, group, etc.,
#
    species = sample(species, nLogs, replace = TRUE)

#
#   generate log centers within the sample rectangle...
#
    if(is.null(sampleRect) || is.na(sampleRect)) {
      sampleRect = matrix(rep(0:1, 2), nrow=2,                       #dimensionless center points
                          byrow = TRUE,
                          dimnames = list(c('x','y'),c('min','max'))
        
                         )
      if(!runQuiet)
        cat('\nNote: logs generated within [0,1] bbox!\n')
    }
    stopifnot(bboxCheck(sampleRect))     #make sure sampleRect is a bbox
    
    cx = runif(nLogs, sampleRect['x','min'], sampleRect['x','max'])
    cy = runif(nLogs, sampleRect['y','min'], sampleRect['y','max'])
    
      
#
#   wrap it up and send it home...
#
    logs = data.frame(species = species,
                      logLen = logLen,
                      buttDiam = buttDiam,
                      topDiam = topDiam,
                      solidType = solidType,
                      x = cx,
                      y = cy,
                      logAngle = logAngle,
                      logAngle.D = logAngle*180/pi,   #log angle in degrees
                      stringsAsFactors = FALSE
                     )
    return(logs)
}   #sampleLogs

   
