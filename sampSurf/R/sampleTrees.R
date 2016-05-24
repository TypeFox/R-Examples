sampleTrees = function(nTrees = 2,
                       dbhs = c(8, 40),   #cm or inches in default
                       topDiams = c(0.4, 0.9),   #proportion multiplier
                       heights = c(5, 15),
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
#   This routine will draw random trees from uniform distributions with
#   the specs given by the arguments. This routine should work equally
#   for English and metric units, just specify the ranges desired that make
#   sense in either system. In addition, you can assume any units for the
#   dimensions based on what you specify for logical ranges.
#
#   There is little error checking done here, please follow the guidelines
#   in the arguments or things could be a bit unexpected.
#
#   Note that tree volumes will be calculated in standingTree() or similar
#   constructor.
#
#   Arguments...
#     nTrees = the number of sample trees to draw in the population
#     dbhs = the breast height diameter range as a vector of (min,max)
#     topDiams = a range of proportional multipliers to the breast hgt diameters
#                which will have been previously drawn: (min,max) in [0,1]
#     heights = the range in tree heights as (min,max)
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
#       --height = tree heights
#       --dbh = breast height diameters
#       --topDiam = top diameters
#       --solidType = solid of revolution type
#       --x = x center point
#       --y = y center point
#
#
#Author...									Date: 25-Oct-2011
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
    if(!is.numeric(nTrees) || length(nTrees)>1 || nTrees<=0 || nTrees>10000 )
      stop('Number of logs is specified incorrectly: ',nTrees)
    nTrees = trunc(nTrees)
    if(length(dbhs)!=2 || diff(dbhs)<0 || any(dbhs<0) )
      stop('Breast height diameters specified incorrectly: ', dbhs[1],',',dbhs[2])
    if(length(topDiams)!=2 || diff(topDiams)<0 || any(topDiams<0) || any(topDiams>1) )
      stop('Top diameter proportions specified incorrectly: ', topDiams[1],',',topDiams[2])
    if(length(solidTypes)!=2 || diff(solidTypes)<0 || any(solidTypes < 1) || any(solidTypes > 10) )
      stop('Solid types specified incorrectly: ', solidTypes[1],',',solidTypes[2])
              
  
#
#   initialize the random seed...
#
    ranSeed = initRandomSeed(startSeed)
  
#
#   generate the logs...
#
    height = round( runif(nTrees, heights[1], heights[2]), 2 )
    solidType = round( runif(nTrees, solidTypes[1], solidTypes[2]), 1 )
    dbh = round( runif(nTrees, dbhs[1], dbhs[2]), 2)
    topDiam = round( runif(nTrees, topDiams[1], topDiams[2])*dbh, 2)

#
#   species codes, names, group, etc.,
#
    species = sample(species, nTrees, replace = TRUE)

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
    
    cx = runif(nTrees, sampleRect['x','min'], sampleRect['x','max'])
    cy = runif(nTrees, sampleRect['y','min'], sampleRect['y','max'])
    
      
#
#   wrap it up and send it home...
#
    trees = data.frame(species = species,
                       height = height,
                       dbh = dbh,
                       topDiam = topDiam,
                       solidType = solidType,
                       x = cx,
                       y = cy,
                       stringsAsFactors = FALSE
                      )
    return(trees)
}   #sampleTrees

   
