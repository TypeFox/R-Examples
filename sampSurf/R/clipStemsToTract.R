#---------------------------------------------------------------------------
#
#   Three sections are delineated in the code below...
#
#     A. Base routines--works for logs or trees...
#        1. StemContainer and SpatialPolygons -- does quick check "all-in" check
#        2. StemContainer and Tract -- main interface routine for logs & trees
#        3. StemContainer and bufferedTract -- main interface routine for logs & trees
#     B. Routine for "standingTrees", signatures...
#        1. standingTrees and SpatialPolygons
#     C. Routine for "downLogs", signatures...
#        1. downLogs and SpatialPolygons
#
#   The way the above work is: the user calls A.2 or A.3, these, in turn, call B.1 or C.1
#   for trees or logs, respectively; the first thing either of these routines does is
#   call A.1 for the overall conservative check on whether all stems are in the SP.
#
#   Arguments...
#     stems = a collection of stems, either "downLogs" or "standingTrees"
#     tract = a "Tract" or "bufferedTract" object
#     checkOnly = TRUE: just check the status of logs, don't delete or clip any;
#                 FALSE: remove and clip as needed
#     runQuiet = TRUE: no reporting; FALSE: report status
#     showPlot = (for logs only) TRUE: show a plot as it develops of clipping;
#                FALSE: no plot -- this is good for seeing what is going on, but
#                is of little use in the overall system
#   Returns...
#     A list with the following from the methods...
#       -- df = a data frame noting which stems are in, out or intersect (logs)
#               and where they intersect (N,S,E,W border); note especially that
#               this df always has as many rows as the ORIGINAL input stem
#               container has stems (df=NA from base A.1 routine)
#       -- status = a vector with total number of stems, number inside, outside
#                   and intersecting (logs)
#       -- stems = checkOnly: the original collection; !checkOnly: The new collection
#                   with any stems outside removed; and for logs, any that intersected
#                  are now new log objects clipped to the boundary
#                  (stems=NA from base A.1 routine)
#
#   Note: if logs are clipped to the boundary, the taper data frame was
#         clipped at the spot where the log needle intersects the boundary and
#         no new data other than the intersection diameter,length is
#         added. This means that the resulting taper data can be only two
#         measurements if the log gets truncated to something short. The reason
#         for this is that the taper data could come from measurements or another
#         taper equation (solidType=NULL) and it gets involved to try to
#         construct interpolated data with splines that would add more points
#         to the log. Fortunately, if a clipped log results in only a few
#         taper points, it is not very long so it does not matter.
#
#         I made the above decision because I felt it was better to just preserve
#         the original taper data in all cases rather than try to do too much
#         here without the user knowing what is going on. One can always add more
#         points to clipped logs if desired as all of the routines are available
#         to do this in the package and then re-insert these back into the stem
#         collection.
#
#   Note: Once a new stem collection has been made that is different from the
#         original one passed as an argument, we must recreate the collection
#         using the appropriate constructor so that the overall bbox and statistics
#         can be correctly adjusted to the clipped population.
#
#Author...									Date: 6-Sept-2013
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   generic definition...
#
setGeneric('clipStemsToTract',  
           function(stems, tract, ...) standardGeneric('clipStemsToTract'),
             signature = c('stems', 'tract')
            )



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#    Section A.  "StemContainer" methods
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


       
#================================================================================
#  A. StemContainer,SpatialPolygons...
#
setMethod('clipStemsToTract',
          signature(stems = 'StemContainer', tract = 'SpatialPolygons'),
function(stems, tract, checkOnly=TRUE, runQuiet=TRUE, ... )
{
#------------------------------------------------------------------------------
#
#  This base method does some simple checking for the others; e.g., an
#  overall simple check for all stems inside the SP. Note that the stems bbox
#  contains the overall bbox for the log polygons and the tree diameter polygons
#  Thus, if this passes with all inside, we can be sure the needle and tree
#  center will also be inside. Likewise, this can fail, but a further test
#  can show that all the needles and centers do in fact lie within. So it's
#  just something to make a quick getaway if it works.
#
#  Please don't call this method directly unless you are adding to the system,
#  use the wrapper methods that call this correctly to avoid problems.
#
#------------------------------------------------------------------------------  
#
#   check for stem collection bbox lying all inside the tract boundary polygon;
#   if they all lie within here, then the needle or center will also lie within
#   and nothing else would be required...
#
#   note: we are checking any portion here w/t the bbox, but in the looking at each
#         stem we will check the needle (log) or centerpoint (tree) in the other
#         methods as required...
#
    bb.stem = bbox(stems)
    bb.tr = bbox(tract)
    
    if(extends(class(stems), 'downLogs'))      #are they logs?
      nStems = length(stems@logs)
    else                                       #must be standingTrees
      nStems = length(stems@trees)

    #note that this will check x-with-x and y-with-y correctly...
    #check is.na(inside) [failed, not all inside] on return...
    if(all(bb.stem[,'min']>=bb.tr[,'min']) && all(bb.stem[,'max']<=bb.tr[,'max'])) 
      status = c(total=nStems, inside=nStems, outside=0, intersect=NA) #good, done!
    else
      status = c(total=nStems, inside=NA, outside=NA, intersect=NA)    #more checking needed

    return(list(df = NA, status = status, stems=NA))
  
}   #signature: StemContainer,SpatialPolygons method
)   #setMethod



       
#================================================================================
#  A.2. StemContainer,Tract...
#
#       Interface for "downLogs" or "standingTrees" plus "Tract"
#
setMethod('clipStemsToTract',
          signature(stems = 'StemContainer', tract = 'Tract'),
function(stems, tract, checkOnly=FALSE, runQuiet=TRUE, ... )
{
#------------------------------------------------------------------------------
#
#   first make sure the units are compatible...
#
    if(stems@units != tract@units)
      stop('Units for Stem collection and Tract are not compatible!')

    resList = clipStemsToTract(stems, perimeter(tract), checkOnly=checkOnly,
                               runQuiet=runQuiet, ...)

    return(resList)
  
}   #signature: StemContainer,Tract method
)   #setMethod



       
#================================================================================
#  A.3. StemContainer, bufferedTract...
#
#       Interface for "downLogs" or "standingTrees" plus "bufferedTract"
#
setMethod('clipStemsToTract',
          signature(stems = 'StemContainer', tract = 'bufferedTract'),
function(stems, tract, checkOnly=FALSE, runQuiet=TRUE, clipToBuffer=TRUE, ... )
{
#------------------------------------------------------------------------------
#
#   first make sure the units are compatible...
#
    if(stems@units != tract@units)
      stop('Units for Stem collection and Tract are not compatible!')

#
#   This will clip trees to the internal buffer by default, or to the
#   tract boundary if desired...
#
    if(clipToBuffer)
      resList = clipStemsToTract(stems, tract@spBuffer, checkOnly=checkOnly,
                                 runQuiet=runQuiet, ...)
    else
      resList = callNextMethod(stems, tract, checkOnly, runQuiet, ...)

    return(resList)
  
}   #signature: StemContainer,bufferedTract method
)   #setMethod








#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#    Section B.  "standingTrees" methods
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


       
#================================================================================
#  B.1. standingTrees,SpatialPolygons...
#
setMethod('clipStemsToTract',
          signature(stems = 'standingTrees', tract = 'SpatialPolygons'),
function(stems, tract, checkOnly=FALSE, runQuiet=TRUE, ... )
{
#------------------------------------------------------------------------------
#
#   For trees, just check that the centerpoint is in the tract...
#
#------------------------------------------------------------------------------
#
#   general check first (A.1 StemContainer,SpatialPolygons) to see if there's
#   anything to do...
#
    resList = callNextMethod(stems, tract, checkOnly, runQuiet, ...)
    if(!is.na(resList$status['inside'])) #done
      return(resList)

#
#   okay, get to work now...
#
    nStems = length(stems@trees)
    bb.tr = bbox(tract)

#
#   simply check to see if centerpoint location of each stem is inside
#   the bounding box for the SpatialPolygon object passed; if a tree is on
#   the polygon, it is counted as being inside...
#
    colNames = c('id', 'inside', 'outside')                             #intersect=NA below
    df = as.data.frame(matrix(FALSE, nrow=nStems, ncol=length(colNames)))
    colnames(df) = colNames
    rownames(df) = names(stems@trees)
    for(i in seq_len(nStems)) {
      stem = stems@trees[[i]]
      df[i, 'id'] = getID(stem)
      cp = coordinates(stem@location)
      if(cp[1,'x'] <= bb.tr['x','max'] &&   #east
         cp[1,'x'] >= bb.tr['x','min'] &&   #west
         cp[1,'y'] <= bb.tr['y','max'] &&   #north
         cp[1,'y'] >= bb.tr['y','min'])     #south
        df[i, 'inside'] = TRUE
      else
        df[i, 'outside'] = TRUE
    }
 
    nOutside = sum(df[, 'outside'])
    nInside = sum(df[, 'inside'])
    status = c(total=nStems, inside=nInside, outside=nOutside, intersect=NA)
    if(!runQuiet) {
      cat('\nTotal stems =', nStems)
      cat('\nStems inside =', nInside)
      cat('\nStems outside =', nOutside)
      cat('\n')
    }

#    
#   if desired, see if we need to delete/remove any "outside" trees...
#
    if(!checkOnly && nOutside > 0) {    #keep df as is, don't delete for the record
      ddx = which(df$outside)           #indexes of those trees outside
      stems@trees[ddx] = NULL           #remove these trees
    }

#
#   we must re-create the collection to get the statistics, bbox, etc. correct...
#
    stems = standingTrees(stems@trees)
    
    return(invisible(list(df = df, status = status, stems = stems)))      
  
}   #signature: standingTrees,SpatialPolygons method
)   #setMethod











#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#    Section C.  "downLogs" methods
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


 
#================================================================================
#  C.1. downLogs, SpatialPolygons...
#
setMethod('clipStemsToTract',
          signature(stems = 'downLogs', tract = 'SpatialPolygons'),
function(stems, tract, checkOnly=FALSE, runQuiet=TRUE, showPlot=FALSE, ... )
{
#------------------------------------------------------------------------------
#
#   For logs, main check for the "needle" all within the tract...
#
#------------------------------------------------------------------------------
#
#   general check first (A.1 StemContainer,SpatialPolygons) to see if there's
#   anything to do...
#
    resList = callNextMethod(stems, tract, checkOnly, runQuiet, ...)
    if(!is.na(resList$status['inside'])) #done
      return(resList)
    
#
#   use rgeos for the log intersection tests in what follows, might as well load it now...
#
    require(rgeos)

#
#   okay, get to work now...
#
    nStems = length(stems@logs)
    bb.tr = bbox(tract)

#
#   we want to know which side of the tract the log intersects...
#
    tract.sl = as(tract, "SpatialLines")    #cast from SP to SL
    xy = coordinates(tract.sl)[[1]][[1]]    #now the coordinates for each side are rows in a matrix
    west = xy[1:2,]
    West = SpatialLines(list(Lines(Line(west),'west')))
    north = xy[2:3,]
    North = SpatialLines(list(Lines(Line(north),'north')))
    east = xy[3:4,]
    East =  SpatialLines(list(Lines(Line(east),'east')))
    south = xy[4:5,]                        #tract is a closed polygon, so 5==1
    South = SpatialLines(list(Lines(Line(south),'south')))


#
#   check loop for logs--no intersection computations yet...
#
    cardinal = c('north','south','east','west')
    colNames = c('id', 'inside', 'outside', 'intersect', cardinal)
    df = as.data.frame(matrix(FALSE, nrow=nStems, ncol=length(colNames)))
    colnames(df) = colNames
    rownames(df) = names(stems@logs)
    for(i in seq_len(nStems)) {
      stem = stems@logs[[i]]
      df[i, 'id'] = getID(stem)
      df[i, 'intersect'] = FALSE
      
      logNeedle = bbox(stem@slNeedleAxis)           #needle coords bbox for convenience 
      #see if it is completely outside the tract (another check required after intersection below)...
      if(logNeedle['x','min'] > bb.tr['x','max'] || #east
         logNeedle['x','max'] < bb.tr['x','min'] || #west
         logNeedle['y','min'] > bb.tr['y','max'] || #north
         logNeedle['y','max'] < bb.tr['y','min']) { #south
        df[i, 'outside'] = TRUE
        next
      }

      #see if it intersects a boundary...
      df[i, 'east'] = gIntersects(stem@slNeedleAxis, East)
      df[i, 'west'] = gIntersects(stem@slNeedleAxis, West)
      df[i, 'north'] = gIntersects(stem@slNeedleAxis, North)
      df[i, 'south'] = gIntersects(stem@slNeedleAxis, South)
      if(any(unlist(df[i, cardinal])))
        df[i, 'intersect'] = TRUE
      #required extra check to catch logs straddling (diagonal) just outside corners
      #but not intersecting, these will fail in the above "out" check...
      else if(logNeedle['x','max'] > bb.tr['x','max'] || #east
              logNeedle['x','min'] < bb.tr['x','min'] || #west
              logNeedle['y','max'] > bb.tr['y','max'] || #north
              logNeedle['y','min'] < bb.tr['y','min'])   #south
        df[i, 'outside'] = TRUE
      else
        df[i, 'inside'] = TRUE        
    } #check logs loop

    nOutside = sum(df[, 'outside'])
    nInside = sum(df[, 'inside'])
    nIntersect = sum(df[, 'intersect'])
    status = c(total=nStems, inside=nInside, outside=nOutside, intersect=nIntersect)
    if(!runQuiet) {
      cat('\nTotal stems =', nStems)
      cat('\nStems inside =', nInside)
      cat('\nStems outside =', nOutside)
      cat('\nStems intersecting =', nIntersect)
      cat('\n')
    }

    if(checkOnly)    
      return(invisible(list(df = df, status = status, stems = stems)))

    if(showPlot) {
      plot(tract, axes=TRUE)
      plot(stems, add=TRUE)
    }

#
#   logs are much more complicated, we can first remove those that lie completely outside
#   the tract, but then we must clip those that intersect the tract boundaries; and
#   a given log can intersect more than one boundary, so that we may need to clip more
#   than once; thus, most of the code in the for loop must be within the directional sections...
#
    
    #first, see if we need to delete any logs...
    if(nOutside > 0) {    #keep df as is, don't delete for the record
      ddx = which(df$outside)           #indexes of those logs outside
      stems@logs[ddx] = NULL            #remove these logs
    }

    #now, see if we need to clip any logs--if not, just return...
    if(nIntersect == 0) {
      stems = downLogs(stems@logs)           #recreate the collection
      return(invisible(list(df = df, status = status, stems = stems)))
    }

    #main loop through all logs that intersect any boundary...
    logNames = rownames(df)                  #names for remaining logs
    idx = which(df$intersect)                #indexes for stems that intersect from full df
    for(i in idx) {                          #clip each log that intersects            
      stem = stems@logs[[ logNames[i] ]]
      logAngle = stem@logAngle               #stays constant, regardless of clipping, unlike length, etc.
      direction = names(which(unlist(df[i, cardinal]))) #vector: can intersect more than once
      
      #loop through NSEW, and operate on those that match direction intersections...
      for(nsew in cardinal) {
        if(!any(direction == nsew) || is.null(stem))  #any intersection in this direction?
          next       
                
        needle = stem@slNeedleAxis
        xy.n = coordinates(needle)[[1]][[1]]                  #needle coordinates==endpoints
        if(nsew == 'north') {
          xy = coordinates(gIntersection(needle, North))[1,]  #coordinates returns a matrix
          xy.in = xy.n[which.min(xy.n[,'y']),]                #inside needle coords == min 'y' -- drop to vector
          testCondition = logAngle > pi                       #tests for butt outside the tract
        }
        else if(nsew == 'south') {
          xy = coordinates(gIntersection(needle, South))[1,]
          xy.in = xy.n[which.max(xy.n[,'y']),]                #inside needle coords == max 'y' -- drop to vector
          testCondition = logAngle < pi
        }
        else if(nsew == 'east') {
          xy = coordinates(gIntersection(needle, East))[1,]
          xy.in = xy.n[which.min(xy.n[,'x']),]                #inside needle coords == min 'x' -- drop to vector
          testCondition = cos(logAngle) < 0
        } 
        else if(nsew == 'west') {
          xy = coordinates(gIntersection(needle, West))[1,]
          xy.in = xy.n[which.max(xy.n[,'x']),]                #inside needle coords == max 'x' -- drop to vector
          testCondition = cos(logAngle) > 0
        }
        #get the length of the log for the portion inside the tract...
        len = sqrt( (xy.in['x']-xy['x'])^2 + (xy.in['y']-xy['y'])^2 ) 
        
        if(testCondition) {                            #butt section is outside
          buttIn = FALSE
          len = stem@logLen - len                      #get outside length--we will clip the butt off
        }
        else                                           #butt section is inside
          buttIn = TRUE
        diam = taperInterpolate(stem, 'diameter', len) #get the diameter for the length
        tdx = findInterval(len, stem@taper$length)
        nTaper = nrow(stem@taper)
        if(buttIn) 
          taper = rbind(stem@taper[1:tdx,], c(diam,len))
        else {
          taper = rbind(c(diam,len), stem@taper[(tdx+1):nTaper,])
          taper$length = taper$length - len  #make length at butt==0
        }
        #first, we need a new center point for the new length, then reconstruct the log...
        dx = (xy['x'] - xy.in['x'])          #x-offset 
        dy = (xy['y'] - xy.in['y'])          #y-offset
        x = xy['x'] - dx/2               
        y = xy['y'] - dy/2
        if(nrow(taper) >=2)   #this should not fail...
          stem = downLog(taper, stem@solidType, logAngle, centerOffset = c(x, y),
                         species = stem@species,
                         logID = getID(stem),
                         vol2wgt = stem@conversions['volumeToWeight'],
                         wgt2carbon = stem@conversions['weightToCarbon'],
                         description = stem@description,
                         units = stem@units,
                         spUnits = stem@spUnits,
                         userExtra = stem@userExtra
                        )
        else
          stem = NULL
        
        if(showPlot) {
          if(nsew %in% cardinal[1:2])
            plot(stem, logCol=transparentColorBase('red',0.5),add=TRUE)
          else
            plot(stem, logCol=transparentColorBase('green',0.5),add=TRUE)
        }
        
      } #nsew loop

      stems@logs[[ logNames[i] ]] = stem
      
    } #clip intersect [i] loop

#
#   we must re-create the collection to get the statistics, bbox, etc. correct...
#
    stems = downLogs(stems@logs)

    return(invisible(list(df = df, status = status, stems = stems)))
  
}   #signature: downLogs,SpatialPolygons method
)   #setMethod

