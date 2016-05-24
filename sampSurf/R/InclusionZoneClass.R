#---------------------------------------------------------------------------
#
#   This file holds the S4 class definitions for the Inclusion Zone method
#   related classes. Inclusion zones are properties of the particular
#   'Stem' subclass, and 'ArealSampling' subclass.
#
#   The following section holds all of the definitions for down logs, those
#   for standing trees are after these in the file.
#
#   Classes...
#     1. InclusionZone: virtual for all 'ArealSampling' methods and Stem
#                       classes
#     2. downLogIZ: virtual for all 'ArealSampling' methods and 'downLog'
#                   subclasses
#     3. standUpIZ: the standup method of sampling down cwd
#     4. chainSawIZ; the chain saw method for sampling down cwd
#     5. sausageIZ: the sausage method of sampling down cwd
#     6. fullChainSawIZ: the FCS has sausage IZ
#     7. pointRelascopeIZ: the point relascope method for sampling down cwd
#     8. perpendicularDistanceIZ: normal PDS
#     9. omnibusPDSIZ: for omnibus estimation under normal PDS
#    10. distanceLimitedIZ: for canonical distance limited sampling (dls)
#    11. distanceLimitedMCIZ: DLMC
#    12. distanceLimitedPDSIZ: for canonical DLPDS
#    13. omnibusDLPDSIZ: combines omnibus PDS component with DLMC component
#    14. hybridDLPDSIZ: DLMC & canonical PDS
#
#   Search for "standingTreeIZ" to find the definitions for standing trees at
#   the end of this file.
#
#Author...									Date: 20-Aug-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#




#=================================================================================================
#
#  1. define the virtual InclusionZone class...
#
setClass('InclusionZone',
         
#
#  slots for the class and its subclasses...
#
    representation(description = 'character',
                   units = 'character',              #English or metric units
                   bbox = 'matrix',                  #overall bounding box
                   spUnits = 'CRS',                  #spatial units
                   puaBlowup = 'numeric',            #per unit area expansion factor
                   puaEstimates = 'list',            #per unit area estimates
                   userExtra = 'ANY'                 #anything else the user wants to include--no checks
                  ),
    prototype = list(description = 'Areal sampling inclusion zone',  #some defaults for validity checking
                     units = .StemEnv$msrUnits$metric,
                     bbox = matrix(rep(0,4), nrow=2, dimnames=list(c('x','y'), c('min','max'))),
                     spUnits = CRS(projargs=as.character(NA)),
                     puaBlowup = 1,                  #does no harm
                     puaEstimates = list(),          #empty list
                     userExtra = NULL
                    ),
    contains = 'VIRTUAL',
    validity = function(object) {
                 if(!(object@units %in% .StemEnv$msrUnits))
                   return('units of measure must be "English" or "metric"')

                 if(length(object@puaBlowup) > 1)
                   return('puaBlowup must be a scalar!')
                 if(!is.na(object@puaBlowup) && object@puaBlowup < 0)  #it can be NA for DLPDS total
                   return('puaBlowup factor must be non-negative!')

                 if(length(object@puaEstimates > 0)) {
                   estNames = match(names(object@puaEstimates), unlist(.StemEnv$puaEstimates))
                   if(any(is.na(estNames)))
                     return('invalid slots in puaEstimates list object!')
                 }

                 #essentially the same checks as in bboxCheck()...
                 if(!class(object@bbox) == 'matrix')
                   return('bbox slot must be a 2x2 matrix')
                 bboxNames = match(rownames(object@bbox), c('x','y'))
                 if(any(is.na(bboxNames)))
                   return('slot bbox rownames must be "x", "y"!')
                 bboxNames = match(colnames(object@bbox), c('min','max'))
                 if(any(is.na(bboxNames)))
                   return('slot bbox colnames must be "min", "max"!')
                 if(any( apply(object@bbox,1,function(x) if(x['min'] >= x['max']) TRUE else FALSE) ))
                   return('in slot bbox, "min" must be less than "max" for x and y!')

                 if(!is.na(object@spUnits@projargs) && object@spUnits@projargs == '+proj=longlat')
                   return(paste('spUnits must be commensurate with units,',
                                'please convert to non-geographic coordinate system!')
                         )
                   
                 return(TRUE)
               } #validity check
) #class InclusionZone 
                     
                   



#=================================================================================================
#
#  2. the down log inclusion zone class is just a direct descendant of 'InclusionZone'...
#
#
setClass('downLogIZ',
    representation(downLog = 'downLog'              #downLog object
                  ),
    #prototype = list(#downLog = new('downLog'),         #some defaults for validity checking
    #                 units = .StemEnv$msrUnits$metric
    #                ),
    contains = c('InclusionZone', 'VIRTUAL'),             #a subclass of the virtual 'InclusionZone' class
    validity = function(object) {

                 if(!validObject(object@downLog))
                   return('object in downLog slot is an invalid object!')

                 if(!identical(object@downLog@units, object@units))
                   return('measurement units do not match between downLog and IZ objects!')

                 if(!identical(object@downLog@spUnits, object@spUnits))
                   return('Spatial units do not match between downLog and IZ objects!')
                 
                 return(TRUE)
               } #validity check
) #class downLogIZ 


 





#=================================================================================================
#
#  3. the standup class is just a direct descendant of 'downLogIZ'...
#
#
setClass('standUpIZ',
    representation(circularPlot = 'circularPlot'              #circularPlot object
                  ),
    #prototype = list(#downLog = new('downLog'),         #some defaults for validity checking
    #                ),
    contains = 'downLogIZ',                             #a subclass of the virtual 'downLogIZ' class
    validity = function(object) {

                 if(!validObject(object@circularPlot))
                   return('object in circularPlot slot is an invalid object!')

                 if(!identical(object@downLog@units, object@circularPlot@units))
                   return('measurement units do not match between downLog and circularPlot objects!')

                 if(!identical(object@downLog@spUnits, object@circularPlot@spUnits))
                   return('Spatial units do not match between downLog and circularPlot objects!')
                 
                 return(TRUE)
               } #validity check
) #class standUpIZ 




#=================================================================================================
#
#  4. the 'chainSawIZ' class is just a direct descendant of 'downLogIZ'...
#
#
setClass('chainSawIZ',
    representation(circularPlot = 'circularPlot',              #circularPlot object
                   sliver = 'SpatialPolygons',                 #the intersection sliver w/in the log
                   bolt = 'list'                               #the bolt/sliver information
                  ),
    #prototype = list(
    #                ),
    contains = 'downLogIZ',                        #a subclass of virtual 'downLogIZ' class
    validity = function(object) {

                 if(!validObject(object@circularPlot))
                   return('object in circularPlot slot is an invalid object!')

                 if(!identical(object@downLog@units, object@circularPlot@units))
                   return('measurement units do not match between downLog and circularPlot objects!')

                 if(!identical(object@downLog@spUnits, object@circularPlot@spUnits))
                   return('Spatial units do not match between downLog and circularPlot objects!')
                 
                 return(TRUE)
               } #validity check
) #class chainSawIZ 


 





#=================================================================================================
#
#  5. the sausage class is a direct descendant of 'downLogIZ'...
#
#
setClass('sausageIZ',
    representation(sausage = 'matrix',              #holds the sausage incl zone in matrix form
                   radius = 'numeric',              #plot radius for sausage
                   area = 'numeric',                #exact area of the inclusion zone
                   perimeter = 'SpatialPolygons',   #sausage perimeter in 'SpatialPolygons' form
                   pgSausageArea = 'numeric'        #polygon sausage area approximation
                  ),
    #prototype = list(#downLog = new('downLog'),         #some defaults for validity checking
    #                ),
    contains = 'downLogIZ',                                  #a subclass of the 'downLogIZ' class
    validity = function(object) {
                 if(object@radius <= 0 || is.na(object@radius))
                   return('plot radius must be positive non-missing!')

                 if(object@area < 0 || is.na(object@area))
                   return('object has negative or missing inlusion zone area!')

                 if(object@pgSausageArea < 0 || is.na(object@pgSausageArea))
                   return('object has negative or missing polygon area!')
                 
                 return(TRUE)
               } #validity check
) #class sausageIZ 





#=================================================================================================
#
#  6. the fullChainSaw class is a direct descendant of 'sausageIZ'...
#
#     This class was formalized with associated changes in the constructors, izGrid constructors,
#     and sampSurf 30-Sept-2013. This was done to get rid of the special code (since it just uses
#     a sausage IZ) and make the package totally OOP in this regard.
#
setClass('fullChainSawIZ',
    contains = 'sausageIZ',                                  #a subclass of the 'sausageIZ' class
) #class fullChainSawIZ 








#=================================================================================================
#
#  7. the pointRelascope class is a direct descendant of 'downLogIZ'...
#
#     I'll call the "mastercard" double/dual circle (blob) just "dual"...
#
#     below: hc = homogeneous coordinates
#
#     added: 13-Jan-2011
#
setClass('pointRelascopeIZ',
    representation(prs = 'pointRelascope',          #point relascope sampling object
                   dualCircle = 'matrix',           #holds the blob perimeter in matrix form w/ hc
                   radius = 'numeric',              #radius for each of the half blobs
                   area = 'numeric',                #exact area of the inclusion zone
                   perimeter = 'SpatialPolygons',   #blob perimeter in 'SpatialPolygons' form
                   pgDualArea = 'numeric',          #polygon blob area approximation
                   dualCenters = 'matrix'           #centers of each dual circle, w/o hc
                  ),
    contains = 'downLogIZ',                         #a subclass of the 'downLogIZ' class
    validity = function(object) {
                 if(object@radius <= 0 || is.na(object@radius))
                   return('plot radius must be positive non-missing!')

                 if(object@area < 0 || is.na(object@area))
                   return('object has negative or missing inclusion zone area!')

                 if(object@pgDualArea < 0 || is.na(object@pgDualArea))
                   return('object has negative or missing polygon area!')

                 if(!identical(object@prs@units, object@units))
                   return('object units do not match the pointRelascope object units')
                 
                 return(TRUE)
               } #validity check
) #class pointRelascopeIZ 









#=================================================================================================
#
#  8. the perpendicularDistanceIZ class is a direct descendant of 'downLogIZ'...
#
#
#     below: hc = homogeneous coordinates
#
#     added: 26-Jan-2011
#
setClass('perpendicularDistanceIZ',
    representation(pds = 'perpendicularDistance',   #perpendicular distance sampling object
                   izPerim = 'matrix',              #holds the iz perimeter in matrix form w/ hc
                   area = 'numeric',                #exact area of the inclusion zone
                   perimeter = 'SpatialPolygons',   #izone perimeter in 'SpatialPolygons' form
                   pgArea = 'numeric',              #polygon izone area approximation
                   pdsType = 'character'            #protocol method: see validity below
                  ),
    contains = 'downLogIZ',                         #a subclass of the 'downLogIZ' class
    validity = function(object) {
                 if(object@area < 0 || is.na(object@area))
                   return('object has negative or missing inclusion zone area!')

                 if(object@pgArea < 0 || is.na(object@pgArea))
                   return('object has negative or missing polygon area!')

                 if(!identical(object@pds@units, object@units))
                   return('object units do not match the perpendicularDistance object units')

                 if(!object@pdsType %in% .StemEnv$pdsTypes)
                   return('illegal pdsType requested!')

                 if(object@pds@units != object@downLog@units)
                   return('units for pds and downLog components incompatible!')
                 
                 return(TRUE)
               } #validity check
) #class perpendicularDistanceIZ 









#=================================================================================================
#
#  9. the omnibusPDSIZ class is a direct descendant of 'perpendicularDistanceIZ'...
#
#
#     added: 4-Feb-2011
#
setClass('omnibusPDSIZ',
    contains = 'perpendicularDistanceIZ'
) #class omnibusPDSIZ 






#=================================================================================================
#
#  10. the distanceLimitedIZ class is a direct descendant of 'downLogIZ'...
#
#
#     below: hc = homogeneous coordinates
#
#     added: 25-May-2011
#
setClass('distanceLimitedIZ',
    representation(dls = 'distanceLimited',         #the distanceLimited sampling object
                   izPerim = 'matrix',              #holds the iz perimeter in matrix form w/ hc
                   area = 'numeric',                #exact area of the inclusion zone
                   perimeter = 'SpatialPolygons',   #izone perimeter in 'SpatialPolygons' form
                   pgArea = 'numeric'               #polygon izone area approximation
                  ),
    contains = 'downLogIZ',                         #a subclass of the 'downLogIZ' class
    validity = function(object) {
                 if(object@area < 0 || is.na(object@area))
                   return('object has negative or missing inclusion zone area!')

                 if(object@pgArea < 0 || is.na(object@pgArea))
                   return('object has negative or missing polygon area!')
                 
                 return(TRUE)
               } #validity check
) #class distanceLimitedIZ 




#=================================================================================================
#
#  11. the distanceLimitedMCIZ class is a direct descendant of 'distanceLimitedIZ'; the dates
#      below are off compared with the parent because I actually developed the MC version first,
#      then added the distanceLimtedIZ class, this necessitated switching the defintion so
#      that the parent was the base for this class of objects
#
#     added: 22&23-Mar-2011
#
setClass('distanceLimitedMCIZ',
    contains = 'distanceLimitedIZ'
) #class distanceLimitedMCIZ 






#=================================================================================================
#
#  12. the distanceLimitedPDSIZ class is a direct descendant of 'perpendicularDistanceIZ'...
#
#      Two class unions are defined below to allow the slots to be either filled with something
#      of the correct class, or be empty with NULL
#
#     added: 8&10-Mar-2011
#
setClassUnion('pdsIZNull', c('perpendicularDistanceIZ', 'omnibusPDSIZ', 'NULL')) #allow for omnibus
setClassUnion('dlsIZNull', c('distanceLimitedIZ', 'NULL')) 
setClass('distanceLimitedPDSIZ',
    representation(dls = 'distanceLimited',               #the distanceLimited ArealSampling object
                   dlsDiameter = 'numeric',               #limiting diameter
                   pdsPart = 'pdsIZNull',                 #pds object for pds section of the log
                   dlsPart = 'dlsIZNull',                 #distance limited component of the log
                   pdsFull = 'pdsIZNull'                  #as if this were a full PDS object  (never NULL?)
                  ),
    contains = 'perpendicularDistanceIZ',
         
    validity = function(object) {

                 if(object@dlsDiameter <= 0)
                   return('distance limit diameter must be positive')
                 
                 return(TRUE)
               } #validity check
                  
) #class distanceLimitedPDSIZ 






#=================================================================================================
#
#  13. the omnibusDLPDSIZ class is a direct descendant of 'distanceLimitedPDSIZ'...
#
#      combines DL with omnibus PDS
#
#     added: 15-Mar-2011
#
setClass('omnibusDLPDSIZ',
    contains = 'distanceLimitedPDSIZ'
) #class omnibusDLPDSIZ 





#=================================================================================================
#
#  14. the hybridDLPDSIZ class is a direct descendant of 'distanceLimitedPDSIZ'...
#
#      combines DL with canonical PDS
#
#     added: 26-July-2011
#
setClass('hybridDLPDSIZ',
    contains = 'distanceLimitedPDSIZ'
) #class hybridDLPDSIZ 


         



         






#---------------------------------------------------------------------------
#
#   This section contains the InclusionZone class definitions for standing
#   trees...
#
#   Classes...
#
#     1. standingTreeIZ: virtual parent
#     2. circularPlotIZ: fixed-area circular plots
#     3. horizontalPointIZ: horizontal point sampling (6-Dec-2011)
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



#=================================================================================================
#
#  1. the standing tree inclusion zone class is just a direct descendant of 'InclusionZone'...
#
#
setClass('standingTreeIZ',
    representation(standingTree = 'standingTree'              #standingTree object
                  ),
    contains = c('InclusionZone', 'VIRTUAL'),             #a subclass of the virtual 'InclusionZone' class
    validity = function(object) {

                 if(!validObject(object@standingTree))
                   return('object in standingTree slot is an invalid object!')

                 if(!identical(object@standingTree@units, object@units))
                   return('measurement units do not match between standingTree and IZ objects!')

                 if(!identical(object@standingTree@spUnits, object@spUnits))
                   return('Spatial units do not match between standingTree and IZ objects!')
                 
                 return(TRUE)
               } #validity check
) #class standingTreeIZ 


 





#=================================================================================================
#
#  2. the circularPlotIZ class is just a direct descendant of 'standingTreeIZ'...
#
#
setClass('circularPlotIZ',
    representation(circularPlot = 'circularPlot'       #circularPlot object
                  ),
    contains = 'standingTreeIZ',                      #a subclass of the virtual 'standingTreeIZ' class
    validity = function(object) {

                 if(!validObject(object@circularPlot))
                   return('object in circularPlot slot is an invalid object!')

                 if(!identical(object@standingTree@units, object@circularPlot@units))
                   return('measurement units do not match between standingTree and circularPlot objects!')

                 if(!identical(object@standingTree@spUnits, object@circularPlot@spUnits))
                   return('Spatial units do not match between standingTree and circularPlot objects!')
                 
                 return(TRUE)
               } #validity check
) #class circularPlotIZ 






#=================================================================================================
#
#  3. the horizontalPointIZ class is a direct descendant of 'circularPlotIZ'...
#
#     what that means is that it saves coding further down the line as, for example, it will use
#     the methods for InclusionZoneGrid, Plot, etc. generics.
#
#
setClass('horizontalPointIZ',
    representation(angleGauge = 'angleGauge'         #angleGauge object
                   
                  ),
    contains = 'circularPlotIZ',                     #a subclass of the 'circularPlotIZ' class
    validity = function(object) {

                 if(!validObject(object@angleGauge))
                   return('object in angleGauge slot is an invalid object!')

                 if(!identical(object@standingTree@units, object@angleGauge@units))
                   return('measurement units do not match between standingTree and angleGauge objects!')
                 
                 return(TRUE)
               } #validity check
) #class horizontalPointIZ 
