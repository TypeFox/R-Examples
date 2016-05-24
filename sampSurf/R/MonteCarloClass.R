#---------------------------------------------------------------------------
#
#   This file contains a range of definitions for implementing the class
#   structure for traditional Monte Carlo Methods for the estimation of
#   integrals (volume).
#
#   The classes contained here are...
#
#   A. "Stem" class...
#
#   1. "MonteCarloSampling" -- the base virtual class
#   2. "crudeMonteCarlo" -- as the name implies, inherits from (1)
#   3. "importanceSampling" -- inherits from (2)
#   4. "controlVariate" -- inherits from (3)
#   5. "antitheticSampling" -- the odd duck, combines (2)--(4)
#
#   B. "InclusionZone" class...
#
#     1. "MonteCarloSamplingIZ" -- virtual to be mixed with other IZ classes
#
#Author...									Date: 2-May-2013
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#


###########################################################################################
#
#   A. "Stem" class...
#
###########################################################################################




#=================================================================================================
#
# A.1. the MonteCarloSampling class...
#
#    note that this class is virtual, so it is the base class for all others...
#
#
setClass('MonteCarloSampling',
    representation(stem = 'Stem',                 #a "Stem" subclass object  ****Allow NULL or NA too???*********
                   segBnds = 'numeric',           #lower & upper height/length bounds for desired stem segment
                                                  #everything that follows is with respect to the above segment...
                   n.s = 'numeric',               #the number of Monte Carlo samples
                   startSeed = 'numeric',         #the start seed used in initRandomSeed, can be NA
                   u.s = 'numeric',               #the uniform random numbers used in sampling
#                  other...                   
                   description = 'character',     #descriptive comment
                   userArgs = 'list'              #dotted args passed, e.g., to proxy functions
                  ),
    contains = 'VIRTUAL',                         #note!
    prototype = list(stem = downLog(logLen=5),    #some defaults for validity checking
                     segBnds = c(0, 5),
                     startSeed = NA_real_,
                     u.s = runif(1),
                     n.s = 1,
                     description = "Monte Carlo Base",
                     userArgs = list()
                    ),         
    validity = function(object) {
                 if(is(object@stem, 'downLog'))
                   H = object@stem@logLen
                 else
                   H = object@stem@height
                 #segment bounds checks...
                 segBnds = object@segBnds
                 if(!length(segBnds) == 2)
                   return('Illegal segBnds length, must be of length 2!')
                 if(segBnds[2] <= segBnds[1] || segBnds[1] < 0 ||
                    segBnds[1] > H || segBnds[2] > H)
                   return('Illegal segment heights/lengths specified in segBnds slot!')

                 n.s = object@n.s
                 u.s = object@u.s
                 if(n.s <= 0)
                   return(paste('Number of Monte Carlo samples=',n.s,'!'))
                 if(length(n.s) > 1)
                   return('Number of Monte Carlo samples must be a scalar!')
                 if(n.s%%1 > 0)
                   return('Number of Monte Carlo samples must be a positive whole number!')
                 if(length(u.s) != n.s )
                   return(paste('Random number vector not of length',n.s))
                 
                 return(TRUE)
               } #validity check
) #class MonteCarloSampling 




#=================================================================================================
#
# A.2. the crudeMonteCarlo class...
#
#    note that this class is a direct descendent of the "MonteCarloSampling" class...
#
#
setClass('crudeMonteCarlo',
    representation(proxy = 'character',           #the proxy taper function used--always uniform here!
                   diam.s = 'numeric',            #sampled diameters at hgt.s
                   rho.s = 'numeric',             #corresponding cross-sectional area
                   hgt.s = 'numeric',             #sampled heights at s points
                   vol.s = 'numeric',             #corresponding volume estimates
                   volEst = 'numeric',            #mean volume for the stem
                   volVar = 'numeric',            #volume variance estimate within the stem
                   ci.lo = 'numeric',             #lower confidence interval for estimate within stem
                   ci.up = 'numeric',             #upper confidence interval for estimate within stem
                   alphaLevel = 'numeric',        #two-tailed alpha level for CIs
                   trueVol = 'numeric',           #true stem volume within segBnds
                   relErrPct = 'numeric'          #relative error in percent
                  ),
    contains = 'MonteCarloSampling',              #a subclass of the virtual base class
    #prototype = list(diam.s = NA_real_,
    #                 rho.s = NA_real_,
    #                 hgt.s = NA_real_,
    #                 vol.s = NA_real_,
    #                 volEst = NA_real_,
    #                 volVar = NA_real_,
    #                 ci.lo = NA_real_,
    #                 ci.up = NA_real_,
    #                 alphaLevel = 0.05,
    #                 trueVol = NA_real_,
    #                 relErrPct = NA_real_,
    #                ),
    validity = function(object) {

                 n.s = object@n.s
                 if(length(object@diam.s) != n.s )
                   return(paste('Sampled diameter vector not of length',n.s))
                 if(length(object@rho.s) != n.s )
                   return(paste('Sampled cross-sectional area vector not of length',n.s))
                 if(length(object@hgt.s) != n.s )
                   return(paste('Sampled height/length vector not of length',n.s))
                 if(length(object@vol.s) != n.s )
                   return(paste('Sampled volume estimate vector not of length',n.s))
                 if(length(object@volEst) != 1)
                   return('Mean volume estimate not of length one!')
                 if(length(object@volEst) != 1)
                   return('Variance of volume estimate not of length one!')
                 if(length(object@ci.lo) != 1)
                   return('Lower confidence interval not of length one!')
                 if(length(object@ci.up) != 1)
                   return('Upper confidence interval not of length one!')
                 if(length(object@trueVol) != 1)
                   return('True segment volume not of length one!')
                 if(length(object@relErrPct) != 1)
                   return('Relative error percent not of length one!')
                 alphaLevel = object@alphaLevel
                 if(length(object@alphaLevel) != 1)
                   return('alpha level not of length one!')
                 if(alphaLevel >= 0.5 || alphaLevel <= 0)
                   return(paste('Illegal alpha probability level:',alphaLevel))
      
                 return(TRUE)
               } #validity check
) #class crudeMonteCarlo 





#=================================================================================================
#
# A.3. the importanceSampling class...
#
#    note that this class is a direct descendent of the "crudeMonteCarlo" class...
#
#
setClass('importanceSampling',
    #representation(),
    contains = 'crudeMonteCarlo',
    validity = function(object) {
     
                 return(TRUE)
               } #validity check
) #class importanceSampling 




#=================================================================================================
#
# A.4. the controlVariate class...
#
#    note that this class is a direct descendent of the "importanceSampling" class...
#
#
setClass('controlVariate',
    representation(diff.s = 'numeric'        #the CV differences 
                  ),
    contains = 'importanceSampling',
    validity = function(object) {
      
                 n.s = object@n.s
                 if(length(object@diff.s) != n.s )
                   return(paste('Control variate difference vector not of length',n.s))
     
                 return(TRUE)
               } #validity check
) #class controlVariate






#=================================================================================================
#
# A.5. the antitheticSampling class...
#
#
setClass('antitheticSampling',
    representation(mcsObj = 'MonteCarloSampling',    #subclass object
                   mcsAnti = 'MonteCarloSampling',   #subclass object antithetic sample
                   volEst = 'numeric',               #mean volume for the stem
                   volVar = 'numeric',               #volume variance estimate within the stem
                   ci.lo = 'numeric',                #lower confidence interval for estimate within stem
                   ci.up = 'numeric',                #upper confidence interval for estimate within stem
                   alphaLevel = 'numeric',           #two-tailed alpha level for CIs
                   trueVol = 'numeric',              #true stem volume within segBnds
                   relErrPct = 'numeric',            #relative error in percent
                   description = 'character'
                  ),
    validity = function(object) {
      
     
                 return(TRUE)
               } #validity check
) #class antitheticSampling









###########################################################################################
#
#   B. "InclusionZone" class...
#
###########################################################################################



#=================================================================================================
#
#  B.1. the MonteCarloSamplingIZ class is a virtual class that can be mixed with one of the
#       normal subclasses of "InclusionZone" in a 'contains=' to get slots of both...
#
#
setClass('MonteCarloSamplingIZ',
    representation(mcsObj = 'MonteCarloSampling',           #dummy argument with all correct info for
                                                            #use in izGrid later
                   antithetic = 'logical',                  #will this be antithetic or not?
                   proxy = 'character'                      #proxy function name
                  ),
    contains = 'VIRTUAL', 
    validity = function(object) {

                 #mcsObj slots will be checked when we do the dummy call to the constructor in
                 #some subclass of 'MonteCarloSamplingIZ' object's constructor
      
                 #could be logical NA, but that's not allowed!...
                 if(is.na(object@antithetic))
                   return('antithetic slot must be TRUE or FALSE!')
                 #similar...
                 if(is.na(object@proxy))
                   return('aproxy slot must contain the name of a proxy function!')

                 return(TRUE)
               } #validity check
) #class MonteCarloSamplingIZ 
