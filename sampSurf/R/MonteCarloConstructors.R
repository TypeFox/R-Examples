#################################################################################
#
#   Generic with constructor method definitions for MonteCarloSampling subclasses
#   including...
#
#     A. 'crudeMonteCarlo' -- only 'cmcProxy' allowed to construct this class
#     B. 'importanceSampling' -- more general, any proxy
#     C. 'controlVariate'
#     D. 'antitheticSampling' -- composites one of the above
#
#   The actual constructors are...
#
#     1a. crudeMonteCarlo for "Stem" subclass objects
#     1b. crudeMonteCarlo for "StemContainer" subclass objects
#         Note that the above two are wrappers for 1a==2b&2c, 1b==2d&2e
#
#     2a. importanceSampling for "list" objects
#     ...Individual Stem objects...
#     2b. importanceSampling for "downLog" objects
#     2c. importanceSampling for "standingTree" objects
#     ...Collections of Stem objects...
#     2d. importanceSampling for "downLogs" objects
#     2e. importanceSampling for "standingTrees" objects
#
#     3a. controlVariate  for "Stem" subclass objects
#     3b. controlVariate for "StemContainer" subclass objects
#         Note that the above two are wrappers for 4a==2b&2c, 4b==2d&2e
#
#     4a. antitheticSampling for any of the above "MonteCarloSampling" subclasses
#     4b. antitheticSampling for "mcsContainer" subclass objects
#
#              -----------------------------------------
#
#   A note on the use of the n.s, startSeed and u.s arguments in the
#   constructor functions...
#
#     1. single stems...
#
#        a. if any u.s are NA or NULL, then n.s and startSeed are used
#        b. else u.s is used, and n.s is set based on its length
#
#     2. collections of stems...
#
#        a. if any u.s are NA or NULL, then n.s and startSeed are used--
#           ***>this results in DIFFERENT u.s for each stem in the collection
#        b. else u.s is used, and n.s is set based on its length within
#           the sigle stem method--
#           ***>this results in the SAME u.s for each stem in the collection
#
#        In both cases, startSeed is set to NA in the single stem method,
#        so the current stream (new (a) or u.s (b)) is used regardless.
#
#     Please note in both 1. or 2. that either n.s & startSeed are used as a
#     pair to specify what happens with u.s ignored; or u.s is used with
#     n.s & startSeed ignored.
#
#
#Author...									Date: 2-May-2013 & beyond...
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#


#================================================================================
#
#   generic definitions...
#
if(!isGeneric("crudeMonteCarlo")) 
  setGeneric('crudeMonteCarlo',  
             function(object, ...) standardGeneric('crudeMonteCarlo'),
             signature = c('object')
            )

if(!isGeneric("importanceSampling")) 
  setGeneric('importanceSampling',  
             function(object, ...) standardGeneric('importanceSampling'),
             signature = c('object')
            )

if(!isGeneric("controlVariate")) 
  setGeneric('controlVariate',  
             function(object, ...) standardGeneric('controlVariate'),
             signature = c('object')
            )


if(!isGeneric("antitheticSampling")) 
  setGeneric('antitheticSampling',  
             function(object, ...) standardGeneric('antitheticSampling'),
             signature = c('object')
            )


      


#################################################################################
#
#                 Crude Monte Carlo Sampling section...
#
#################################################################################

       
#================================================================================
#  1a. crudeMonteCarlo constructor method for 'Stem' subclasses...
#
setMethod('crudeMonteCarlo',
          signature(object = 'Stem'),
          
function(object,
         segBnds = NULL,                          #segment length bounds
         n.s = 1,                                 #no of samples
         startSeed = NA,
         u.s = NA,
         alphaLevel = 0.05,                       #two-tailed alpha for CIs
         description = 'crude Monte Carlo',
         ...
        )
{
#---------------------------------------------------------------------------
#
#   note that this simply calls the IS constructor with CMC proxy...
#
    cmcObj = importanceSampling(object = object,
                                segBnds = segBnds,
                                n.s = n.s,
                                startSeed = startSeed,
                                u.s = u.s,
                                proxy = 'cmcProxy',      #always use this proxy for CMC
                                alphaLevel = alphaLevel,
                                description = description,
                                ...)
    
    return( cmcObj )
}   #crudeMonteCarlo for "Stem"
)   #setMethod



#================================================================================
#  1b. method for class 'StemContainer' and all subclasses; this will end up
#      calling either the 'downLogs' or 'standingTrees' method...
#
setMethod('crudeMonteCarlo',
          signature(object = 'StemContainer'),
          
function(object,
         segBnds = NULL,                          #segment length bounds
         n.s = 1,                                 #no of samples
         startSeed = NA,
         u.s = NA,
         alphaLevel = 0.05,                       #two-tailed alpha for CIs
         description = 'crude Monte Carlo',
         ...
        )
{
#---------------------------------------------------------------------------
#
#   note that this simply calls the IS constructor with CMC proxy...
#
    cmc = importanceSampling(object = object,
                             segBnds = segBnds,
                             n.s = n.s,
                             startSeed = startSeed,
                             u.s = u.s,
                             proxy = 'cmcProxy',      #always use this proxy for CMC
                             alphaLevel = alphaLevel,
                             description = description,
                             ...)    
    return( cmc )
}   #crudeMonteCarlo for "StemContainer"
)   #setMethod








#################################################################################
#
#                 Importance Sampling section...
#
#################################################################################




#================================================================================
#  2a. method for class 'list' -- please do not use this method, it is the "guts"
#      for the other two wrapper methods which call it and are more convenient...
#
#  Note that this will create a 'controlVariate' object if the argument below is
#  TRUE, otherwise it will create either a 'crudeMonteCarlo' or 'importanceSampling'
#  object, depending on what proxy is used. This is a small stretch from exact
#  oop, but since CMC is a special case of IS, it was simpler to do things
#  this way--and CV is just a slight twist on IS...
#
setMethod('importanceSampling',
          signature(object = 'list'),
          
function(object,
         segBnds=c(low=0, up=object$height),       #height/length bounds
         n.s = 1,                                  #no of imp samples
         startSeed = NA,
         u.s = NA,
         proxy = 'gvProxy',
         alphaLevel = 0.05,                        #two-tailed alpha for CIs
         controlVariate = FALSE,                   #TRUE for CV method
         description = 'Monte Carlo Sampling',
         ...
        )
{
#---------------------------------------------------------------------------
#
#   deconstruct the list passed, otherwise everything is the same for
#   standingTrees or downLogs...
#
#---------------------------------------------------------------------------
#
    height = object$height
    taper = object$taper
    units = object$units
    stemVol = object$stemVol
    stem = object$stem

#
#   a couple checks...
#
    if(any(is.null(segBnds)) || any(is.na(segBnds))) {
      segBnds = rep(NA,2)
      segBnds[1] = 0
      segBnds[2] = height
    }
    if(segBnds[1] < 0 || segBnds > height || segBnds[2] <= segBnds[1])
      stop('Illegal height bounds on stem in importance sampling!')
    n.s = as.integer(round(n.s))                   #integer to be consistent with length() below
    if(n.s < 1)
      stop('You must specify an integer number of importance sampling points!')


    #cross-sectional area (ba) conversion factor and height to dbh...
    csaFactor =  ifelse(units==.StemEnv$msrUnits$English, .StemEnv$baFactor['English'],
                       .StemEnv$baFactor['metric'])

#
#   get the random number(s) used in IS for s=1,...n.s points; two options...
#     (1) use n.s and startSeed to specify the random number stream
#     (2) use the u.s stream passed as an argument, this will be useful
#         in multiple calls for the same stem under antithetic sampling
#
    if(any(is.null(u.s)) || any(is.na(u.s))) {    #use n.s and startSeed
      initRandomSeed(startSeed)
      u.s = runif(n.s)
    }
    else                                          #useful for antithetic re-call of method
      n.s = as.integer(length(u.s))               #use u.s, which determines n.s 


#
#   use the desired proxy and get the importance heights...
#
    proxyFun = getProxy(proxy)                            #retrieve the actual proxy function
    p = proxyFun(stem, u.s, segBnds, ...)                 #and apply it
    g = p$g                                               #closure -- we need this for CMC, IS and CV
    G = p$G                                               #integral total--i.e., segment volume
    hgt.s = p$hgt.s
    
    #we also need segment range and heights from the CMC proxy if CV sampling...
    if(controlVariate) {                                  #draw heights using CMC
      cmcFun = getProxy('cmcProxy')                       #retrieve the CMC proxy
      pc = cmcFun(stem, u.s, segBnds, ...)                #and apply it
      finv = 1/pc$g()                                     #the segBnds difference/range
      hgt.s = pc$hgt.s                                    #CMC random heights for diameters and x-sec
    }

    #calculate the true measurements...
    diam.s = taperInterpolate(stem, 'diameter', hgt.s)    #vector of true diameters returned
    rho.s = csaFactor*diam.s*diam.s                       #cross-sectional area for diameters
    #and the estimates using the proxy returns...
    if(controlVariate) {                                  #CV
      diff.s = rho.s - g(hgt.s)                           #CV differences
      vol.s = G + finv*(diff.s)                           #CV volume estimates
    }
    else                                                  #CMC or IS
      vol.s = G/g(hgt.s) * rho.s                          #individual volume estimates
    vol.s = ifelse(is.nan(vol.s), NA, vol.s)              #NaN if zero cross-section at tip

#
#   some proxies could cause negative volumes (see e.g. the gvProxy)...
#
    if( any(vol.s[!is.na(vol.s)]<0) )
      warning('Negative volume estimates at one or more sample points have been detected!')

    
#
#   volume estimate, variance and CIs...
#
    volEst = mean(vol.s, na.rm=TRUE)
    signifLevel = 1-alphaLevel/2
    if(n.s > 1) {
      #volVar = mean( (G*rho.s/g(hgt.s) - volEst)^2, na.rm=TRUE )/(n.s-1) #note /n.s is in the mean()
      volVar = mean( (vol.s - volEst)^2, na.rm=TRUE )/(n.s-1) #note /n.s is in the mean()
      tv = qt(signifLevel, n.s-1)
      ci.half = tv*sqrt(volVar)
      ci.lo = volEst - ci.half
      ci.up = volEst + ci.half
    }
    else {
      volVar = NA_real_
      ci.lo = NA_real_
      ci.up = NA_real_
    }
    relErrPct = (volEst-stemVol)/stemVol*100

#
#   create the return object based on the proxy used...
#
    if(proxy == 'cmcProxy')
      objType = 'crudeMonteCarlo'
    else                                                #for IS and CV (for CV see below)
      objType = 'importanceSampling'

    iso = new(objType,
              stem = stem,
              segBnds = segBnds,
              u.s = u.s,
              n.s = n.s,
              description = description,
              userArgs = list(...),
              diam.s = diam.s,
              rho.s = rho.s,
              hgt.s = hgt.s,
              vol.s = vol.s,
              volEst = volEst,
              volVar = volVar,
              ci.lo = ci.lo,
              ci.up = ci.up,
              alphaLevel = alphaLevel,
              trueVol = stemVol,
              relErrPct = relErrPct,
              proxy = proxy,
              startSeed = ifelse(is.na(startSeed), as.double(NA), startSeed)
    )
    
#
#   coerce to 'controlVariate' from IS if applicable...
#
    if(controlVariate) {
      iso = as(iso, 'controlVariate')
      iso@diff.s = diff.s                              #add the differences
    }
    
    return(invisible(iso))        
}   #importanceSampling for "list"
)   #setMethod

      





#================================================================================
#  2b. method for class 'downLog'...
#
setMethod('importanceSampling',
          signature(object = 'downLog'),
          
function(object,
         segBnds=c(low=0,  up=object@logLen),     #segment length bounds
         n.s = 1,                                 #no of imp samples
         startSeed = NA,
         u.s = NA,
         proxy = 'gvProxy',
         alphaLevel = 0.05,                       #two-tailed alpha for CIs
         description = 'Importance Sampling',
         ...
        )
{
#---------------------------------------------------------------------------
#
#   just set the list up, and call the version for that signature...
#
    stem = list()
    stem$height = object@logLen
    stem$taper = object@taper
    stem$units = object@units
    stem$stemVol = object@logVol
    stem$stem = object                 #also require the object for interpolation
    

#    
#   handle the occasion where some bolt or section volume is desired rather than total...
#
    stem$stemVol = segmentVolume(object, segBnds)
    #segVol = segmentVolume(object, segBnds)
    #if(!is.na(segVol))
     # stem$stemVol = segVol
 
    isObj = importanceSampling(object = stem,
                               segBnds = segBnds,
                               n.s = n.s,
                               startSeed = startSeed,
                               u.s = u.s,
                               proxy = proxy,
                               alphaLevel = alphaLevel,
                               description = description,
                               ...)   
    return( isObj )
}   #importanceSampling for "downLog"
)   #setMethod
      



       
#================================================================================
#  2c. method for class 'standingTree'...
#
setMethod('importanceSampling',
          signature(object = 'standingTree'),
          
function(object,
         segBnds=c(low=0,  up=object@height),  #height bounds
         n.s = 1,                              #no of imp samples
         startSeed = NA,
         u.s = NA,
         proxy = 'gvProxy',
         alphaLevel = 0.05,                    #two-tailed alpha for CIs
         description = 'Importance Sampling',
         ...
        )
{
#---------------------------------------------------------------------------
#
#   just set the list up, and call the version for that signature...
#
    stem = list()
    stem$height = object@height
    stem$taper = object@taper
    stem$units = object@units
    stem$stemVol = object@treeVol
    stem$stem = object                 #also require the object for interpolation

#    
#   handle the occasion where some bolt or section volume is desired rather than total...
#
    stem$stemVol = segmentVolume(object, segBnds)
    #segVol = segmentVolume(object, segBnds)
    #if(!is.na(segVol))
     # stem$stemVol = segVol
 
    isObj = importanceSampling(object = stem,
                               segBnds = segBnds,
                               n.s = n.s,
                               startSeed = startSeed,
                               u.s = u.s,
                               proxy = proxy,
                               alphaLevel = alphaLevel,
                               description = description,
                               ...)   
    
    return( isObj )
}   #importanceSampling for "standingTree"
)   #setMethod
      




       
#================================================================================
#  2d. method for class 'downLogs'...
#
setMethod('importanceSampling',
          signature(object = 'downLogs'),
          
function(object,
         segBnds = NULL,
         n.s = 1,                             #no of imp samples
         startSeed = NA,
         u.s = NA,
         proxy = 'gvProxy',
         alphaLevel = 0.05,                    #two-tailed alpha for CIs
         description = 'Importance Sampling',
         ...
        )
{
#---------------------------------------------------------------------------
#
#   set the seed here for all trees, and keep the stream going...
#
    if(any(is.null(u.s)) || any(is.na(u.s)))   #set using startSeed
      initRandomSeed(startSeed)
    
    
    lis = lapply(object@logs, importanceSampling,
                 segBnds = segBnds,
                 n.s = n.s,
                 startSeed = NA,              #note
                 u.s = u.s,
                 proxy = proxy,
                 alphaLevel = alphaLevel,
                 description = description,
                 ...)

    cont = mcsContainer(lis)
    return(cont)    
}   #importanceSampling for "downLogs"
)   #setMethod
      
      



       
#================================================================================
#  2e. method for class 'standingTrees'...
#
setMethod('importanceSampling',
          signature(object = 'standingTrees'),
          
function(object,
         segBnds = NULL,
         n.s = 1,                             #no of imp samples
         startSeed = NA,
         u.s = NA,
         proxy = 'gvProxy',
         alphaLevel = 0.05,                    #two-tailed alpha for CIs
         description = 'Importance Sampling',
         ...
        )
{
#---------------------------------------------------------------------------
#
#   set the seed here for all trees, and keep the stream going...
#
    if(any(is.null(u.s)) || any(is.na(u.s)))   #set using startSeed
      initRandomSeed(startSeed)
    
    lis = lapply(object@trees, importanceSampling,
                 segBnds = segBnds,
                 n.s = n.s,
                 startSeed = NA,              #note
                 u.s = u.s,
                 proxy = proxy,
                 alphaLevel = alphaLevel,
                 description = description,
                 ...)

    cont = mcsContainer(lis)
    return(cont)
}   #importanceSampling for "standingTrees"
)   #setMethod


      







#################################################################################
#
#                 Control Variate Sampling section...
#
#################################################################################
      

       
#================================================================================
#  3a. controlVariate constructor method for 'Stem' subclasses...
#
setMethod('controlVariate',
          signature(object = 'Stem'),
          
function(object,
         segBnds = NULL,                          #segment length bounds
         n.s = 1,                                 #no of samples
         startSeed = NA,
         u.s = NA,
         proxy = 'gvProxy',
         alphaLevel = 0.05,                       #two-tailed alpha for CIs
         description = 'Control Variate Sampling',
         ...
        )
{
#---------------------------------------------------------------------------
#
#   note that this simply calls the IS constructor with appropriate CV=TRUE...
#
    cvObj = importanceSampling(object = object,
                               segBnds = segBnds,
                               n.s = n.s,
                               startSeed = startSeed,
                               u.s = u.s,
                               proxy = proxy,
                               alphaLevel = alphaLevel,
                               controlVariate = TRUE,
                               description = description,
                               ...)
    
    return( cvObj )
}   #controlVariate for "Stem"
)   #setMethod



#================================================================================
#  3b. method for class 'StemContainer' and all subclasses; this will end up
#      calling either the 'downLogs' or 'standingTrees' method...
#
setMethod('controlVariate',
          signature(object = 'StemContainer'),
          
function(object,
         segBnds = NULL,                          #segment length bounds
         n.s = 1,                                 #no of samples
         startSeed = NA,
         u.s = NA,
         proxy = 'gvProxy',
         alphaLevel = 0.05,                       #two-tailed alpha for CIs
         description = 'Control Variate Sampling',
         ...
        )
{
#---------------------------------------------------------------------------
#
#   note that this simply calls the IS constructor with appropriate CV=TRUE...
#
    cv = importanceSampling(object = object,
                            segBnds = segBnds,
                            n.s = n.s,
                            startSeed = startSeed,
                            u.s = u.s,
                            proxy = proxy,
                            alphaLevel = alphaLevel,
                            controlVariate = TRUE,
                            description = description,
                            ...)    
    return( cv )
}   #controlVariate for "StemContainer"
)   #setMethod









#################################################################################
#
#                 Antithetic Sampling section...
#
#################################################################################



#================================================================================
#
#   4a. Antithetic sampling on any of the MonteCarloSampling subclass objects...
#
setMethod('antitheticSampling',
          signature(object = 'MonteCarloSampling'),
          
function(object,
         alphaLevel = 0.05,                        #two-tailed alpha for CIs
         description = 'Antithetic Sampling',
         ...
        )
{
#---------------------------------------------------------------------------
#
#   n.s and startSeed are not needed in the call below, just the anithetic
#   (1-u.s) points...
#
#---------------------------------------------------------------------------
#
    n.s = object@n.s
    trueVol = object@trueVol
    
    if(is(object, 'controlVariate'))
      controlVariate = TRUE
    else
      controlVariate = FALSE

#
#   Note: there may be proxy-specific arguments in object@userArgs that need to be
#         passed on to the proxy for estimation to be correct (e.g., the solidTypeProxy
#         argument for wbProxy), the way to do this is to use do.call with a list of
#         arguments that include the object@userArgs list if it contains anything; this
#         slot will contain the "..." args from the original call...
#
    args = list(object = object@stem,
                segBnds = object@segBnds,
                #n.s = object@n.s,
                #startSeed = object@startSeed,
                u.s = 1 - object@u.s,          #antithetic points
                proxy = object@proxy,
                alphaLevel = alphaLevel,
                controlVariate = controlVariate,
                description = paste('antithetic: ', object@description, sep='')
                )
    if(length(object@userArgs)>0)           #any extras to be passed?
      args = c(args, object@userArgs)       #include them if so
   
    anti = do.call(importanceSampling, args)
    
#
#   volume estimate, variance and CIs...
#
    vol.s = (object@vol.s + anti@vol.s)/2
    volEst = (object@volEst + anti@volEst)/2
    signifLevel = 1-alphaLevel/2
    if(n.s > 1) {
      volVar = mean( (vol.s - volEst)^2, na.rm=TRUE )/(n.s-1) #note /n.s is in the mean()
      tv = qt(signifLevel, n.s-1)
      ci.half = tv*sqrt(volVar)
      ci.lo = volEst - ci.half
      ci.up = volEst + ci.half
    }
    else {
      volVar = NA_real_
      ci.lo = NA_real_
      ci.up = NA_real_
    }
    relErrPct = (volEst-trueVol)/trueVol*100

    antiObj = new('antitheticSampling',
                  mcsObj = object,
                  mcsAnti = anti,
                  volEst = volEst,
                  volVar = volVar,
                  ci.lo = ci.lo,
                  ci.up = ci.up,
                  alphaLevel = alphaLevel,
                  trueVol = trueVol,
                  relErrPct = relErrPct,
                  description = description
                 )
       
    return(invisible(antiObj))        
}   #antitheticSampling for "MonteCarloSampling" subclass
)   #setMethod






#================================================================================
#
#   4b. Antithetic sampling on mcsContainer objects...
#
setMethod('antitheticSampling',
          signature(object = 'mcsContainer'),
          
function(object,
         alphaLevel = 0.05,                        #two-tailed alpha for CIs
         description = 'Antithetic Sampling collection',
         ...
        )
{
#---------------------------------------------------------------------------
#
#   just apply the constructor to the list within the container and then
#   create the new container...
#    
    lis = lapply(object@mcsObjs, antitheticSampling,
                 alphaLevel = alphaLevel,
                 description = description,
                 ...)
    
    cont = antitheticContainer(lis)
    return(cont)
}   #antitheticSampling for "mscContainer"
)   #setMethod
