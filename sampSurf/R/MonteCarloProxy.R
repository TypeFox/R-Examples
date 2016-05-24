#################################################################################
#
#   Proxy routines for Monte Carlo sampling...
#
#   Routines in this file include...
#
#     1. getProxy = find and retrieve a proxy function from the character name
#     2. cmcProxy = the crude Monte Carlo proxy
#     3. gvProxy = the simple proxy from Gregoire and Valentine (2008) p. 112
#     4. wbProxy = the default taper function in sampSurf as a proxy
#
#   All proxy functions must take the following arguments as the first three,
#   others are permitted following these; getProxy will check to make sure this
#   is correct...
#
#     stem = a "Stem" object
#     u.s = the uniform random numbers (vector) for MC sampling
#     segBnds = the lower & upper height/length bounds for sampling, normally
#               c(low=0, up=total hgt or length)
#
#   All proxy functions must return at least the following...
#
#   Returns...  (see notation in Gregoire and Valentine, p. 106)
#     g = the proxy function itself (as an R function/closure) always as a
#         function of hgt & returning either cross-sectional area or something
#         proportional to it when called
#     G = the integral of g == normalizing factor == segment volume
#     hgt.s = the sampled heights using the u.s argument vector passed--see
#             Gregoire & Valentine (2008) p. 109 for the root equation that
#             yields these heights under IS; for CMC (p. 98) and CV (p 115.)
#
#   The proxy functions will not be exported or readily available to the user
#   though they can be retrieved with getFromNamespace() if desired; nor is        ****????maybe
#   getProxy exported.
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



#================================================================================
#
#   1. This routine will find and return the desired proxy function from either:
#      one of the defaults in the namespace, or elsewhere, e.g., in the user's
#      workspace...
#
#   We could have used getAnywhere here, but I don't want to be looking in other
#   namespaces, etc. 
#
#   Note that the actual function itself is retrieved and returned...
#
getProxy = function(proxy = c('cmcProxy', 'gvProxy', 'wbProxy'),
                    ...)
{       
#--------------------------------------------------------------------------------
#
#   proxy = the possible proxy functions provided in sampSurf
#   
#--------------------------------------------------------------------------------
#
#   must pass a character, not default proxy...
#
    if(!is.character(proxy) || missing(proxy))
      stop('The proxy function must be specified as a character variable')
    
    ssProxies = eval(formals()$proxy)                   #just gives the choices above
    if(!is.na(match(proxy, ssProxies)))                 #want one of these?
      proxyFun = getFromNamespace(proxy, 'sampSurf')    #extract it
    else
      proxyFun = get(proxy)                             #must want something else

#
#   check it to make sure it has the correct arguments...
#
    proxyFormals = names(formals(proxy))
    if(proxyFormals[1] != 'stem' || proxyFormals[2] != 'u.s' || proxyFormals[3] != 'segBnds' ||
       is.na(match('...', proxyFormals)) ) #dots can be anywhere else, after other arguments
      stop(paste('The desired proxy (',proxy,') must have arguments: (stem, u.s, segBnds, ...)',sep=''))

#
#   and the correct returns...
#
    body = as.list(body(proxyFun))
    lb = length(body)
    retList = as.character(body[[lb]])[2]
    retList = gsub(' ', '', retList)
    allPresent = length(grep('g=', retList)) + length(grep('G=', retList)) + length(grep('hgt.s=', retList))
    if(allPresent != 3)
      stop(paste('The desired proxy (',proxy,') must return a list with components: g, G and hgt.s at a minimum!'))
   

    return(proxyFun)
}   #getProxy






       
#================================================================================
#
#   The proxy functions will not be exported or readily available to the user
#   though they can be retrieved with getFromNamespace() if desired...
#
#--------------------------------------------------------------------------------



#--------------------------------------------------------------------------------
#
#   2. The crude Monte Carlo proxy choses heights uniformly within the bole 
#      segment specified. Therefore, the hgt argument in g(), for example, is 
#      simply a dummy argument. The integral of g() is one over the bounds in
#      segBnds and is always one.
#
#--------------------------------------------------------------------------------
cmcProxy =  function(stem, u.s, segBnds, ...)
{
    g = function(hgt = NA) 1/(segBnds[2] - segBnds[1])
    G = 1
    hgt.s = segBnds[1] + u.s/g()
    return(list(g=g, G=G, hgt.s=hgt.s))
}  #cmcProxy



#--------------------------------------------------------------------------------
#
#   3. Proxy from Gregoire & Valentine (2008) p.112: (H-h)...
#
#      ***> This proxy can produce negative volumes under control variate sampling
#           because g() is only proportional to cross-sectional area, returning
#           a height that, when differenced, causes negative volumes (differences)
#
#--------------------------------------------------------------------------------
#
gvProxy = function(stem, u.s, segBnds, ...)
{
    if(is(stem, 'downLog'))
      height = stem@logLen
    else if(is(stem, 'standingTree'))
      height = stem@height
    else
      stop('Illegal \"Stem\" object passed!')
       
    g = function(hgt) height - hgt                 #cross-section proxy 
    G = (g(segBnds[1])^2 - g(segBnds[2])^2)/2      #total segment volume == normalizing constant
    #and the importance heights...
    hgt.s = height - sqrt( (1-u.s)*g(segBnds[1])^2 + u.s*g(segBnds[2])^2 ) #vector
    return(list(g=g, G=G, hgt.s=hgt.s))
}  #gvProxy



#--------------------------------------------------------------------------------
#
#   4. This proxy is developed from the default taper function used in sampSurf.
#
#   Please note that this routine makes use of R's closure in the g() function
#   passed back: g will always use the modified version of stem found here when
#   being evaluated (not the version of stem that might be in the calling routine
#   someplace else and passed in the "stem" argument), so it will always use the
#   proxy determined by the solidType settings below. This is the essence of 
#   closures, read up on them if it does not make sense (also lexical scoping).
#   (I don't mean to sound condescending here, there's just no room for explanation
#    and there is plenty published about them--see e.g., Chambers or Gentleman's
#    books.)
#
#   Arguments...
#
#     solidTypeProxy = (0,10] with (0,1) meaning that we take that percentage of
#                      the solidType for the stem passed as the proxy; [1,10] means
#                      use that value as the proxy, so any shape within the standard
#                      wbTaper function can be used for testing/comparison;
#                      NA means the proxy == the stem shape, unless it is user-defined
#                      taper, in which case the default proxy is always solidType=3
#
#     truncateProxyStem = TRUE: truncate the proxy top to avoid inflating estimates
#                               under IS--see below for details;
#                         FALSE: leave it alone, which is better for CV
#                         Note that this is only used if the stem tapers to essentially
#                         zero; otherwise, it left alone, as in FALSE
#                         This could be made numeric==the truncation diameter...
#
#     wbProxySolve = 'uniroot' for uniroot solution (default); 'nlminb' for solver--
#                    note that uniroot is about 5 times faster than nlminb.
#
#     warningsOn = TRUE: produce warnings; FALSE: silent
#
#   Please note with soldTypeProxy above that if the value specified equals the
#   actual value for the stem object, importance sampling will be exact, and
#   there will be no error in the estimate.
#
#   Proxy tip truncation...
#
#   Please also note that when the real tree tapers to the tip (zero diameter)
#   low solidTypeProxy values (less than 5) can cause sampling at the tip which
#   will generate either NaNs or inflated volume estimates. We fix this below
#   by making sure that the actual proxy has a top diameter of around an inch
#   (see below for details). In this case, then even when solidTypeProxy=NA, there
#   will be a difference between the true and proxy taper, so the volume extimate
#   will not be exact (because the stem tapers to the tip, but the proxy does not).
#
#   The use of warning()s below is unfortunate, but I think necessary so that
#   one is reminded of what is going on here. I will forget myself, so it is
#   as much for me. One can wrap the call to the MonteCarloSampling class
#   constructor in a suppressWarnings(...) or change options(warn) if desired.
#
#   29-May-2013: I decided to make warnings optional--see above--as they are a pain!
#
#--------------------------------------------------------------------------------
#
wbProxy = function(stem, u.s, segBnds,
                   solidTypeProxy = 3,                           #default is parabolic
                   truncateProxyStem = TRUE,
                   wbProxySolve = c('uniroot','nlminb'),
                   warningsOn = FALSE,
                   ...)
{
#
#   checks, etc...
#
    wbProxySolve = match.arg(wbProxySolve)
    if(!is.na(solidTypeProxy) && (solidTypeProxy<=0 || solidTypeProxy > .StemEnv$solidTypes[2]))
      stop(paste('Illegal value for solidTypeProxy (',solidTypeProxy,') in wbProxy!',sep=''))

    if(!is(stem, 'Stem'))
      stop('Illegal \"Stem\" object passed!')

#
#   define the shape of the proxy stem through the solidTypeProxy argument...
#
#***>Note that the stem object below is the proxy stem from here on...
#
#   A. If solidTypeProxy=NA, then...
#
#     a. If stem@solidType!=NULL, then leave alone, the proxy shape is the same
#        as the stem, so volume is exact.
#     b. If stem@solidType=NULL, then splines will be used in g() and G.
#        Just leave "stem" as is for this. This should be equivalent to using the 
#        exact solidType when known for the default proxy; i.e., the volume
#        estimate should be exact, except perhaps for a little rounding...
#
#   B. If solidTypeProxy!=NA, then...
#
#     a. if solidTypeProxy is in the legal range, use it
#        ...at this point solidTypeProxy can only be in (0,1)...
#     b. solidTypeProxy is less than one, and the stem uses the default taper model,
#        then take the proportion of stem@solidType
#     c. if the stem is user-supplied taper, then (0,1) does not make sense,
#        set to 3 as a compromise
#
    if(!is.na(solidTypeProxy)) {
      if(solidTypeProxy >= .StemEnv$solidTypes[1])                  #also applies to stem@solidType=NULL
        stem@solidType = solidTypeProxy               
      else if(!is.null(stem@solidType))                             #must be in (0,1)
        stem@solidType = max(1.0, stem@solidType*solidTypeProxy)    #can be made as close as desired to true taper
      else if(is.null(stem@solidType)) {                            #only possibilty left is stem@solidType=NULL
        if(warningsOn)
          warning(paste('User-defined taper and solidTypeProxy=',solidTypeProxy,'==> solidTypeProxy=3 used!'))
        stem@solidType = 3
      }
      else                                                          #should never get here
        stop(paste('Bad solidTypeProxy =',solidTypeProxy,'passed!'))
    }

#
#   we don't want this proxy to taper to zero diameter or we'll get inflated volumes
#   (see above), so set the top diameter to one inch in this case -- i.e., if it is
#   about a half inch, we call it zero as it will indeed inflate for diameters
#   less that this based on sample runs...
#
#   note that this is all that is required, we don't have to rebuild the stem object
#   (though we could) because taperInterpolate and segmentVolume will use the either
#   the default taper equation (solidType!=NULL) or the taper data (solidType=NULL)
#   which is set below. We could rebuild the stem to get the correct profile and
#   spatial components, but it is doubtful anyone would ever want this...
#
    units = stem@units
    #if(identical(stem@topDiam, 0)) {  #remember: in feet or meters!
    if(stem@topDiam < 0.01 && truncateProxyStem) {  #remember: in feet or meters!
      nt = nrow(stem@taper)
      penulDiam = stem@taper[nt-1, 'diameter']
      if(units == 'metric')
        truncDiam = 0.0254
      else
        truncDiam = 0.083333
      if(penulDiam > truncDiam)                                          #topDiam for solidType!=NULL 
        stem@topDiam = truncDiam
      else
        stem@topDiam = 0.5*penulDiam
      stem@taper[nt, 'diameter'] = stem@topDiam                          #for solidType=NULL stems
      if(warningsOn)
        warning(paste('Stem tapers to ~zero dib, proxy topDiam set to:',stem@topDiam,'to avoid inflated estimates!'))
      #re-build the stem to get graphics to be the same????...
    }
      
    csaFactor =  ifelse(units==.StemEnv$msrUnits$English, .StemEnv$baFactor['English'],
                        .StemEnv$baFactor['metric'])
      
    #cross-sectional area (closure, will always use the proxy stem object from above)... 
    g = function(hgt) csaFactor * taperInterpolate(stem, 'diameter', hgt)^2 
    G = segmentVolume(stem, segBnds)  #segment volume
      
    #function for uniroot...
    fh.uni = function(hgt, u.h) segmentVolume(stem, c(segBnds[1], hgt)) - u.h*G
    #function for nlminb...
    fh.nlm = function(hgt, u.h) {if(identical(segBnds[1], hgt)) #penalty--keep away from lower bound
                                   return(1e+10)                #which gives an error in segmentVolume
                                 ff = abs( segmentVolume(stem, c(segBnds[1], hgt)) - u.h*G )
                                 return(ff)
                                } #abs needed for nlminb

#
#   get the importance heights...
#
    n.u = length(u.s)
    hgt.s = rep(NA, n.u)
    for(i in seq_len(n.u)) {
      if(wbProxySolve == 'uniroot') {
        fh.lower = fh.uni(segBnds[1]+1e-5, u.s[i]) #note the offset from lower bound here is needed for segmentVolume
        fh.upper = fh.uni(segBnds[2], u.s[i])
        hgt.s[i] = uniroot(f=fh.uni, interval=segBnds, u.h=u.s[i], f.lower=fh.lower, f.upper=fh.upper)$root
      }
      else
        hgt.s[i] = nlminb(mean(segBnds), fh.nlm, lower=segBnds[1], upper=segBnds[2], u.h=u.s[i])$par
    }
    return(list(g=g, G=G, hgt.s=hgt.s))
} #wbProxy


