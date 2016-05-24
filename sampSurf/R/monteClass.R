#---------------------------------------------------------------------------
#
#
#Author...									Date: 16-Feb-2012
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
#  'montePop' contains the population of elements from which we will be drawing samples...
#
setClass('montePop',
         
#
#  slots for the class and its subclasses...
#
    representation(mean = 'numeric',                    #mu
                   var = 'numeric',                     #sigma^2
                   stDev = 'numeric',                   #sigma
                   N = 'numeric',                       #population size
                   total = 'numeric',                   #population total
                   popVals = 'numeric' ,                #vector population
                   zeroTruncated = 'logical',           #TRUE: zeros gone; FALSE: zeros present
                   n = 'numeric',                       #sample sizes for the below slots
                   fpc = 'numeric',                     #finite population correction by n
                   varMean = 'numeric',                 #variance of the mean by n
                   stErr = 'numeric',                   #standard error of the mean by n
                   description = 'character'            #more descriptive name
                  ),
    prototype = prototype(description = '',             #some defaults for validity checking--
                          mean = NA_real_,              #return an invalid object to force
                          var = NA_real_,               #the user to enter something good!
                          stDev = NA_real_,
                          N = NA_integer_,
                          total = NA_real_,
                          popVals = NA_real_,
                          zeroTruncated = FALSE,
                          fpc = NA_real_,        
                          varMean = NA_real_,
                          stErr = NA_real_,
                          n = NA_real_ 
                         ),
    #contains = 'VIRTUAL',
    validity = function(object) {
                 if(is.na(object@mean) || is.na(object@var) || is.na(object@stDev) ||
                    is.na(object@N) || is.na(object@total) || all(is.na(object@popVals)) )
                   return('Population object has missing slots--illegal object!')

                 #cast popVals in case the population has integer zeros...
                 if(object@zeroTruncated && any(sapply(as.double(object@popVals),identical,0.0)))
                   return('Population is supposed to be zero-truncated, but contains zeros!')

                 if(any(!is.na(object@n)) && any(object@n >= object@N))
                   return('At least on sample size \"n\" greater than population size \"N\'!')

#                check for valid names in object slots...                 
                 if(!all(is.na(object@n))) { 
                   n.x = 'slot must be in the form \"n.x\", where x is the respective sample size'
                   n.names = paste('n', object@n, sep='.')
                   if(!identical(names(object@n), n.names))
                     return(paste('Names for \"n\"',n.x)) 
                   if(!identical(names(object@fpc), n.names))
                     return(paste('Names for \"fpc\"', n.x))
                   if(!identical(names(object@varMean), n.names))
                     return(paste('Names for \"varMean\"', n.x)) 
                   if(!identical(names(object@stErr), n.names))
                     return(paste('Names for \"stErr\"', n.x))
                 }

                 return(TRUE)
               } #validity check
) #class montePop


#
#   some defaults for validity checking--return an invalid object to force
#   the user to enter something good!
#
#  ...or we could set this up to be a draw from some distribution in general...
#
setMethod('initialize', 'montePop',
          function(.Object, ...) {
             .Object@mean = NA_real_
             .Object@var = NA_real_
             .Object@stDev = NA_real_
             .Object@N = NA_integer_
             .Object@total = NA_real_
             .Object@popVals = NA_real_
             .Object@zeroTruncated = FALSE
             .Object@n = NA_real_
             .Object@fpc = NA_real_
             .Object@varMean = NA_real_
             .Object@stErr = NA_real_
             callNextMethod(.Object, ...)
}) #initialize
            
#=================================================================================================




         



#=================================================================================================
#
#  'monteSample' contains on set of MC samples for a given mcSample size and differing sample
#  sizes (n)...
#
#  This class is VIRTUAL, use on of the subclasses or make one of your own...
#
setClass('monteSample',
         
#
#  slots for the class and its subclasses...
#
    representation(mcSamples = 'numeric',               #number of Monte Carlo samples
                   n = 'numeric',                       #vector of sample sizes
                   alpha = 'numeric',                   #two-tailed alpha level
                   replace = 'logical',                 #TRUE: with replacement; FALSE: w/o
                   ranSeed = 'numeric',                 #initial random seed value
                                                        #following data frames dim: (mcSamples x length(n))
                   fpc = 'numeric',                     #finite population correction factor by sample size
                   means = 'data.frame',                #sample means by sample size 
                   vars = 'data.frame',                 #sample variances by sample size
                   stDevs = 'data.frame',               #sample standard deviations by sample size
                   varMeans = 'data.frame',             #sample variance of the mean by sample size
                   stErrs = 'data.frame',               #sample standard error(s) of the mean
                   lowerCIs = 'data.frame',             #normal theory lower CI
                   upperCIs = 'data.frame',             #normal theory upper CI
                   caught = 'data.frame',               #normal theory caught pop mean (TRUE/FALSE)
                   caughtPct = 'numeric',               #normal theory percent caught
                   stats = 'data.frame'                 #summary statistics
                  ),
         
    contains = 'VIRTUAL',
 
    validity = function(object) {
                 if( !all(apply(object@caught,2,is.logical)) )
                   return('monteSample: caught slot must be a data frame with all logical values!')

                 if(object@alpha <=0 || object@alpha >= 1)  #awfully wide range though
                   return(paste('monteSample: illegal alpha value (',object@alpha,') specified!',sep=''))

                 if(!all(object@fpc >= 0) || !all(object@fpc <= 1))
                   return('monteSample: all finite population corrections must be 0<=fpc<=1')

#                check for valid names in object slots...
                 n.x = 'slot must be in the form \"n.x\", where x is the respective sample size'
                 n.names = paste('n', object@n, sep='.')
                 if(!identical(names(object@n), n.names))
                   return(paste('Names for \"n\"', n.x )) 
                 if(!identical(names(object@fpc), n.names))
                   return(paste('Names for \"fpc\"', n.x))
                 if(!identical(colnames(object@means), n.names))
                   return(paste('Column names for \"means\"', n.x))
                 if(!identical(colnames(object@vars), n.names))
                   return(paste('Column names for \"vars\"', n.x))
                 if(!identical(colnames(object@stDevs), n.names))
                   return(paste('Column names for \"stDevs\"', n.x))
                 if(!identical(colnames(object@varMeans), n.names))
                   return(paste('Column names for \"varMeans\"', n.x))
                 if(!identical(colnames(object@stErrs), n.names))
                   return(paste('Column names for \"stErrs\"', n.x))
                 if(!identical(colnames(object@lowerCIs), n.names))
                   return(paste('Column names for \"lowerCIs\"', n.x))
                 if(!identical(colnames(object@upperCIs), n.names))
                   return(paste('Column names for \"upperCIs\"', n.x))
                 if(!identical(colnames(object@caught), n.names))
                   return(paste('Column names for \"caught\"', n.x))
                 if(!identical(names(object@caughtPct), n.names))
                   return(paste('Names for \"caughtPct\"', n.x )) 

                 return(TRUE)
               } #validity check
) #class monteSample
#=================================================================================================


#=================================================================================================
#
#  'monteNTSample': Normal Theory samples subclass...
#
setClass('monteNTSample',
         
#
#  slots for the class and its subclasses...
#
    representation(t.values = 'numeric'                #two-tailed t-value for each level n
                   ),
         
    contains = 'monteSample',
 
    validity = function(object) {

                 return(TRUE)
               } #validity check
) #class monteNTSample
#=================================================================================================



#=================================================================================================
#
#  'monteBSSample': Bootstrap samples subclass...
#
setClass('monteBSSample',
         
#
#  slots for the class and its subclasses...
#
    representation(degenerate = 'numeric',      #number of degenerate bs samples by sample size
                                                #can occur in small n on some attributes and sampling methods
                   R = 'numeric'                #number of bootstrap replicates
                   ),
         
    contains = 'monteSample',
 
    validity = function(object) {

                 return(TRUE)
               } #validity check
) #class monteBSSample
#=================================================================================================












setClassUnion('monteNTSampleOrNULL', c('monteNTSample', 'NULL'))
setClassUnion('monteBSSampleOrNULL', c('monteBSSample', 'NULL'))
#=================================================================================================
#
#  'monte' is the general class for CLT stuff...
#
setClass('monte',
         
#
#  slots for the class and its subclasses...
#
    representation(pop = 'montePop',                    #population object
                   estimate = 'character',              #attribute for sampling surface
                   NTsamples = 'monteNTSampleOrNULL',   #normal theory results
                   BSsamples = 'monteBSSampleOrNULL',   #bootstrap results
                   description = 'character'            #more descriptive name
                  ),
    prototype = prototype(description = ''        #some defaults for validity checking--
                    ),
    #contains = 'VIRTUAL',
    validity = function(object) {

                 #check for corresponding n here and in pop...
                 ns = 'Sample sizes (n) between \"montePop\" object and'
                 if(any(is.na(object@pop@n)))
                   return('Sample sizes (n) must be present with no NAs in \"montePop\" object for monte')
                 #also for NTsamples object, and no NAs...
                 if(!is.null(object@NTsamples)) {
                   if(any(is.na(object@NTsamples@n)))
                     return('Sample sizes (n) must be present with no NAs in \"NTsamples\" slot for monte')
                   if(length(intersect(object@pop@n, object@NTsamples@n)) != length(object@NTsamples@n))
                     return(paste(ns,'\"monteNTSample\" object must match exactly!'))
                 }
                 #also for BSsamples object, and no NAs...
                 if(!is.null(object@BSsamples)) {
                   if(any(is.na(object@BSsamples@n)))
                     return('Sample sizes (n) must be present with no NAs in \"BSsamples\" slot for monte')
                   if(length(intersect(object@pop@n, object@BSsamples@n)) != length(object@BSsamples@n))
                     return(paste(ns,'\"monteBSSample\" object must match exactly!'))
                 }
               

                 return(TRUE)
               } #validity check
) #class monte
#=================================================================================================
