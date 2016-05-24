#---------------------------------------------------------------------------
#
#   This file holds the S4 definition for the constructor methods of the
#   monte class & subclasses and related classes...
#
#   The methods include...
#     1. "montePop" for general numeric vector
#     2. "montePop" for "sampSurf" objects
#     3. "monteNTSample" for normal theory samples for objects of class "montePop"
#     4. "monteBSSample" for bootstrap samples for objects of class "montePop"
#     5. "monte" for objects of class "montePop"
#     6. "monte" for objects of class "sampSurf"
#     7. "monte" for objects of class "numeric"
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
#   generic definition...
#
if(!isGeneric("montePop")) 
  setGeneric('montePop',  
             function(popVals, ...) standardGeneric('montePop'),
             signature = c('popVals')
            )


if(!isGeneric("monteNTSample")) 
  setGeneric('monteNTSample',  
             function(population, ...) standardGeneric('monteNTSample'),
             signature = c('population')
            )


if(!isGeneric("monteBSSample")) 
  setGeneric('monteBSSample',  
             function(population, ...) standardGeneric('monteBSSample'),
             signature = c('population')
            )


if(!isGeneric("monte")) 
  setGeneric('monte',  
             function(object, ...) standardGeneric('monte'),
             signature = c('object')
            )








#================================================================================
#  1. method for class montePop--numeric vector population...
#
setMethod('montePop',
          signature(popVals = 'numeric'),
function(popVals,
         zeroTruncate = FALSE,
         n = NA,
         description = 'Monte Carlo Population Object',
         ...
        )
{
#------------------------------------------------------------------------------
#
#   well we really should have a valid population here to begin with...
#
    if(length(popVals) <= 1)
      stop('Need a population of size > 1!')

#
#   no need for NAs in the population, this will just make things more difficult...
#
    if(any(is.na(popVals))) {
      popVals = popVals[!is.na(popVals)]
      cat('\n Missing values (NA) have been removed from the population\n')
    }
   
    if(zeroTruncate)
      popVals = popVals[popVals > 0]

#
#   population stats/parameters are simple enough...
#
    N = length(popVals)
    mean = mean(popVals)
    var = var(popVals)             #*(N-1)/N #alternative population variance
    stDev = sqrt(var)
    total = sum(popVals)

#
#   if sample size information is available, use it...
#
    if(!all(is.na(n))) {      
      nn = length(n)
      if(nn < 1 || !is.numeric(n))
        stop('Please specify \"n\" as an integer vector!')
      n = sort(round(n))                   #integer sample sizes only, sort drops NAs by default
      n.names = paste('n',n,sep='.')
      names(n) = n.names
      if(any(n >= N))
        stop('Sample sizes \"n\" can not be larger than the population size \"N\"!')
      fpc = (N - n)/N
      names(fpc) = n.names
      varMean = var/n * fpc
      stErr = sqrt(varMean)
    }
    else {                                #otherwise, it's okay for these slots to be NA
      n = NA_real_
      fpc = NA_real_
      varMean = fpc
      stErr = fpc
    }

    pop = new('montePop',
              mean = mean,
              var = var,
              stDev = stDev,
              N = N,
              total = total,
              popVals = popVals,
              description = description,
              zeroTruncated = zeroTruncate,
#
              n = n,
              fpc = fpc,
              varMean = varMean,
              stErr = stErr
             )
    
    return(pop)
    
}   #montePop constructor for 'numeric'
)   #setMethod






#================================================================================
#  2. method for class montePop--sampSurf signature...
#
setMethod('montePop',
          signature(popVals = 'sampSurf'),
function(popVals,
         zeroTruncate = TRUE,
         n = NA,
         description = 'Monte Carlo Population Object: sampSurf',
         ...
        )
{
#------------------------------------------------------------------------------
#
#   first, extract the population values==the object@tract grid values...
#
    vals = getValues(popVals@tract)

    pop = montePop(vals, zeroTruncate = zeroTruncate, n = n, description = description, ...)

    return(pop)

}   #montePop constructor for 'sampSurf'
)   #setMethod










#================================================================================
#  3. method for functions and class monteNTSample...
#
# 
#
#
setMethod('monteNTSample',
          signature(population = 'montePop'),
function(population,
         n = c(10),
         mcSamples = 100,   #number of MC samples
         alpha = 0.05,      #two-tailed alpha level
         replace = TRUE,
         startSeed = NA,
         runQuiet = TRUE,
         ...
        )
{
#------------------------------------------------------------------------------
#
#   a few checks...
#
    nn = length(n)
    if(nn < 1 || !is.numeric(n))
      stop('Please specify \"n\" as an integer vector!')
    n = sort(round(n))                   #integer sample sizes only, sort drops NAs by default
    n.names = paste('n',n,sep='.')
    names(n) = n.names
    if(mcSamples < 2 || !is.numeric(mcSamples) || length(mcSamples)>1) #may eventually allow vector
      stop('Please specify \"mcSamples\" as an integer scalar >= 2!')

    if(alpha <=0 || alpha >=1)
      stop('Please specify 0<alpha<1!')
    
    
#
#   get population info...
#
    popVals = population@popVals
    popMean = population@mean
    N = population@N
    if(any(n >= N))
      stop('Sample sizes \"n\" can not be larger than the population size \"N\"!')
    fpc = (N - n)/N
    names(fpc) = n.names
    
    
#
#   first we draw the samples and compute the stats on each--stored in data frames...
#
    if(!runQuiet)
      cat('\nCalculating Normal theory intervals...')
    ranSeed = initRandomSeed(startSeed)
    alpha = alpha/2
    t.values = qt(1-alpha, n-1)
    names(t.values) = n.names
    means = data.frame(matrix(NA, nrow=mcSamples, ncol=nn))
    colnames(means) = n.names
    vars = means
    stDevs = means
    varMeans = means
    stErrs = means
    lowerCIs = means
    upperCIs = means
    caught = means
    caught[] = FALSE
    for(i in seq_len(mcSamples)) {
      for(j in seq_len(nn)) {
        n.j = n[j]
        samp = sample(popVals, n.j, replace=replace, ...)
        means[i, j] = mean(samp)
        vars[i, j] = var(samp)
        stDevs[i, j] = sqrt(vars[i, j])
        varMeans[i, j] = vars[i, j]/n.j * fpc[j]
        stErrs[i, j] = sqrt(varMeans[i, j])
        lowerCIs[i, j] = means[i, j] - t.values[j]*stErrs[i, j]
        upperCIs[i, j] = means[i, j] + t.values[j]*stErrs[i, j]
        if(lowerCIs[i, j] <= popMean && popMean <= upperCIs[i, j])
          caught[i, j] = TRUE
      }
    }
    caughtPct = colMeans(caught) * 100
    names(caughtPct) = n.names

#
#   normal theory summary statistics--grand means over each sample for each n...
#
    nt.means = colMeans(means)
    nt.vars = colMeans(vars)
    nt.stDevs = colMeans(stDevs)
    nt.varMeans = colMeans(varMeans)
    nt.stErrs = colMeans(stErrs)
    nt.lowerCIs = colMeans(lowerCIs)
    nt.upperCIs = colMeans(upperCIs)
    nt.Stats = data.frame(rbind(nt.means, nt.vars, nt.stDevs, nt.varMeans, nt.stErrs, nt.lowerCIs, nt.upperCIs))
    rownames(nt.Stats) = c('mean', 'var', 'stDev', 'VarMean', 'stErr', 'lowerCI', 'upperCI')
    colnames(nt.Stats) = n.names

    if(!runQuiet) cat('\n')   

    ms = new('monteNTSample', 
             mcSamples = mcSamples,
             n = n,
             fpc = fpc,
             means = means,
             vars = vars,
             stDevs = stDevs,
             varMeans = varMeans,
             stErrs = stErrs,
             lowerCIs = lowerCIs,
             upperCIs = upperCIs,
             caught = caught,
             caughtPct = caughtPct,
             stats = nt.Stats,
             alpha = alpha*2,       #two-tailed
             t.values = t.values,
             replace = replace,
             ranSeed = ranSeed
            )

    return(ms)

}   #monteNTSample constructor
)   #setMethod








#================================================================================
#  4. method for functions and class monteBSSample...
#
#  note that the way this is set up currently, the samples drawn here are not
#  the same as those drawn in monteNTsample even with starting the RNG at the
#  right seed because of the bootstrap iterations w/in the loop. One way to
#  change this would be to draw the samples and save them first from the start
#  seed, then use them in an independent loop for bootstrapping.
#
#
setMethod('monteBSSample',
          signature(population = 'montePop'),
function(population,
         n = c(10),
         mcSamples = 100,   #number of MC samples
         R = 100,           #number of bootstrap samples/replicates
         alpha = 0.05,      #two-tailed alpha level
         replace = TRUE,
         startSeed = NA,
         runQuiet = TRUE,
         ...
        )
{
#------------------------------------------------------------------------------
#
#   a few checks...
#
    nn = length(n)
    if(nn < 1 || !is.numeric(n))
      stop('Please specify \"n\" as an integer vector!')
    n = sort(round(n))                   #integer sample sizes only, sort drops NAs by default
    n.names = paste('n',n,sep='.')
    names(n) = n.names
    mcSamples = round(mcSamples)
    if(mcSamples < 2 || !is.numeric(mcSamples) || length(mcSamples)>1) #may eventually allow vector
      stop('Please specify \"mcSamples\" as an integer scalar >= 2!')
    R = round(R)
    if(R < 10 || !is.numeric(R) || length(R)>1) 
      stop('Please specify \"R\" as an integer scalar >= 10!')

    if(alpha <=0 || alpha >=1)
      stop('Please specify 0<alpha<1!')

    
#
#   get population info...
#
    popVals = population@popVals
    popMean = population@mean
    N = population@N
    if(any(n >= N))
      stop('Sample sizes \"n\" can not be larger than the population size \"N\"!')
    fpc = (N - n)/N
    names(fpc) = n.names
    
    
#
#   first we draw the samples and compute the stats on each--stored in data frames...
#
    if(!runQuiet)
      cat('\nCalculating bootstrap intervals...')
    ranSeed = initRandomSeed(startSeed)
    means = data.frame(matrix(NA, nrow=mcSamples, ncol=nn))
    colnames(means) = n.names
    vars = means
    stDevs = means
    varMeans = means
    stErrs = means
    lowerCIs = means
    upperCIs = means
    caught = means
    caught[] = FALSE
    meanboot = function(x,idx) mean(x[idx])  #calculate the mean for each BS sample
    degenerate = rep(0, nn)                  #keep tract of number of degenerate samples
    names(degenerate) = n.names
    for(i in seq_len(mcSamples)) {
      for(j in seq_len(nn)) {
        n.j = n[j]
        samp = sample(popVals, n.j, replace=replace) 
        boot.samp = boot(samp, meanboot, R, ...)
        means[i, j] = mean(boot.samp$t0)
        vars[i, j] = var(samp)                #note: original sample--not bootstrap...
        stDevs[i, j] = sqrt(vars[i, j])
        varMeans[i, j] = vars[i, j]/n.j * fpc[j]
        stErrs[i, j] = sqrt(varMeans[i, j])
        if(length(unique(samp)) == 1) {   
          degenerate[j] = degenerate[j] + 1
          next
        }
        suppressWarnings({                    #can be noisy with warnings...
                         boot.cis = boot.ci(boot.samp, 1-alpha, type='bca')
                         })
        lowerCIs[i, j] = boot.cis$bca[1, 4]
        upperCIs[i, j] = boot.cis$bca[1, 5]
        if(lowerCIs[i, j] <= popMean && popMean <= upperCIs[i, j])
          caught[i, j] = TRUE
      }
    }
    caughtPct = colMeans(caught) * 100
    names(caughtPct) = n.names

#
#   bootstrap summary statistics--grand means over each sample for each n;
#   we can have NAs in bootstrap CIs if there were degenerate samples...
#
    m.means = colMeans(means, na.rm = TRUE)        #remove NAs just in case here too...
    m.vars = colMeans(vars, na.rm = TRUE)
    m.stDevs = colMeans(stDevs, na.rm = TRUE)
    m.varMeans = colMeans(varMeans, na.rm = TRUE)
    m.stErrs = colMeans(stErrs, na.rm = TRUE)
    m.lowerCIs = colMeans(lowerCIs, na.rm = TRUE)
    m.upperCIs = colMeans(upperCIs, na.rm = TRUE)
    bs.Stats = data.frame(rbind(m.means, m.vars, m.stDevs, m.varMeans, m.stErrs, m.lowerCIs, m.upperCIs))
    rownames(bs.Stats) = c('mean', 'var', 'stDev', 'varMean', 'stErr', 'lowerCI', 'upperCI')
    colnames(bs.Stats) = n.names

    if(!runQuiet) cat('\n')   

    ms = new('monteBSSample', 
             mcSamples = mcSamples,
             n = n,
             fpc = fpc,
             means = means,
             vars = vars,
             stDevs = stDevs,
             varMeans = varMeans,
             stErrs = stErrs,
             lowerCIs = lowerCIs,
             upperCIs = upperCIs,
             caught = caught,
             caughtPct = caughtPct,
             stats = bs.Stats,
             R = R,
             alpha = alpha,
             replace = replace,
             degenerate = degenerate,
             ranSeed = ranSeed
            )

    return(ms)

}   #monteBSSample constructor
)   #setMethod











#================================================================================
#  5. constructor method for class monte from montePop...
#
setMethod('monte',
          signature(object = 'montePop'),
function(object,
         zeroTruncate = TRUE,
         n = object@n,                                      #default to montePop values
         mcSamples = 100,                                   #number of MC samples
         type = c('both', 'normalTheory', 'bootstrap'),
         R = 100,                                           #number of bootstrap samples/replicates
         alpha = 0.05,                                      #two-tailed alpha level
         description = 'Monte Carlo Object',
         replace = TRUE,
         startSeed = NA,
         runQuiet = TRUE,
         ...
        )
{
#------------------------------------------------------------------------------
#
#   all other methods must call this one, so do the checks here...
#
#   a few checks-- sample size (n) checks are handled in montePop but n can be NA,
#                  there, so do it again here...
#
    nn = length(n)
    if(nn < 1 || !is.numeric(n))
      stop('Please specify \"n\" as an integer vector!')
    n = sort(round(n))                   #integer sample sizes only, sort drops NAs by default
 #   n.names = paste('n',n,sep='.')
 #   names(n) = n.names
    if(mcSamples < 1 || !is.numeric(mcSamples) || length(mcSamples)>1) #may eventually allow vector
      stop('Please specify \"mcSamples\" as an integer scalar!')

    if(alpha <=0 || alpha >=1)
      stop('Please specify 0<alpha<1!')

    type = match.arg(type)

#
#   check that sample sizes agree between objects...
    pop = object
    if(any(is.na(pop@n)))
      stop('Sample sizes must be present (no NAs) in \"montePop\" object for monte')
    if(length(intersect(pop@n, n)) != length(n))
      stop('Sample sizes between \"montePop\" object and argument \"n\" must match exactly!')
    
#
#   normal theory samples...
#
    if(type == 'both' || type == 'normalTheory')
      mnts = monteNTSample(pop, n=n, mcSamples=mcSamples,
                           alpha=alpha, replace=replace, startSeed=startSeed,
                           runQuiet=runQuiet, ...
                          )
    else
      mnts = NULL

#
#   bootstrap samples...
#    
    if(type == 'both' || type == 'bootstrap')
      mbss = monteBSSample(pop, n=n, mcSamples=mcSamples, R=R,
                           alpha=alpha, replace=replace, startSeed=startSeed,
                           runQuiet=runQuiet, ...
                          )
    else
      mbss = NULL
    
#
#   create the object...
#
    mo = new('monte',
             pop = pop,
             estimate = NA_character_,      #only meaningful for sampSurf objects
             NTsamples = mnts,
             BSsamples = mbss,
             description = description
            )

    return(mo)
}   #monte constructor for 'montePop' objects
)   #setMethod
    
      
    





#================================================================================
#  6. constructor method for class monte from sampSurf...
#
setMethod('monte',
          signature(object = 'sampSurf'),
function(object,
         zeroTruncate = TRUE,
         n = c(10),
         mcSamples = 100,                                   #number of MC samples
         type = c('both', 'normalTheory', 'bootstrap'),
         R = 100,                                           #number of bootstrap samples/replicates
         alpha = 0.05,                                      #two-tailed alpha level
         description = 'Monte Carlo Object',
         replace = TRUE,
         startSeed = NA,
         runQuiet = TRUE,
         ...
        )
{
#------------------------------------------------------------------------------
#
#   first, create the population object...
#
    pop = montePop(object, zeroTruncate = zeroTruncate, n = n, ...)

#
#   then just call the montePop constructor...
#
    mo = monte(object = pop,
               zeroTruncate = zeroTruncate,
               n = n,
               mcSamples = mcSamples,
               type = type,
               R = R,      
               alpha = alpha,
               description = description,
               replace = replace,
               startSeed = startSeed,
               runQuiet = runQuiet,
               ... )

    mo@estimate =  object@estimate

    return(mo)
}   #monte constructor for 'sampSurf' objects
)   #setMethod





#================================================================================
#  7. constructor method for class monte from numeric...
#
setMethod('monte',
          signature(object = 'numeric'),
function(object,
         zeroTruncate = TRUE,
         n = c(10),
         mcSamples = 100,                                   #number of MC samples
         type = c('both', 'normalTheory', 'bootstrap'),
         R = 100,                                           #number of bootstrap samples/replicates
         alpha = 0.05,                                      #two-tailed alpha level
         description = 'Monte Carlo Object',
         replace = TRUE,
         startSeed = NA,
         runQuiet = TRUE,
         ...
        )
{
#------------------------------------------------------------------------------
#
#   first, create the population object...
#
    pop = montePop(object, zeroTruncate = zeroTruncate, n = n, ...)

#
#   then just call the montePop constructor...
#
    mo = monte(object = pop,
               zeroTruncate = zeroTruncate,
               n = n,
               mcSamples = mcSamples,
               type = type,
               R = R,      
               alpha = alpha,
               description = description,
               replace = replace,
               startSeed = startSeed,
               runQuiet = runQuiet,
               ... )

    return(mo)
}   #monte constructor for 'numeric' objects
)   #setMethod
    
