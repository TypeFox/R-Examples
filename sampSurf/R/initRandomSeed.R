initRandomSeed = function(startSeed=NA, runQuiet=TRUE, ...)
{
#---------------------------------------------------------------------------
# 
#   This routine sets the random number "state", originally developed for
#   particle filtering.
#
#   The default (startSeed=NULL or NA) is to continue on with whatever stream
#   is currently defined in ".Random.seed". But if no seed is passed, and
#   .Random.seed does not exist, it is initialized to the current clock time.
#
#   Arguments...
#     startSeed = a single value, preferably integer; NULL or NA implies 
#                 default, see comments above
#     runQuiet = T: print some info; F: nothing
#     ... = gobbled
#
#   Returns...
#     the current seed either as newly set or continuing on (see comments below)
#
#   The default uniform RNG in R is "Mersenne-Twister" which produces a "seed
#   set" when initialized. This seed set seems to be initialized to postion
#   624 so that the next draw will change the seed set (and hence the position
#   to 1). Subsequently it will take 624 more draws to change the seed set,
#   to a new one &c. So what is returned by default from this function is
#   a vector of length 626.
#
#   Note: NULL added to NA as a default on 3-Aug-2010 (JHG).
#
#Author...	                                     Date: 18-Nov-2008
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   this will catch someone passing startSeed=.Random.seed, for example...
#
    if(any(!is.null(startSeed)) && any(!is.na(startSeed)) && length(startSeed)>1)
      stop(paste('Length of startSeed = ',length(startSeed),' must be 1!',sep=''))
       
#
#   use the seed passed if any...
#
    if(!is.null(startSeed) && !is.na(startSeed)) {
      set.seed(startSeed)
      if(!runQuiet)
        cat('\nSeed created from startSeed =',startSeed)
    }

#
#   if none exists, create a seed; otherwise, keep using what's there...
#
    if(!exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE)) { #create a seed if need be
      set.seed(as.numeric(Sys.time()))
      if(!runQuiet)
        cat('\nSeed created from system clock time.')
    }
    else if(!runQuiet && (is.null(startSeed) || is.na(startSeed)))
      cat('\nCurrent seed used.')

    if(!runQuiet) {
      rng = RNGkind()
      cat('\nCurrent RNG types: ',rng)
      if("Mersenne-Twister" %in% rng) {
        rng.pos = .Random.seed[2]
        cat('\nCurrent position in the seed set =', rng.pos)
        cat('\nWhich is:', .Random.seed[rng.pos+2]) #note, this is NOT the same as the seed itself!
        cat('\nNote: seed set changes next draw after the postion = 624')
      }
      cat('\n\n')
    }

    return(invisible(.Random.seed))
}   #initRandomSeed

      
