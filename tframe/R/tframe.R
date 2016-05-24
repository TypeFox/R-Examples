
classed <- function(x, cls) {class(x) <- cls; x}
# structure would work to replace classed (but adds some overhead).
#classed <- function(x, cls) structure(x, class=cls)

###########################################################################

#internal utility
# Use this with "for (i in seq(length=m) )" as m==0 returns NULL and for does no loops
seqN <- function(N) {if (0==length(N)) NULL else if (N<=0) NULL else seq(N)}



# start, end, frequency, time need to be masked from R base so that
# tframe methods can work on the tframe attribute rather than class(x)

# The .tframe methods are "default" methods for tframes. Other more specific
#  methods can be defined (see eg start.tstframe for tframes from ts objects). 


# Tobs should give the number of data points in the time direction.
# for consistency check Tobs needs to look at the data not the tframe,
# i.e. the number of (vector) observations.
Tobs <- function(x) UseMethod("Tobs")
Tobs.default <- function(x) {if (is.array(x)) dim(x)[1] else length(x) }

# the functions start, end, and frequency in tframe and dse do not 
#  need "...", but the generic in R has it, so it is added here.

start.tframed     <- function(x, ...) tfstart(tframe(x)) 
end.tframed       <- function(x, ...) tfend(tframe(x)) 
frequency.tframed <- function(x, ...) tffrequency(tframe(x)) 
time.tframed      <- function(x, ...) tftime(tframe(x)) 
Tobs.tframed   <- function(x)      Tobs(tframe(x)) 

start.tframe     <- function(x, ...) tfstart(x)
end.tframe       <- function(x, ...) tfend(x) 
frequency.tframe <- function(x, ...) tffrequency(x) 
time.tframe      <- function(x, ...) tftime(x) 
Tobs.tframe   <- function(x) # formerly default for tfTobs
  {if (is.null(x)) return(NULL) else
   if (!is.tframe(x)) x <- tframe(x)
   1+round((x[2]-x[1])*x[3])
  }

tfstart     <- function(x) UseMethod("tfstart")
tfend       <- function(x) UseMethod("tfend")
tffrequency <- function(x) UseMethod("tffrequency")
tftime <- function(x) UseMethod("tftime")


# these server two purposes. One is a method for tframe's. Two is a consistent
#programming method with tfstart(NULL) returning NULL (which start does not).
tfstart.default     <- function(x) 
  {if (is.null(x)) return(NULL) else
   #if (!is.tframe(x)) x <- tframe(x)
   #c(floor(x[1]), round(1 +(x[1]%%1)*x[3]))
   if (is.tframe(x)) c(floor(x[1]), round(1 +(x[1]%%1)*x[3]))
   else start(x)
  }
tfend.default       <- function(x)
  {if (is.null(x)) return(NULL) else
   #if (!is.tframe(x)) x <- tframe(x)
   #c(floor(x[2]), round(1 + (x[2]%%1)*x[3]))
   if (is.tframe(x)) c(floor(x[2]), round(1 + (x[2]%%1)*x[3]))
   else end(x)
  }
tffrequency.default <- function(x)
  {if (is.null(x)) return(NULL) else
   #if (!is.tframe(x)) x <- tframe(x)
   #x[3]
   if (is.tframe(x)) x[3] else frequency(x)
  }
tftime.default      <- function(x)
  {if (is.null(x)) return(NULL) else
   #if (!is.tframe(x)) x <- tframe(x)
   #tframed(x[1]+(seq(Tobs(x))-1)/x[3], tf=x)
   if (is.tframe(x)) tframed(x[1]+(seq(Tobs(x))-1)/x[3], tf=x)
   else time(x)
  }


tfL <- function(x, p=1) UseMethod("tfL") 
 
tfL.tframe <- function(x, p=1){ x + (p/x[3]) * c(1, 1, 0)}

tfL.default <- function(x, p=1){
    tframe(x) <- tfL(tframe(as.ts(x)), p=p)
    x} 


diff.tframed <- function(x, lag=1,   differences=1, ...)
   {tf <- diff(tframe(x), lag=lag, differences=differences) 
    tframe(x) <- NULL
    tframed(diff(x, lag=lag, differences=differences), tf)
    }

diff.tframe <- function (x,lag=1, differences=1, ...) 
 {d <- lag * differences
  tfTruncate(x, start=if(d >= 0) 1+d else NULL, 
                  end=if(d <  0) Tobs(x)-d else NULL)
 }


tfspan <- function(x, ...)
  {others <- list(...)
   tfspan <- x
   #this is a kludge to get the overall time span from the result of tbind.
   if (0 != length(others)) for (d in others) tfspan <- tbind(tfspan , d)
   tframe(tfspan)
  }


# Note tfprint prints the data. tframePrint  prints the tframe info. 

tfprint <- function(x, ...)  UseMethod("tfprint")

tfprint.default <- function(x,...)
 {xx <- x
  if(1 == nseries(xx)) xx <- matrix(xx, length(xx), 1)
  dimnames(xx) <- list(format(time(tframe(x))), seriesNames(x))
  tframe(xx) <- NULL
  seriesNames(xx) <- NULL
  print(xx, ...)
  invisible(x)
 }



tfwindow <- function(x, tf=NULL, start=tfstart(tf), end=tfend(tf), warn=TRUE)
  UseMethod("tfwindow")

tfwindow.default <- function(x, tf=NULL, start=tfstart(tf), end=tfend(tf), warn=TRUE)
  {# With the default warn=TRUE warnings will be issued if no truncation takes
   #  place because start or end is outside the range of data.
   # kludge
   x <- ts(x, start=tfstart(x), end=tfend(x), frequency=tffrequency(x))
   if (!warn) 
     {opts <- options(warn = -1)
      on.exit(options(opts))
     }
   y <- window(x, start=start, end=end)
   if (is.matrix(x) && !is.matrix(y) ) y <- matrix(y, length(y), ncol(x))
   y <- tframed(unclass(y), tframe(y))
   seriesNames(y) <- seriesNames(x)
   y
  }


# window a tframe
tfwindow.tframe <- function(x, tf=NULL, start=tfstart(tf), end=tfend(tf), warn=TRUE)
      tframe(tfwindow(time(x), tf=tf, start=start, end=end, warn=warn))

###############################################

#  tframe  methods   <<<<<<<<<<<<

################################################
is.tframe <- function(x) inherits(x, "tframe")
#is.tframed <- function(x) inherits(tframe(x), "tframe")
is.tframed <- function(x) UseMethod("is.tframed")
is.tframed.default <- function(x) {!is.null(tsp(x))}


# above does not distinguish "true" tframed objects since tframe(x) needs
# to try very hard to return tframes from ts and old tsp objects.
#is.Ttframed <- function(x) {!is.null(attr(x, "tframe"))}


tframe <- function(x) UseMethod("tframe")

tframe.default <- function(x){ #extract the tframe
  if(is.null(x)) NULL
  else if(is.tframe(x)) x   
  else if(!is.tframed(x)) NULL   
  #else if(is.tframed(x)) tframe(x)  this causes recursion. instead use
  else if (!is.null(attr(x, "tframe"))) attr(x, "tframe") # extractor
  else if (!is.null(tsp(x)))	classed(tsp(x), "tframe") # extractor
  else if(is.vector(x)) classed(c(1,length(x),1), "tframe") # extractor
  else if(is.matrix(x)) classed(c(1,  nrow(x),1), "tframe") # extractor
  else if(is.array(x) ) classed(c(1,dim(x)[1],1), "tframe") # extractor
  #else NULL
  # to catch possible lingering old representations
  else stop("Cannot extract tframe from tframed object")
}
# using classed(tsp(as.ts(x)), "tframe") in the last line above 
# makes too many things into tframes (eg lists)

as.tframed <- function(x) # guarantee x has a tframe 
 {# tframe(x) generates a default
  if (is.tframed(x)) x else tframed(x, tf=tframe(x))
 }
 
as.tframe <- function(...) #constructor
 {#expecting a combination of start, end, frequency, delta, Tobs,
  #which has enough info to calculate periods. (defaults work for other things.)
  # This is not very generic. The list must define a ts.
  lst <- list(...)
  if(is.null(lst$Tobs) & (is.null(lst$start) | is.null(lst$end)) ) 
     stop("Must supply Tobs or start and end.")

  if(is.null(lst$start) & is.null(lst$end)) lst$start <- c(1,1)
  if (! is.null(lst$frequency))   f <- lst$frequency
  else if (! is.null(lst$deltat)) f <- 1/lst$deltat
  else f <- 1
  #more generic date calc. would be nice here
  if (! is.null(lst$Tobs)) p <- lst$Tobs
  else p <- 1 + f * (lst$end[1] - lst$start[1]) + (lst$end[2] - lst$start[2])

  # ts seems to want missing values rather than null.
  if (is.null(lst$start)) z <- ts(rep(0,p), end=lst$end,   frequency=f) 
  else                    z <- ts(rep(0,p), start=lst$start, frequency=f)
  tframe(z)
  }
 

"tframe<-" <- function(x, value) 
  {if(is.null(value)) tfUnSet(x) else tfSet(value, x) 
  }

tfUnSet <- function(x) UseMethod("tfUnSet") # for NULL value
tfUnSet.default <- function(x) {
     # this is for unusual cases that actually are tframed class
     attr(x, "tframe") <- NULL
     class(x) <- class(x)[class(x) != "tframed"]
     # this is for old tsp cases
     tsp(x) <- NULL
     x
    } 

tfSet <- function(value, x) UseMethod("tfSet") # dispatch on value

# tfSet.default <- function(value, x) {
#   # It is tempting in the next to try and make a ts if value is from a ts, 
#   #  but that will not work for cases were x does not fit the ts model, so
#   #  that would break  tframe(x) <- tframe(y) 
#   if(!is.tframe(value)) {
#       # do.call does not seem to work when x is passed as NULL 
#       if(is.null(value$start) & is.null(value$end))
#                         stop("Could not determine a tframe from value.")
#       value <- as.tframe(start=value$start, end=value$end, 
#                         frequency=value$frequency, Tobs=Tobs(x))
#       }
#   if(! is.tframe(value)) stop("Could not determine a tframe from value.")
#   # next is checking after the fact, but value may just be start and freq
#   #  which is not enough to know Tobs
#   attr(x, "tframe") <- value
#   if((!is.null(value)) && !checktframeConsistent(tframe(x), x))
#      stop("time frame value in tframe assignment is not consistent with data.")
#   classed(x, c(class(x), "tframed"))
# }

  # It is tempting in the next to try and make a ts if value is from a ts, 
  #  but that will not work for cases were x does not fit the ts model, so
  #  that would break  tframe(x) <- tframe(y) 

tfSet.list <- function(value, x) {
  if(!is.tframe(value)) {
      # If value is not a tframe then only ts is attempted
      # do.call does not seem to work when x is passed as NULL 
      if(is.null(value$start) & is.null(value$end))
                        stop("Could not determine a tframe from value.")
      #value <- as.tframe(start=value$start, end=value$end, 
      #                  frequency=value$frequency, Tobs=Tobs(x))
      return(tfSet.tstframe(value, x))
      }
  #  if(! is.tframe(value)) stop("Could not determine a tframe from value.")
  #  # next is checking after the fact, but value may just be start and freq
  #  #  which is not enough to know Tobs
  #  attr(x, "tframe") <- value
  #  if((!is.null(value)) && !checktframeConsistent(tframe(x), x))
  #     stop("time frame value in tframe assignment is not consistent with data.")
  stop("stopped in tfSet.list")  #classed(x, c(class(x), "tframed"))
}

tfSet.default <- function(value, x) {
  if(  is.tframed(value)) return( tfSet(tframe(value), x)) #recall  
  if(is.numeric(value) && (length(value) == 3)) {
      # assuming tsp 
      tsp(x) <- value
      return(x)
      }
#  if(! is.tframe(value))  value <- as.tframe(value)  
#  if(! is.tframe(value)) stop("Could not determine a tframe from value.")
#  # next is checking after the fact, but value may just be start and freq
#  #  which is not enough to know Tobs
#  attr(x, "tframe") <- value
#  if((!is.null(value)) && !checktframeConsistent(tframe(x), x))
#     stop("time frame value in tframe assignment is not consistent with data.")
  stop("stopped in tfSet.default. Should never be here.")  #  classed(x, c(class(x), "tframed"))
}

# making tframed generic allows tframed.TSdata to specify input and output names

tframed <- function(x, tf=NULL, names = NULL, ...) UseMethod("tframed")

tframed.default <- function(x, tf=NULL, names = NULL, start=NULL, end=NULL, ...)
 {# return x as a tframed object with tframe tf
  if (!is.null(names))  seriesNames(x) <-  names
  if (is.null(tf))
     if ((!is.null(start)) | (!is.null(end))) 
           tf <- as.tframe(start=start, end=end, Tobs=Tobs(x), ...)
     else  tf <- tframe(x) # this generates a default
  tframe(x) <- tf
  x
 }


###############################################

#  Generic .tframe methods (these act on the tframe not on the data)

###############################################


#tfprint.tframe <- function(x, ...) UseMethod("tframePrint")
tfprint.tframe <- function(x, ...) UseMethod("print")
#tframePrint <- function(x, ...) UseMethod("tframePrint")

#tframePrint.default <- function(x, digits=NULL, quote=TRUE, prefix="", ...) 
#  {if (! is.tframe(x)) x <- tframe(x)
#   invisible(print(unclass(x), quote=quote, prefix=prefix, ...)) }

print.tframe <- function(x, ...) invisible(print(unclass(x), ...))


tfTruncate.tframe <- function(x, start=NULL, end=NULL)
    {# like window but uses indexes rather than dates 
     if (!is.null(end))   x[2] <- x[1] + (end-1)/x[3]
     if (!is.null(start)) x[1] <- x[1] + (start-1)/x[3]
     x
    }



tfExpand.tframe <- function(x, add.start=0, add.end=0)
    {x[2] <- x[2] + add.end/x[3]
     x[1] <- x[1] - add.start/x[3]
     x
    }


checktframeConsistent <- function(tf, x) UseMethod("checktframeConsistent")

checktframeConsistent.default <- function(tf, x) Tobs(tf) == Tobs(x)

testEqualtframes <- function(tf1, tf2) UseMethod("testEqualtframes")

testEqualtframes.default <- function(tf1, tf2) { all(tf1==tf2)}



# Following could be used to do date comparisons like tfstart() < tfend()


earliestStart <- function(x, ...)
    tfstart(append(list(x),list(...))[[earliestStartIndex(x, ...)]])

earliestStartIndex <- function(x, ...) UseMethod("earliestStartIndex")

earliestStartIndex.default <- function(x, ...)
  {tf <- list(tframe(x))
   for (i in list(...)) tf <- append(tf, list(tframe(i)))
   do.call("earliestStartIndex", tf) #dispatch on 1st element of tf
  }

earliestStartIndex.tframe <- function(x, ...) 
    {r <- 1
     fr <- tffrequency(x)
     args <- list(x, ...)
     for (i in seq(length(args)))
         {tf <- args[[i]]
          if (tf[3] != fr) stop("frequencies must be that same.")
          if (tf[1] < args[[r]][1]) r <- i
         }           
     r
    }




earliestEnd <- function(x, ...)
    tfend(append(list(x),list(...))[[earliestEndIndex(x, ...)]])

earliestEndIndex <- function(x, ...) UseMethod("earliestEndIndex")

earliestEndIndex.default <- function(x, ...)
  {tf <- list(tframe(x))
   for (i in list(...)) tf <- append(tf, list(tframe(i)))
   do.call("earliestEndIndex", tf) #dispatch on 1st element of tf
  }

earliestEndIndex.tframe <- function(x, ...) 
    {r <- 1
     fr <- tffrequency(x)
     args <- list(x, ...)
     for (i in seq(length(args)))
         {tf <- args[[i]]
          if (tf[3] != fr) stop("frequencies must be that same.")
          if (tf[2] < args[[r]][2]) r <- i
         }           
     r
    }



latestStart <- function(x, ...)
    tfstart(append(list(x),list(...))[[latestStartIndex(x, ...)]])

latestStartIndex <- function(x, ...) UseMethod("latestStartIndex")

latestStartIndex.default <- function(x, ...)
  {tf <- list(tframe(x))
   for (i in list(...)) tf <- append(tf, list(tframe(i)))
   do.call("latestStartIndex", tf)
  }


latestStartIndex.tframe <- function(x, ...) 
    {r <- 1
     fr <- tffrequency(x)
     args <- list(x, ...)
     for (i in seq(length(args)))
         {tf <- args[[i]]
          if (tf[3] != fr) stop("frequencies must be that same.")
          if (tf[1] > args[[r]][1]) r <- i
         }           
     r
    }



latestEnd <- function(x, ...)
    tfend(append(list(x),list(...))[[latestEndIndex(x, ...)]])

latestEndIndex <- function(x, ...) UseMethod("latestEndIndex")

latestEndIndex.default <- function(x, ...)
  {tf <- list(tframe(x))
   for (i in list(...)) tf <- append(tf, list(tframe(i)))
   do.call("latestEndIndex", tf)
  }

latestEndIndex.tframe <- function(x, ...) 
    {r <- 1
     fr <- tffrequency(x)
     args <- list(x, ...)
     for (i in seq(length(args)))
         {tf <- args[[i]]
          if (tf[3] != fr) stop("frequencies must be that same.")
          if (tf[2] > args[[r]][2]) r <- i
         }           
     r
    }




###############################################

#  stamped specific methods   <<<<<<<<<<<<
#  stamped class TS have a date/time stamp associated with each time point
################################################

#checktframeConsistent.stamped <- function(tf, x)
#  {Tobs(x) == Tobs(tf)}

testEqualtframes.stamped <- function(tf1, tf2)
  {all(tf1$stamp == tf2$stamp)}

Tobs.stamped <- function(x) length(tframe(x))

###############################################

testEqual <- function(obj1, obj2, fuzz = 0) UseMethod("testEqual")

testEqual.default <- function(obj1, obj2, fuzz=1e-16) 
  {if      (is.null(obj1)) is.null(obj2)
   else if (is.array(obj1)) testEqual.array(obj1, obj2, fuzz=fuzz)
   else if (is.numeric(obj1)) testEqual.numeric(obj1, obj2, fuzz=fuzz)
   else if (is.list(obj1)) testEqual.list(obj1, obj2, fuzz=fuzz)
   else is.logical(all.equal(obj1, obj2, tolerance=fuzz))
  }

testEqual.array <- function(obj1, obj2, fuzz=1e-16) 
  {if(!is.array(obj2))                     r <-FALSE
   else if (any(dim(obj1) != dim(obj2)))   r <- FALSE
   else if ("character" == mode(obj1))     r <- all(obj1 == obj2)
   else if ("numeric" == mode(obj1))
              r <- testEqual.numeric(obj1, obj2, fuzz=fuzz)
   else stop(paste("matrix of mode ", mode(obj1), " not testable."))
   if (is.na(r))  r <- FALSE
    r
  }

testEqual.matrix <- testEqual.array

testEqual.numeric <- function(obj1, obj2, fuzz=1e-16) 
  {r <- all(is.infinite(obj1) == is.infinite(obj2))
   if (r) 
          {nna <- !is.na(c(obj1))
           r <- fuzz >= max(abs(c(obj1)[nna] - c(obj2)[nna]))
          }
   if (is.na(r))  r <- FALSE
   r
  }

testEqual.list <- function(obj1, obj2, fuzz=1e-16) 
  {r <- length(obj1) == length(obj2)
   if (r) for (i in seq(length(obj1)))
        {if(r) r <- testEqual(obj1[[i]], obj2[[i]], fuzz=fuzz) }
   r
  }

#if (!exists("lag")) lag <- function(x, ...) UseMethod("lag")

#if (!exists("lag.default"))  lag.default <- function(x, ...) {stop("no lag function") }



###############################################

# Time dimension methods for data manipulation

###############################################



splice <- function(mat1, mat2, ...) UseMethod("splice")

splice.default <- function(mat1, mat2, ...)
{#  (... further arguments, currently disregarded)
 # splice together 2 time series matrices. If data  is provided in both for
 #  a given period then mat1 takes priority.
 # The result starts at the earlier of mat1 and mat2 and ends at the later.
 # dimnames are taken from mat1.
 # The frequencies should be the same.
 if (is.null(mat1)) return(mat2)
 if (is.null(mat2)) return(mat1)
 freq <- tffrequency(mat1)
 if (freq != tffrequency(mat2)) stop("frequencies must be the same.")
 p <- nseries(mat1)
 if (p != nseries(mat2))   stop("number of series must be the same.")
 fr <- c(freq,1)
 st <- min(fr %*% tfstart(mat1), fr %*% tfstart(mat2))
 strt <- c(st %/% freq, st %% freq)
 en <- max(fr %*% tfend(mat1), fr%*% tfend(mat2))
 r1 <-r2 <-tframed(matrix(NA, 1+en-st, p), list(start=strt, frequency=freq))
 r1[c((fr %*% tfstart(mat1))-st) + 1:Tobs(mat1),] <- mat1
 r2[c((fr %*% tfstart(mat2))-st) + 1:Tobs(mat2),] <- mat2
 na <- is.na(r1)
 r1[na] <- r2[na] # put mat2 only in na locations of mat1
 #dimnames(r1)<-list(round(time(r1),digits=3),dimnames(mat1)[[2]])
 seriesNames(r1)<- seriesNames(mat1)
 r1 <- tframed(r1, list(start=earliestStart(mat1,mat2), 
                        end =latestEnd(mat1,mat2), frequency=freq))
 r1
}


tfTruncate <- function(x, start=NULL, end=NULL) 
    if(is.null(x)) return(NULL) else UseMethod("tfTruncate")
  # similar to window but start and end specify periods relative to the 
  #   beginning (eg x[start:end] for a vector).
  #   NULL means no truncation.


tfTruncate.default <- function(x, start=NULL, end=NULL)
    {tf <- tfTruncate(tframe(x), start, end)
     if (is.null(start)) start <- 1
     if (is.matrix(x)) 
        {if (is.null(end)) end <- dim(x)[1]
         z <- x[start:end,,drop=FALSE]
        }
     else 
        {if (is.null(end)) end <- length(x)
         z <- x[start:end]
        }
     tframe(z) <- tf
     z
    }

tfExpand <- function(x, add.start=0, add.end=0)  
    if(is.null(x)) return(NULL) else UseMethod("tfExpand")
  # expand (a tframe) by add.start periods on the beginning
  # and add.end Tobs on the end

tfExpand.default <- function(x, add.start=0, add.end=0)
    {tf <- tfExpand(tframe(x), add.start=add.start, add.end=add.end)
     selectSeries(tbind(x, time(tf)), series = -(nseries(x)+1))
    }


trimNA <- function(x, startNAs= TRUE, endNAs= TRUE) UseMethod("trimNA") 

trimNA.default <- function(x, startNAs= TRUE, endNAs= TRUE)
{# trim NAs from the ends of a ts matrix or vector.
 # (Observations for all series are dropped in a given period if any 
 #  one contains an NA in that period.)
 # if startNAs=F then beginning NAs are not trimmed.
 # If endNAs=F   then ending NAs are not trimmed.
 sample <- ! if (is.matrix(x)) apply(is.na(x),1, any) else is.na(x)
 if (!any(sample)) warning("data is empty after triming NAs.")
 s <- if (startNAs) min(time(x)[sample]) else tfstart(x)
 e <- if (endNAs)   max(time(x)[sample]) else tfend(x)
 tfwindow(x, start=s, end=e, warn=FALSE)
}

aggregate.tframed <- function (x, ...)
   {tf <- tframe(x)
    nm <- seriesNames(x)
    # this is assuming tf is actual a ts tframe
    r <- aggregate(ts(unclass(x), start=tf[1], end=tf[2], frequency=tf[3]), ...)
    tframed(r, tf=tframe(r), names=nm)
   }

###############################################

# Non-time dimension methods

###############################################


nseries <- function(x) UseMethod("nseries") 
nseries.default <- function(x)   if (is.null(x)) 0 else NCOL(x)

   

 seriesNames <- function(x)       UseMethod("seriesNames")
"seriesNames<-" <- function(x, value)UseMethod("seriesNames<-")

 seriesNames.default <- function(x) {
    if (is.null(x)) return(NULL)
    else names <- if (is.matrix(x)) dimnames(x)[[2]] else attr(x, "seriesNames")
    if (is.null(names)) names <- paste("Series", seq(NCOL(x)))
    names
    }

"seriesNames<-.default" <- function(x, value) {
   if (is.null(value)) {
      if (is.matrix(x)) dimnames(x)[[2]] <- NULL
      attr(x,"seriesNames") <- NULL
      }      
   else {
      if (mode(value) != "character") value <- seriesNames(value)
      if (length(value) != nseries(x))
         stop("length of names (",length(value),
	      ") does not match number of series(",nseries(x),").")
      if (is.matrix(x)) dimnames(x) <- list(dimnames(x)[[1]], value)
      else attr(x,"seriesNames") <- value
      }
   x
   }



selectSeries <- function(x, series=seqN(nseries(x))) UseMethod("selectSeries")

selectSeries.default <- function(x, series=seqN(nseries(x))) {
  names <- seriesNames(x)
  if (is.character(series)) series <- match(names,series, nomatch=0) > 0
  if(all(0==series) | is.null(series)) r <- NULL
  else if (!is.matrix(x)) r <- x  # vector case
  else {
#    r <- classed(tframed(x[, series, drop = FALSE], tframe(x)), class(x))# reconstructor
#   tframe assignment cannot guarantee that the object has the right structure
#   for a class, so above can give a deformed object in the class.
    r <- tframed(x[, series, drop = FALSE], tframe(x))
    seriesNames(r) <- names[series]
    }
  r
  }


tbind <- function(x, ..., pad.start=TRUE, pad.end=TRUE, warn=TRUE)
   UseMethod("tbind")


tbind.default <- function (x, ..., pad.start = TRUE, pad.end = TRUE, warn = TRUE) 
{
    if (is.null(x)) {
	#stop("first argument cannot be NULL.")
        r <- list(...)
	if(1 != length(r))
	    stop("If first argument is NULL then only one other series should be supplied.")
	return(r[[1]])
	}
    fr <- tffrequency(x)
    for (i in list(...)) {
        if (!is.null(i) && (fr != tffrequency(i))) 
            stop("frequencies must be the same.")
    }
    fr <- c(fr, 1)
    st <- fr %*% tfstart(x)
    for (i in list(...)) if (!is.null(i)) 
        st <- min(st, fr %*% tfstart(i))
    en <- fr %*% tfend(x)
    for (i in list(...)) if (!is.null(i)) 
        en <- max(en, fr %*% tfend(i))
    r <- NULL
    sn <- NULL
    nm <- attr(x, "names")
    attr(x, "names") <- NULL
    for (z in append(list(x), list(...))) {
        if (!is.null(z)) {
            if (is.matrix(z)) {
                if (st == (fr %*% tfstart(z))) 
                  before <- NULL
                else before <- matrix(NA, (fr %*% tfstart(z)) - 
                  st, dim(z)[2])
                if (en == (fr %*% tfend(z))) 
                  aft <- NULL
                else aft <- matrix(NA, en - (fr %*% tfend(z)), 
                  dim(z)[2])
                r <- cbind(r, rbind(before, z, aft))
            }
            else {
                if (st == (fr %*% tfstart(z))) 
                  before <- NULL
                else before <- rep(NA, (fr %*% tfstart(z)) - 
                  st)
                if (en == (fr %*% tfend(z))) 
                  aft <- NULL
                else aft <- rep(NA, en - (fr %*% tfend(z)))
                r <- cbind(r, c(before, z, aft))
            }
            sn <- c(sn, seriesNames(z))
        }
    }
    if (!is.null(nm)) 
        dimnames(r) <- list(nm, NULL)
    if (length(sn) == ncol(r)) 
        seriesNames(r) <- sn
    r <- tframed(r, list(start = c((st - 1)%/%fr[1], 1 + (st - 
        1)%%fr[1]), frequency = fr[1]))
    if (!(pad.start & pad.end)) 
        r <- trimNA(r, startNAs = !pad.start, endNAs = !pad.end)
    if (is.null(r)) 
        warning("intersection is NULL")
    r
}
