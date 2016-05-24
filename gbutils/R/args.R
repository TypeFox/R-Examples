## from package pcts
                                                                  # zapsmall  - chops numbers
# 2013-04-02 - obsolete, used it exclusively with from=1, seq_len() is better.
# seqinc <-                                    # returns a length(0) numeric vector if from>to
#   function(from,to){                         #   useful for "for" loops
#     seq(from,length.out=max(0,to-from+1))    # do we need "by" or other arguments?
#   }

shiftright <- function(x,k=1){                          # works only for 0 <= k <= length(x) !
    if(k == 0)             # 2014-02-02 bug fix was returning numeric(0) when k=0 !
        return(x)
   ind <- 1:(length(x)-k)
   c(x[-ind],x[ind])
}

shiftleft <- function(x,k=1){                           # works only for 0 <= k <= length(x) !
    if(k == 0)             # 2014-02-02 bug fix was returning numeric(0) when k=0 !
        return(x)

   ind <- seq_len(k)   # basically ind <- 1:k, but empty vector if k==0 to work for k=0.
   c(x[-ind],x[ind])   # TODO: however, this line will not work properly for k == 0,
                       #       therefore the "if" remedy above.
}

myouter <-
function(x,y,fun){
   res <- matrix(0, nrow=length(x), ncol=length(y))
   for(i in 1:length(x))
      for(j in 1:length(y))
         res[i,j] <- fun(x[i],y[j])
   res
}

nposargs <- function(x,a=FALSE){# x is expected to be a call, usually the result of sys.call()
  wrk <- as.list(x)
  tot <- length(wrk) - 1
  nameswrk <- names(wrk)
  named <- if(!is.null(nameswrk)) length(nameswrk[nameswrk!=""])
           else                   0
  res <- tot - named

  # a patch needed for my purposes follows
  if( named == 0 && res == 2 && a[1] )
    res <- res-1
  res
}

## 2015-02-12 removed myargnames() since it  has identical functionality to allNames
##
## myargnames <- function(x){                                  # x is expected to be a list
##   nameswrk <- names(x)
##   res  <- if(is.null(nameswrk))
##             character(length(x))
##           else
##             nameswrk
##   res
## }

isargunnamed <- function(x,k){
  if( 0<k && k<=length(x))
    # 2015-02-12 was: identical(myargnames(x)[k],"")
    identical(allNames(x)[k],"")
  else
    FALSE
}
