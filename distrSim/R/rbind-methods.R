### rbind methods

setMethod("rbind", "ANY",  function(x,...){cal0 <- sys.call(sys.parent())[-1]
          f <- base::rbind; do.call(f,cal0)})

setMethod("rbind", "SeqDataFrames",  function(x,...){
    SL <- list(...)
    len <- length(SL)
    if (len >= 1)
       { if (!all(lapply(SL@data, is.data.frame)))
              stop("all elements must be data frames")
          f <- function(y) {list(ncol(y), names(y))}
          g <- function(y) identical(f(y), f(y@data[[1]]))
          if (!all(lapply(SL@data, g)))
               stop("all elements must have the same column structure")
          vLl <- lapply(SL, function(y) length(y@data))
          vLs <- sum(as.numeric(vLl))+length(x@data)
          vL <- vector("list", vLs)
          S <- 0
          for( i in 1:(len+1))
             {xs <- if(i==1) x else SL[[i-1]]
              LS <- length(xs@data)
              for (j in 1: LS)
                   {vL[[j+S]] <- xs@data[[j]]
                    if(!is.null(names(xs@data)))
                       names(vL)[j+S] <- names(xs@data)[j+S]
                   }
              S <- S + LS
              }
          x <- new("SeqDataFrames", data = vL)
       }else  return(x)
       })
