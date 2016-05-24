setMethod("[", "SeqDataFrames", function(x, i, j, k, drop = FALSE){
          if(missing(k)) k <- seq_len(length(x@data))
          kl <- length(k)
          if (kl == 1){
              daf <- "["(x@data[[k]],i,j, drop = drop)
              if (drop)
                   return(daf)
              else return(new("SeqDataFrames", data = list(daf)))
          }else {
              kn <- seq_len(length(x@data))
              if(!is.null(names(x@data)))
                 names(kn) <- names(x@data)
              kl0 <- kn[k]
              kll <- length(kl0)
              lis <- vector("list",kll)
              for (kk in seq_len(kll))
                  {lis[[kk]] <- as.data.frame("["(x@data[[kl0[kk]]], i,j, drop = drop))
                   if(!is.null(names(x@data)))
                      names(lis)[kk] <- names(x@data)[kl0[kk]]}
              return(new("SeqDataFrames", data = lis))
          }})

setReplaceMethod("[", "SeqDataFrames", function(x, i, j, k, value){
          if(missing(k)) k <- seq_len(length(x@data))
          if(length(k)==1){
             if((k<=length(x@data))||!is(try(x@data[[k]],silent=TRUE),"try-error"))
                {zl <- x@data
                 z  <- zl[[k]]
                 if (missing(i))
                    { i <- seq_len(nrow(z))
                     #if(!is.null(dim(value)))
                     #    z <- data.frame(matrix(NA,nrow(value),ncol(x@data[[1]])))
                     #    else z <- data.frame(matrix(NA,length(value),ncol(x@data[[1]])))
                    }
                 z[i,j] <- value
                 zl[[k]] <- z
                 x@data <- zl
               }else{
                 if(!is.null(dim(value)))
                      z <- data.frame(matrix(NA,nrow(value),ncol(x@data[[1]])))
                 else z <- data.frame(matrix(NA,length(value),ncol(x@data[[1]])))
                 z[i,j] <- value
                 x@data <- c(x@data,list(z))
              }
             return(x)}

          if(missing(j)) j <- seq_len(ncol(x@data[[1]]))
          if(missing(i)) i <- lapply(seq_len(length(x@data)),function(y) seq_len(nrow(x@data[[y]])))

          if(is(value, "SeqDataFrames")) value <- value@data

          kn <- seq_len(length(x@data))
          if(!is.null(names(x@data)))
                 names(kn) <- names(x@data)
          kl0 <- kn[k]
          kll <- length(kl0)

## Is the following line correct?
## should it be: if(!is.list(i)) i <- lapply(kl0, function(y) y)
## or: if(!is.list(i)) i <- as.list(kl0)?


          if(!is.list(i)) i <- lapply(kl0, function(y) i)

          if(is(value,"atomic"))
             value <- lapply(seq_len(kll),
                             function(y) data.frame(matrix(
                                         rep(value, length( i[[kl0[y]]])*length(j)),
                                         length(i[[kl0[y]]]),
                                         length(j) )
                                         )
                             )
          if (kll==1) value <- list(value)
          if(is(value,"data.frame"))
                value <- lapply(kl0,function(y) value)

          zl <- x@data
          for(kk in seq_len(kll))
                 {z <- zl[[kl0[kk]]]
                  z[c(unlist(i[[kl0[kk]]])),j] <- value[kk]
                  zl[[kl0[kk]]] <- z
                  }
          x@data <-  zl
          return(x)
          })
