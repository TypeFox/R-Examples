#source("d:/MWP/eAnalysis/rindex/R/rindex.R")

rindexAutobatch <- function(n, batch=64, basis=2, minimum=3){
  as.integer(max(minimum,ceiling(n/(basis^max(ceiling(log(n,basis)-log(batch,basis)), 0)))))
}

rindex <- function(x, uni=NULL, batch=NULL, verbose=FALSE){
  stopifnot(mode(x)=="character")
  n <- length(x)
  nNA <- sum(is.na(x))
  sorttime <- system.time({
    o <- order(x, na.last=TRUE)
    x <- x[o]
    # xx instead of duplicated tune by special C function scanning the already sorted x
    if (is.null(uni))
      uni <- !any(duplicated(x))
  })
  if (is.null(batch))
    batch <- rindexAutobatch(n-nNA)
  obj <- list(
    val   = x
  , pos   = o
  , n     = length(o)
  , nNA   = sum(is.na(x))
  , batch = batch
  , uni   = uni
  , tree = new.env(parent=globalenv())
  )
  class(obj) <- "rindex"
  treetime <- system.time(
    obj <- rindexAddTree(obj)
  )
  if (verbose)
    print(rbind(sort=sorttime, tree=treetime)[,1:3])
  obj
}

c.rindex <- function(...){
  l <- list(...)
  k <- length(l)
  if (k>1){
    l[[1]] <- rindexDelTree(l[[1]])
    pos <- l[[1]]$pos
    val <- l[[1]]$val
    n <- l[[1]]$n
    nNA <- l[[1]]$nNA
    for (i in 2:k){
      l[[i]] <- rindexDelTree(l[[i]])
      l[[i]]$pos <- l[[i]]$pos + n
      n <- n + l[[i]]$n
      nNA <- nNA + l[[i]]$nNA
    }
    pos <- do.call("c", lapply(l,function(i)i$pos))
    val <- do.call("c", lapply(l,function(i)i$val))
    o <- order(val)         # pairewise xx merge sort would be faster
    val <- val[o]
    pos <- pos[o]
    uni <- !any(duplicated(val)) # xx scanning along sorted values would be faster than duplicated()
    obj <- list(
      val   = val
    , pos   = pos
    , n     = n
    , nNA   = nNA
    , batch = rindexAutobatch(n)
    , uni   = uni
    , tree  = new.env(parent=globalenv())
    )
    class(obj) <- "rindex"
    obj <- rindexAddTree(obj)
    obj
  }else if (k==1){
    l[[1]]
  }else{
    NULL
  }
}

rindexAddTree <- function(obj, batch=NULL)
{
  if (is.null(batch)){
    batch <- obj$batch
  }else{
    if (batch<as.integer(3))
      stop("minum batch size is 3")
    obj$batch <- as.integer(batch)
  }
  if (obj$n-obj$nNA){
    val <- obj$val
    if (obj$uni){
      recuni <- function(low,high){
        if ((high-low+1)>batch){
          mid <- as.integer((high+low)/2)
          return(list(mid=val[mid], low=Recall(low, mid), high=Recall(mid+as.integer(1), high)))
        }else{
          return(list(mid=NA, low=low, high=high))
        }
      }
      tree <- recuni(as.integer(1), obj$n-obj$nNA)
    }else{
      rectie <- function(low,high,flexbatch){
        if ((high-low+1)>flexbatch){
          mid <- as.integer((high+low)/2)
          # make sure that all equal values end up in the same leave
          # thus, may need to adjust mid for ties and maybe skip splitting (if everything equal)
          midval <- val[mid]
          goon <- TRUE
          midlow <- mid
          midhigh <- mid
          takewhich <- 0
          while(goon){
            if (midhigh<high){
              midhigh <- midhigh + as.integer(1)
              if (val[midhigh]!=midval){
                takewhich <- as.integer(+1)
                break
              }
            }else{
              goon <- FALSE
            }
            if (midlow>low){
              midlow <- midlow - as.integer(1)
              if (val[midlow]!=midval){
                takewhich <- as.integer(-1)
                break
              }
            }else{
              goon <- FALSE
            }
          }
          if (takewhich>0){
            if (midhigh==mid+as.integer(1))  # the default if no ties are present
              return(list(mid=midval
              , low =Recall(low, midhigh-as.integer(1), flexbatch=batch)
              , high=Recall(midhigh             , high, flexbatch=batch)
              ))
            else
              return(list(mid=midval
              , low =Recall(low, midhigh-as.integer(1), flexbatch=rindexAutobatch(midhigh-low   , batch=obj$batch))  # taking obj$batch here assures that batch does not recursively shrink
              , high=Recall(midhigh             , high, flexbatch=rindexAutobatch(high-midhigh+1, batch=obj$batch))  # taking obj$batch here assures that batch does not recursively shrink
              ))
          }else if (takewhich<0){
            return(list(mid=val[midlow]
            , low =Recall(low                 , midlow, flexbatch=rindexAutobatch(midlow-low+1, batch=obj$batch))    # taking obj$batch here assures that batch does not recursively shrink
            , high=Recall(midlow+as.integer(1), high  , flexbatch=rindexAutobatch(high-midlow , batch=obj$batch))    # taking obj$batch here assures that batch does not recursively shrink
            ))
          }
        }
        # if (low-high+1)<=flexbatch OR all values equal (takewhich==0)
        return(list(mid=NA, low=low, high=high))
      }
      tree <- rectie(as.integer(1), obj$n-obj$nNA, flexbatch=batch)
    }
    assign("tree", tree, envir=obj$tree)
  }
  obj
}

rindexDelTree <- function(obj){
  if (obj$n-obj$nNA){
    if (exists("tree", envir=obj$tree))
      remove("tree", envir=obj$tree)
  }
  obj
}

rindexNodes <- function(obj)
{
  if (obj$n-obj$nNA){
    rec <- function(node){
      if (is.na(node$mid)){
        return(as.integer(1))
      }else{
        return(as.integer(1)+Recall(node$high)+Recall(node$low))
      }
    }
    if (exists("tree", obj$tree))
      return(rec(get("tree", envir=obj$tree)))
  }
  return(as.integer(0))
}

rindexBytes <- function(obj)
{
  if (obj$n-obj$nNA){
    if (exists("tree", obj$tree))
      return(c(val=object.size(obj$val), pos=object.size(obj$pos), tree=object.size(get("tree", envir=obj$tree))))
  }
  return(c(val=object.size(obj$val), pos=object.size(obj$pos), tree=as.integer(0)))
}

print.rindex <- function(x, tree=FALSE, ...){
  nNodes <- rindexNodes(x)
  nBytes <- rindexBytes(x)
  tBytes <- sum(nBytes)
  if (nNodes && tree){
    rec <- function(level, node){
      if (level)
        cat(paste(rep(" ", level), collapse=" ", sep=""))
      if (is.na(node$mid)){
        cat("low=", node$low, "  mid=", NA, "  high=", node$high, "\n", sep="")
      }else{
        cat("low=", NA, "  mid=", node$mid, "  high=", NA, "\n", sep="")
        Recall(level+1,node$low)
        Recall(level+1,node$high)
      }
    }
    if (exists("tree", x$tree))
      rec(0, get("tree", envir=x$tree))
  }
  cmod <- mode(x$val)
  cat(if (x$uni)"unique", cmod, "rindex of length", x$n, "with",x$nNA,"NAs at batch size", x$batch,"\n")
  print(data.frame(Bytes=nBytes, Percent=paste(format(nBytes/tBytes, digits=2), "%", sep="")))
}


length.rindex <- function(x){
  x$n
}

names.index <- function(x){
  stop("names currently not supported on rindex")
}

"names<-.rindex" <- function(x, value){
  stop("assignment of names currently not supported on rindex")
}

str.rindex <- function(object, ...){
  cat("modified str for object of class `rindex'\n")
  object <- unclass(object)
  NextMethod("str")
}


sort.rindex <- function(x, decreasing = FALSE, na.last = NA, ...){
  if (length(list(...)))
    stop("sort.rindex only allows arguments x=rindex, decreasing and na.last")
  if (!inherits(x, "rindex"))
    stop("first argument to sort.rindex must inherit from class rindex")
  if (x$n-x$nNA){
    if (decreasing){
      ret <- x$val[(x$n-x$nNA):as.integer(1)]
    }else{
      ret <- x$val[as.integer(1):(x$n-x$nNA)]
    }
  }else{
    ret <- x$val[0]
  }
  if (is.na(na.last) || !x$nNA){
    return(ret)
  }else if (na.last){
    return(c(ret, rep(NA, x$nNA)))
  }else{
    return(c(rep(NA, x$nNA), ret))
  }
}

#note that order(decreasing=TRUE) does not inverst the positions of NA values !!
order.rindex <- function(..., na.last = TRUE, decreasing = FALSE){
  if (length(list(...))!=1)
    stop("order.rindex requires exactly one dotted argument")
  x <- list(...)[[1]]
  if (!inherits(x, "rindex"))
    stop("first argument to order.rindex must inherit from class rindex")
  if (x$n-x$nNA){
    if (decreasing){
      ret <- x$pos[(x$n-x$nNA):as.integer(1)]
      if (!x$uni)
        warning("order.rindex(..., decreasing=FALSE) handles ties diferent from order(..., decreasing=FALSE)")
    }else{
      ret <- x$pos[as.integer(1):(x$n-x$nNA)]
    }
  }else{
    ret <- integer()
  }
  if (is.na(na.last) || !x$nNA){
    return(ret)
  }else if (na.last){
    return(c(ret, x$pos[(x$n-x$nNA+as.integer(1)):x$n]))
  }else{
    return(c(x$pos[(x$n-x$nNA+as.integer(1)):x$n], ret))
  }
}


"[.rindex" <- function(x, i, ...){
  z <- x$val
  z[x$pos] <- x$val
  x <- z
  NextMethod("[")
}

"[<-.rindex" <- function(x, i, ..., value){
  stop("assignment to rindex not supported by static readonly ramtree")
}


is.na.rindex <- function(x){
  ret <- logical(x$n)
  if (x$nNA)
    ret[x$pos[(x$n-x$nNA+as.integer(1)):x$n]] <- TRUE
  ret
}


rindexFind <- function(obj, val, findlow=TRUE)
{
  if (obj$n){
    if (is.na(val)){
      if (obj$nNA){
        if (findlow){
          return(as.integer(c(0,obj$n-obj$nNA+as.integer(1))))
        }else{
          return(c(as.integer(0),obj$n))
        }
      }
    }else{
      v <- obj$val
      rec <- function(node){
        if (is.na(node$mid)){
          if (findlow){
            for (i in node$low:node$high){
              if (v[i]>=val){
                if (v[i]==val){
                  return(c(as.integer(0),i))
                }else{
                  return(c(as.integer(1),i))
                }
              }
            }
          }else{
            for (i in node$high:node$low){
              if (v[i]<=val){
                if (v[i]==val){
                  return(c(as.integer(0),i))
                }else{
                  return(c(as.integer(1),i))
                }
              }
            }
          }
          return(as.integer(c(-1,-1)))
        }else{
          if (val>node$mid){
            Recall(node$high)
          }else{
            Recall(node$low)
          }
        }
      }
      if (!rindexNodes(obj))
        rindexAddTree(obj)
      return(rec(get("tree", envir=obj$tree)))
    }
  }
  return(as.integer(c(-1,-1)))
}


rindexFindlike <- function(obj, val, findlow=TRUE)
{

  if (obj$n){
    if (is.na(val)){
      if (obj$nNA){
        if (findlow){
          return(as.integer(c(0,obj$n-obj$nNA+as.integer(1))))
        }else{
          return(c(as.integer(0),obj$n))
        }
      }
    }else{
      v <- obj$val
      l <- nchar(val)
      islikehigh <- function(node){
        if (is.na(node$mid)){
          val==substr(obj$val[node$low],1,l)
        }else{
          Recall(node$low)
        }
      }
      rec <- function(node){
        if (is.na(node$mid)){
          if (findlow){
            for (i in node$low:node$high){
              vi <- substr(v[i],1,l)
              if (vi>=val){
                if (vi==val){
                  return(c(as.integer(0),i))
                }else{
                  return(c(as.integer(1),i))
                }
              }
            }
          }else{
            for (i in node$high:node$low){
              vi <- substr(v[i],1,l)
              if (vi<=val){
                if (vi==val){
                  return(c(as.integer(0),i))
                }else{
                  return(c(as.integer(1),i))
                }
              }
            }
          }
          return(as.integer(c(-1,-1)))
        }else{
          if (val>substr(node$mid,1,l) || (!findlow && islikehigh(node$high)) ){
            Recall(node$high)
          }else{
            Recall(node$low)
          }
        }
      }
      if (!rindexNodes(obj))
        rindexAddTree(obj)
      return(rec(get("tree", envir=obj$tree)))
    }
  }
  return(as.integer(c(-1,-1)))
}


  rindexFindInterval <- function(obj
, low=NULL, high=NULL
, low.include=TRUE, high.include=TRUE
, low.exact=FALSE, high.exact=FALSE
, lowFUN=FUN, highFUN=FUN
, FUN=rindexFind
)
{
  n <- obj$n
  if (!n)
    return(integer())

  countNA <- 0

  if (is.null(low)){
    ilow <- as.integer(1)
  }else{
    if (is.na(low))
      countNA <- countNA + 1
    flow <- lowFUN(obj, low, findlow=TRUE)
    if (flow[1]==0){
      if (low.include){
        ilow <- flow[2]
      }else {
        ilow <- lowFUN(obj, low, findlow=FALSE)[2] + as.integer(1)
        if (ilow>(n-obj$nNA))
          return(integer())
      }
    }else if(flow[1]==1 && !low.exact){
        ilow <- flow[2]
    }else{
        return(integer())
    }
  }

  if (is.null(high)){
    ihigh <- n - obj$nNA
  }else{
    if (is.na(high))
      countNA <- countNA + 1
    if (obj$uni && !is.null(low) && low==high && flow[1]==0 && low.include && high.include){
      ihigh <- ilow
    }else{

      fhigh <- highFUN(obj, high, findlow=FALSE)
      if (fhigh[1]==0){
        if (high.include){
          ihigh <- fhigh[2]
        }else {
          ihigh <- highFUN(obj, high, findlow=TRUE)[2] - as.integer(1)
          if (ihigh<as.integer(1))
            return(integer())
        }
      }else if(fhigh[1]==1 && !high.exact){
          ihigh <- fhigh[2]
      }else{
          return(integer())
      }

      if (ihigh<ilow)
        return(integer())

    }
  }

  if ( countNA == 1 )
    stop("NAs and normal values must not be mixed")

  return(ilow:ihigh)
}


rindexMatch <- function(obj, x, findlow=TRUE, what=c("ind", "pos", "val")){
  what <- match.arg(what)
  if (!obj$uni && missing(findlow))
    warning("indexMatch used with non-unique index, returning ", if (findlow) "first" else "last", " match")
  ret <- as.vector(sapply(x, function(i){
    res <- rindexFind(obj, i, findlow=findlow)
    if (res[1]==0)
      res[2]
    else
      as.integer(NA)
  }))
  if (what=="pos")
    obj$pos[ret]
  else if (what=="val")
    obj$val[ret]
  else
    ret
}


rindexEQ <- function(obj, x, what=c("pos", "val", "ind"), ...){
  what <- match.arg(what)
  if (length(x)==1){
    ret <- rindexFindInterval(obj, x, x, ...)
  }else if (length(x)>1){
    return(lapply(x, function(i){
      ret <- rindexFindInterval(obj, i, i, ...)
      if (what=="pos")
        obj$pos[ret]
      else if (what=="val")
        obj$val[ret]
      else
        ret
    }))
  }else{
    ret <- integer()
  }
  if (what=="pos")
    obj$pos[ret]
  else if (what=="val")
    obj$val[ret]
  else
    ret
}

"==.rindex" <- function(e1,e2){
  if (length(e2)==1){
    ret <- logical(e1$n)
    if (is.na(e2)){
      if (e1$nNA)
        ret[e1$pos[(e1$n-e1$nNA+as.integer(1)):e1$n]] <- TRUE
    }else{
      ret[e1$pos[rindexFindInterval(e1, e2, e2)]] <- TRUE
      if (e1$nNA)
        ret[e1$pos[(e1$n-e1$nNA+as.integer(1)):e1$n]] <- NA
    }
    ret
  }else{
    stop("second parameter of rindex operator must have length 1")
  }
}

rindexNE <- function(obj, x, what=c("pos", "val", "ind"), ...){
  what <- match.arg(what)
  nlow <- as.integer(1)
  nhigh <- obj$n-obj$nNA
  if (length(x)==1){
    ret <- (nlow:nhigh)[is.na(match(nlow:nhigh, rindexFindInterval(obj, x, x, ...)))]
  }else if (length(x)>1){
    return(lapply(x, function(i){
      ret <- (nlow:nhigh)[is.na(match(nlow:nhigh, rindexFindInterval(obj, i, i, ...)))]
      if (what=="pos")
        obj$pos[ret]
      else if (what=="val")
        obj$val[ret]
      else
        ret
    }))
  }else{
    ret <- nlow:nhigh
  }
  if (what=="pos")
    obj$pos[ret]
  else if (what=="val")
    obj$val[ret]
  else
    ret
}

"!=.rindex" <- function(e1,e2){
  if (length(e2)==1){
    ret <- logical(e1$n)
    if (is.na(e2)){
      if (e1$nNA)
        ret[-e1$pos[(e1$n-e1$nNA+as.integer(1)):e1$n]] <- TRUE
    }else{
      ret[-e1$pos[rindexFindInterval(e1, e2, e2)]] <- TRUE
      if (e1$nNA)
        ret[e1$pos[(e1$n-e1$nNA+as.integer(1)):e1$n]] <- NA
    }
    ret
  }else{
    stop("second parameter of rindex operator must have length 1")
  }
}



rindexLT <- function(obj, x, what=c("pos", "val", "ind"), ...){
  what <- match.arg(what)
  if (length(x)==1){
    ret <- rindexFindInterval(obj, NULL, x, high.include=FALSE, ...)
  }else if (length(x)>1){
    return(lapply(x, function(i){
      ret <- rindexFindInterval(obj, NULL, i, high.include=FALSE, ...)
      if (what=="pos")
        obj$pos[ret]
      else if (what=="val")
        obj$val[ret]
      else
        ret
    }))
  }else{
    ret <- integer()
  }
  if (what=="pos")
    obj$pos[ret]
  else if (what=="val")
    obj$val[ret]
  else
    ret
}

"<.rindex" <- function(e1,e2){
  if (length(e2)==1){
    ret <- logical(e1$n)
    if (is.na(e2)){
      stop("rindex<NA not allowed")
    }else{
      if (e1$nNA)
        ret[e1$pos[(e1$n-e1$nNA+as.integer(1)):e1$n]] <- NA
      ret[e1$pos[rindexFindInterval(e1, NULL, e2, high.include=FALSE)]] <- TRUE
    }
    ret
  }else{
    stop("second parameter of rindex operator must have length 1")
  }
}


rindexGT <- function(obj, x, what=c("pos", "val", "ind"), ...){
  what <- match.arg(what)
  if (length(x)==1){
    ret <- rindexFindInterval(obj, x, NULL, low.include=FALSE, ...)
  }else if (length(x)>1){
    return(lapply(x, function(i){
      ret <- rindexFindInterval(obj, i, NULL, low.include=FALSE, ...)
      if (what=="pos")
        obj$pos[ret]
      else if (what=="val")
        obj$val[ret]
      else
        ret
    }))
  }else{
    ret <- integer()
  }
  if (what=="pos")
    obj$pos[ret]
  else if (what=="val")
    obj$val[ret]
  else
    ret
}

">.rindex" <- function(e1,e2){
  if (length(e2)==1){
    ret <- logical(e1$n)
    if (is.na(e2)){
      stop("rindex>NA not allowed")
    }else{
      if (e1$nNA)
        ret[e1$pos[(e1$n-e1$nNA+as.integer(1)):e1$n]] <- NA
      ret[e1$pos[rindexFindInterval(e1, e2, NULL, low.include=FALSE)]] <- TRUE
    }
    ret
  }else{
    stop("second parameter of rindex operator must have length 1")
  }
}



rindexLE <- function(obj, x, what=c("pos", "val", "ind"), ...){
  what <- match.arg(what)
  if (length(x)==1){
    ret <- rindexFindInterval(obj, NULL, x, high.include=TRUE, ...)
  }else if (length(x)>1){
    return(lapply(x, function(i){
      ret <- rindexFindInterval(obj, NULL, i, high.include=TRUE, ...)
      if (what=="pos")
        obj$pos[ret]
      else if (what=="val")
        obj$val[ret]
      else
        ret
    }))
  }else{
    ret <- integer()
  }
  if (what=="pos")
    obj$pos[ret]
  else if (what=="val")
    obj$val[ret]
  else
    ret
}

"<=.rindex" <- function(e1,e2){
  if (length(e2)==1){
    ret <- logical(e1$n)
    if (is.na(e2)){
      stop("rindex<=NA not allowed")
    }else{
      if (e1$nNA)
        ret[e1$pos[(e1$n-e1$nNA+as.integer(1)):e1$n]] <- NA
      ret[e1$pos[rindexFindInterval(e1, NULL, e2, high.include=TRUE)]] <- TRUE
    }
    ret
  }else{
    stop("second parameter of rindex operator must have length 1")
  }
}


rindexGE <- function(obj, x, what=c("pos", "val", "ind"), ...){
  what <- match.arg(what)
  if (length(x)==1){
    ret <- rindexFindInterval(obj, x, NULL, low.include=TRUE, ...)
  }else if (length(x)>1){
    return(lapply(x, function(i){
      ret <- rindexFindInterval(obj, i, NULL, low.include=TRUE, ...)
      if (what=="pos")
        obj$pos[ret]
      else if (what=="val")
        obj$val[ret]
      else
        ret
    }))
  }else{
    ret <- integer()
  }
  if (what=="pos")
    obj$pos[ret]
  else if (what=="val")
    obj$val[ret]
  else
    ret
}

">=.rindex" <- function(e1,e2){
  if (length(e2)==1){
    ret <- logical(e1$n)
    if (is.na(e2)){
      stop("rindex>=NA not allowed")
    }else{
      if (e1$nNA)
        ret[e1$pos[(e1$n-e1$nNA+as.integer(1)):e1$n]] <- NA
      ret[e1$pos[rindexFindInterval(e1, e2, NULL, low.include=TRUE)]] <- TRUE
    }
    ret
  }else{
    stop("second parameter of rindex operator must have length 1")
  }
}


