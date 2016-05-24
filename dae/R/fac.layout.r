"factor.list" <- function(generate, order="standard")
#takes a generate list and creates a list of factor names, with levels 
#information, and a list of factor relative replications, both of which are 
#returned as a list of the two parallel lists called factors and reps.
{ n <- length(generate)
  which.ord <- pmatch(casefold(order), c("standard", "yates"), nomatch="")
	if (which.ord == "")	stop("order must be either standard or yates")
# standard order
	if (which.ord == "1")
    counter <- 1:n
  else
# Yates order
    counter <- n:1
  kfac <- 0
  for(i in counter) 
  { if(!(names(generate[i]) == ""))
    { kfac=kfac+1
      if (kfac == 1)
      { fnames <- list(generate[[i]])
        names(fnames) <- names(generate)[i]
        freps <- 1
      }
      else
      { knames <- list(generate[[i]])
        names(knames) <- names(generate)[i]
        fnames <- c(fnames, knames)
        freps <- c(freps, 1)
      }
    }
    else
    { if (kfac == 0)
        if (which.ord == "1")
          stop("Must start with a factor name - set times argument instead")
        else
          stop("Must end with a factor name - set each argument instead")
      freps[kfac] <- generate[[i]]
    }
  }
  if (which.ord == "2") #reverse lists for Yates order
  { fnames <- fnames[kfac:1]
    freps <- freps[kfac:1]
  }
  return(list(factors = fnames,reps = freps))
}

"fac.rand.cross" <- 
function(unrandomized, unr.names, unr.levels, nested.factors=NULL, except=NULL, 
         randomized, seed=NULL)
{ if (is.data.frame(randomized))
    n <- nrow(randomized)
  else
    n <- length(randomized)
  which.unr <- 0
  if(is.data.frame(unrandomized)) #for data.frame
  { which.unr <- 1
    facgen <- unrandomized
  }
  else 
  { facgen <- fac.gen(generate=unr.names)
  }
#generate random factor values
  nunr <- ncol(facgen)
  for(i in 1:nunr)
  { if (!(names(unr.names)[i] %in% except))
      rno <- runif(unr.levels[i])[as.integer(facgen[[i]])]
    else
      rno <- as.integer(facgen[[i]])
    if(i == 1)
      facrand <- data.frame(rno)
    else
      facrand <- data.frame(facrand,rno)
  }
  names(facrand) <- names(unr.names)
  facrand.ord <- 1:n
  facrand.ord[do.call(order, facrand)] <- facrand.ord
  facrand <- fac.divide(facrand.ord, unr.names)
  facrand
}

"fac.rand.nest" <- 
function(unrandomized, unr.names, unr.levels, nested.factors=NULL, except=NULL, 
         randomized, seed=NULL)
{ nnested <- length(nested.factors)
  names.nested <- names(nested.factors)
#check nested factors for transitivity
  for (i in 1:nnested)
  { nnest <- length(nested.factors[[i]])
    for (j in 1:nnest)
    { knested <- match(nested.factors[[i]][j], names.nested)
      if(!is.na(knested))
      { if(!all(nested.factors[[knested]] %in% nested.factors[[i]]))
        {  stop(names.nested[i]," is nested in ",nested.factors[[i]][j],
               " but is not nested in all those that ",nested.factors[[i]][j]," is.")
        }
      }
    }
  }
#order the unrandomized factors so that any subsets of the generalized factor
#formed from the nesting factors of a nested factor are to the left of the 
#nested factor. Also sort numbers of levels in levels and form vector 
#unr.denestord giving location of unrandomized factors in nestord.
#first a list of all non-nested factors
  if (is.data.frame(randomized))
    n <- nrow(randomized)
  else
    n <- length(randomized)
  which.unr <- 0
  if(is.data.frame(unrandomized)) #for data.frame
  { which.unr <- 1
    nunr <- ncol(unrandomized)
  }
  else 
    nunr <- length(unr.names)
  kunr <- 0
  unr.levels.nestord <- rep(1, length=nunr)
  unr.denestord <- rep(NA, length=nunr)
  for(i in 1:nunr)
  { if(is.na(match(names(unr.names)[i], names.nested)))
    { kunr <- kunr+1
      unr.denestord[i] <- kunr
      unr.levels.nestord[kunr] <- unr.levels[i]
      if(kunr == 1)
      { unr.nestord <- list(unr.names[[i]])
        names(unr.nestord) <- names(unr.names)[i]
      }
      else
      { knames <- list(unr.names[[i]])
        names(knames) <- names(unr.names)[i]
        unr.nestord <- c(unr.nestord, knames)
      }
    }
   } 
#now sort nested factors for number of nesting factors and add to list in this order
  if(!is.null(nested.factors))
  { nested.sort <- sort(sapply(nested.factors, FUN=length))
    for (i in 1:nnested)
    { kunr <- kunr+1
      kunr.name <- names(nested.sort)[i]
      kno <- match(kunr.name, names(unr.names))
      unr.denestord[kno] <- kunr
      unr.levels.nestord[kunr] <- unr.levels[kno]
      if(kunr == 1)
      { unr.nestord <- list(unr.names[[kno]])
        names(unr.nestord) <- names(unr.names)[kno]
      }
      else
      { knames <- list(unr.names[[kno]])
        names(knames) <- names(unr.names)[kno]
        unr.nestord <- c(unr.nestord, knames)
      }
    }
  }
#generate random factor values with the factors in nested order
  facgen <- fac.gen(generate=unr.nestord)
  for(i in 1:nunr)
  { if (names(unr.nestord)[i] %in% except)
      rno <- as.integer(facgen[[i]])
    else
    { knest <- match(names(unr.nestord)[i], names.nested)
      if(is.na(knest))   #nonnested factor
      { rno <- runif(unr.levels.nestord[i])[as.integer(facgen[[i]])]
      }
      else     #nested factor
      { kfac <- length(nested.factors[[knest]])+1
        kfacnos <- rep(1, length=kfac)
        kfacnos[1] <- i
        for (j in 1:(kfac-1))
        { kfacnos[j+1] <- match(nested.factors[[knest]][j], names(unr.nestord))
          if(is.na(kfacnos[j+1]))
            stop("Nesting factor not in list of unrandomized factors.")
        }
        sort(kfacnos)
        #determine number of random nos required and 
        #generate radix to expand to n-vector
        radix <- rep(1, length=n)
        each <- 1
        for (j in kfac:1)
        { radix <- radix + (as.integer(facgen[[kfacnos[j]]])-1)*each
          each <- each*unr.levels.nestord[kfacnos[j]]
        }
        rno <- runif(each)[radix]
      }
    }
    if(i == 1)
      facrand <- data.frame(rno)
    else
      facrand <- data.frame(facrand,rno)
  }
  names(facrand) <- names(unr.nestord)
  facrand.ord <- 1:n
  facrand.ord[do.call(order, facrand)] <- facrand.ord
  facrand <- fac.divide(facrand.ord, unr.nestord)
  #reorder facrand so that the factors and their values are in an appropriate 
  #order for the original factor order i.e. in unr.names
  facgen <- fac.gen(generate=unr.names)
  #form radix from facgen, numbers in radix labelling the levels of the factors
  #in order as for unr.nestord
  radix <- rep(1, length=n)
  each <- 1
  #loop over factors in unr.nestord
  for (j in nunr:1)
  { kno <- match(names(unr.nestord)[j], names(unr.names))
    radix <- radix + (as.integer(facgen[[kno]])-1)*each
    each <- each*unr.levels[kno]
  }
  facrand <- facrand[radix,]
  facrand <- facrand[, unr.denestord]
  facrand
}

"fac.layout" <- 
function(unrandomized, nested.factors=NULL, except=NULL, randomized, seed=NULL, unit.permutation = TRUE)
{
#generate a layout for a design consisting of randomized factors that are 
#randomized to the unrandomized factors, taking into account the nesting between
#the unrandomized factors.
#unrandomized is a data.frame or list of factors (no numbers), along with their 
#levels. If it is a list, each component of the list is a factor name that 
#contains either a single numeric value that is the number of levels, a numeric 
#vector that contains the levels of the factor or a character vector that 
#contains the labels of the factor.
#nested.factor is a list of the factors in unrandomized that are nested in other
#factors in unrandomized. The factors within which they are nested are the 
#elements of the component.
#randomized is a factor or a data frame containing the generated factor or 
#factors to be randomized.
#process randomized argument
  if(!is.data.frame(randomized) & !is.factor(randomized))
    stop("randomized must be a factor or data frame.")
  if (is.data.frame(randomized))
    n <- nrow(randomized)
  else
    n <- length(randomized)
#process seed argument
  if (!is.null(seed))
    set.seed(seed)
#process unrandomized argument and form unrandomized factor list 
  if(is.data.frame(unrandomized)) #for data.frame
  { if(!all(sapply(unrandomized, FUN=is.factor)))
    stop("All columns in the unrandomized data.frame must be factors")
    nunr <- ncol(unrandomized)
    if (nrow(unrandomized) != n)
      stop("The number of rows in the unrandomized data frame and the length of the factor(s) must be equal.")
    unr.names <- vector("list", length = nunr)
    names(unr.names) <- as.list(names(unrandomized))
    for (i in 1:nunr)
      unr.names[[i]] <- levels(unrandomized[[i]])
  }
  else
  { if (!is.list(unrandomized))  #for (generate) list
       stop("unrandomized must be a list or a data.frame.")
    if(any(names(unrandomized) == ""))
      stop("all components of unrandomized list must be named.")
    facs.reps <- factor.list(unrandomized, order="standard")
    unr.names <- facs.reps$factors
    nunr <- length(unr.names)
  }
#if except is not NULL, check that it contains only unrandomized factors
  if (!is.null(except))
  { if (!is.character(except))
      stop("except must be a character vector")
    if (!all(except %in% names(unr.names)))
      stop("except must contain the names of only unrandomized factors")
  }
#form vector of numbers of levels
  unr.levels <- rep(1, times=nunr)
  for (i in 1:nunr)
  { if (is.numeric(unr.names[[i]]) | is.character(unr.names[[i]]))
    { if (length(unr.names[[i]]) == 1)
        unr.levels[i] <- unr.names[[i]]
      else
        unr.levels[i] <- length(unr.names[[i]])
    }
    else
    { stop("Levels of factors must be specified using either numeric or character vectors")
    }
  }
  if(n != prod(unr.levels))
    stop("The product of the numbers of levels of the unrandomized factors ", 
         "must equal the length of the randomized factors.")
#process nested.factor argument
  if(is.null(nested.factors))
   facrand <- fac.rand.cross(unrandomized=unrandomized, unr.names=unr.names, 
                             unr.levels=unr.levels, nested.factors=nested.factors, 
                             except=except, seed=seed, randomized=randomized)
 else
  { if(!is.list(nested.factors))
      stop("nested.factors must be a list.")
    names.nested <- names(nested.factors)
    if(any(names.nested == ""))
      stop("all components of nested.factors must be named.")
    if(!all(sapply(nested.factors, FUN=is.character)) )
      stop("All elements of the nested.factors list must be of class character")
    facrand <- fac.rand.nest(unrandomized=unrandomized, unr.names=unr.names, 
                             unr.levels=unr.levels, nested.factors=nested.factors, 
                             except=except, seed=seed, randomized=randomized)
  }
  #form layout which is to be in standard or data frame order for the 
  #unrandomized factors supplied in unrandomized i.e. in unr.names order
  if (is.data.frame(unrandomized)) #get into unrandomized data.frame order
  { perm.derand <- do.call(order, facrand)
    perm.unrand <- do.call(order, unrandomized)
    #on the right of an expression perm.derand (perm.unrand)
    #puts facrand (unrandomized) into standard order.
    #on the left of an expression perm.derand (perm.unrand)
    #puts standard order into the same order as facrand (unrandomized),
    #i.e. a random permutation of standard order (data frame order).
    #put facrand into data.frame order and join with randomized factors
    for (i in 1:nunr)
      attributes(facrand[[i]]) <- attributes(unrandomized[[i]])
    faclay <- facrand
    faclay[perm.unrand, ] <- facrand
    faclay <- data.frame(faclay, randomized)
    #compute randomization for data.frame order
    perm.dat <- vector("numeric", length=n)
    #order to put permuted data frame into standard order
    perm.dat[perm.derand] <- perm.unrand
    #order to get from unrandomized data frame to permuted data frame
    perm.dat[perm.unrand] <- perm.dat
    perm.derand.dat <- vector("numeric", length=n)
    perm.derand.dat[perm.dat] <- 1:n
    if (unit.permutation)
      faclay <- data.frame(.Units = 1:n, .Permutation = order(perm.derand.dat),
                           faclay[perm.derand.dat, ])
    else
      faclay <- faclay[perm.derand.dat, ]
    
  }
  else  #get in standard order for unr.names
  { perm.derand <- do.call(order, facrand)
    faclay <- data.frame(facrand, randomized)
    if (unit.permutation)
      faclay <- data.frame(.Units = 1:n, .Permutation = order(perm.derand), 
                           faclay[perm.derand, ])
    else
      faclay <- faclay[perm.derand, ]
  }
  if (is.factor(randomized))
    names(faclay)[length(faclay)] <- deparse(substitute(randomized))
  rownames(faclay) <- 1:n
  faclay
}

