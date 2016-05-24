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
    { kfac <- kfac+1
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

"fac.gen" <- 
function(generate, each=1, times=1, order="standard")
{
#generate is a list of factors and numbers that specify the pattern in which 
#the levels of the factors are to be generated.
#If the component of the list is a factor name, it should be the name of a list 
#that contains either a single numeric value that is the number of levels, a 
#numeric vector that contains the levels of the factor or a character vector 
#that contains the labels of the factor.
  if(!is.list(generate))
    stop("generate must be a list")
  facs.reps <- factor.list(generate, order)
  fnames <- facs.reps$factors
  freps <- facs.reps$reps
  nfac <- length(fnames)
  levels <- rep(1, times=nfac)
  for (i in 1:nfac)
  { if (is.numeric(fnames[[i]]) | is.character(fnames[[i]]))
    { if (length(fnames[[i]]) == 1)
        levels[i] <- fnames[[i]]
      else
        levels[i] <- length(fnames[[i]])
    }
    else
    { stop("Levels of factors must be specified using either numeric or character vectors")
    }
  }
  n <- prod(levels)*prod(freps)*each*times
	which.ord <- pmatch(casefold(order), c("standard", "yates"), nomatch="")
	if (which.ord == "")	stop("order must be either standard or yates")
# standard order
	if (which.ord == "1") 
    counter <- nfac:1
  else
# Yates order
    counter <- 1:nfac
  genlist <- vector("list", nfac)
  keach <- each
  for (i in counter)
  { lev <- 1:levels[i]
    keach <- keach*freps[i]
    ktimes <- n/(levels[i]*keach)
    { if (is.numeric(fnames[[i]]))
      { if (length(fnames[[i]]) != 1)
          lev <- fnames[[i]]
        genlist[[i]] <- factor(rep(lev, times=ktimes, each=keach))   
      }
      else
      { genlist[[i]] <- factor(rep(lev, times=ktimes, each=keach), labels=fnames[[i]])
      }
    keach <- keach*levels[i] 
    }
  }
  genframe <- data.frame(genlist)
  names(genframe) <- names(fnames) 
 	genframe
}

