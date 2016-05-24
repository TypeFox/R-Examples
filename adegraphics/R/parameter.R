changelatticetheme <- function(...) {
  ## lattice.options(default.theme = adegspecial)
  ## change lattice theme
  ## if a device is open, it would apply only to new devices
  if(try(is.list(...), silent = TRUE) == TRUE)
    changes <- as.list(...)
  else 
    changes <- list(...)
  
  newtheme <- get("lattice.theme", envir = getFromNamespace(".LatticeEnv", ns = "lattice")) 
  adegspecial <- get("adegtheme", envir = .ADEgEnv)
  
  if(length(changes))
    newtheme <- modifyList(newtheme, changes, keep.null = TRUE)
  else
    ## come back at the starting point
    newtheme <- modifyList(newtheme, adegspecial, keep.null = TRUE)   
  
  ## for all new devices
  lattice.options(default.theme = switch(EXPR = .Device, newtheme))
  if(dev.cur() != 1)  ## if a device is open
    trellis.par.set(newtheme)
  
  invisible(newtheme)
}


.mergingList <- function(tomerge) {
  ## merge elements of the list by name recurcively
  lnames <- names(tomerge)
  counter <- 0 ## safety counter
  while(length(lnames) != length(unique(lnames))) {
    ## be sure that there are duplicated values
    indix <- match(lnames, lnames)
    remove <- c()
    for(i in 1:length(indix)) {
      if(i != indix[i]) {
        tomerge[[indix[i]]] <- c(tomerge[[indix[i]]], tomerge[[i]])
        remove <- c(remove, i)
      }
    }
    if(length(remove))
      tomerge[remove] <- NULL
    tomerge <- lapply(tomerge, FUN = function(x) {if(is.list(x) & (length(x) > 1)) .mergingList(x) else x})
    counter <- counter + 1
    if(counter == 50) 
      stop("error in .mergingList", call. = FALSE)
    lnames <- names(tomerge)
  }
  return(tomerge)
}


.replaceList <- function(x, val) {
  ## replaceList: inspired by modifyList but 
  ## replace only previous existing elements and with partial names matching
  ## x: list to modify, val: modications to pass
  ## x structure can not be changed
  ## To be more specific if an element is a list, it cannot be change with a single value
  rest <- list()  
  returned <- list()
  xnames <- names(x)
  
  for (v in names(val)) {
    indix <- pmatch(v, xnames, nomatch = 0)
    if(indix > 0) {
      ## if there is a match
      if(is.list(x[[indix]]) && (!is.list(val[[v]])))
        stop(paste("cannot replace a list: ", xnames[indix], " by a single value element", sep = ""), call. = FALSE)
      else {
        if(is.list(x[[indix]])) {
          ## recursivity
          replace <- .replaceList(x[[indix]], val[[v]])
          returned <- c(returned, list(replace$select))
          rest <- c(rest, replace$rest)
        }
        else 
          returned[[(length(returned) + 1)]] <- val[[v]] ## else replace values
        
        names(returned)[length(returned)] <- xnames[indix]
      }
    }
    else rest <- c(rest, val[v])
  }
  return(list(select = returned, rest = rest))
}


.getlist <- function(keys, values) {
  ## assembles keys and values as list of list
  ## keys: list of characters vectors, the keys splitted, values: the original list
  result <- list()
  for(i in 1:length(keys)) {
    l <- list(values[[i]])
    names(l) <- keys[[i]][length(keys[[i]])] 
    if(length(keys[[i]]) > 1)
      for(j in (length(keys[[i]]) - 1):1) {
        l <- list(l)
        names(l)[1] <- keys[[i]][j]
      }
    result[[i]] <- l
  }
  return(result)
}


separation <- function(... , pattern = 0, split = TRUE) {
  ## separate between the list passed to the function and the one already known
  ## if pattern is 1, compare to trellis parameters
  ## if pattern is 0, compare to 'padegraphic' parameters
  
  ## gets dots
  if(try(is.list(...), silent = TRUE) == TRUE)
    tmp_list <- as.list(...)
  else
    tmp_list <- list(...)
  if(is.null(names(tmp_list)))
    names(tmp_list) <- tmp_list
  if(!length(tmp_list))
    return(list(select = list(), rest = list()))
  
  ## get pattern 
  if(is.list(pattern))
    listpat <- pattern
  else {
    if(pattern > 1)
      stop("pattern must be 0 or 1 in 'separation' function", call. = FALSE)
    else {
      if(pattern == 1)
      	listpat <- trellis.par.get()
      else{
          listpat <- get("padegraphic", envir = .ADEgEnv)
      }
    }
  }
  ## splitting list keys
  if(!is.list(pattern)) {
    if(pattern != 1 && split) {
      ## adegpar, collates keys
      sep <- strsplit(names(tmp_list), split = ".", fixed = TRUE)
      values <- tmp_list
      val <- .getlist(keys = sep, values = values) ## assemblies keys with values, as list of list...
      val <- sapply(val, FUN = function(x) return(x))
      val <- .mergingList(val)
    } else
    	val <- tmp_list
  } else
    val <- tmp_list
  res <- .replaceList(x = listpat, val)
  res[[1]] <- .mergingList(res[[1]])
  return(res)
}


adegpar <- function(...) {
  ## case 0: nothing in parenthesis
  ## case 0 bis: only one key (no indication sublist, "paxes")
  ## case 1: ...= "axes.draw", "sub", "sub.size" # one level, only names
  ## case 2: ...= "axes" = list("draw") # two levels, only names
  ## case 3: ...= "axes.draw" = FALSE, "sub.size" = 12 # one level, key names and matching values
  ## case 4: ...= axes=list(draw=TRUE), sub=list(size=55) # two levels, key names and matching values
  ## case 5 : ... is a complete list
  ## if ... is a list
  
  ## does not assign, only find corresponding element in list patti
  recursfinder <- function(x, patti) {
    result <- list()
    
    okfu <- function(x, patti) {
      ## okfu: retrieve good values and keys (patti)
      if(length(x) > 1) 
        stop("x has length > 1") ## to remove
      idx <- pmatch(names(x), names(patti))
      if(!is.na(idx))
        return( patti[idx])
      else    
        return(NA)
    }
    if(!is.list(x[[1]]))
      result <- c(result, okfu(x, patti))
    else {
      idx <- pmatch(names(x), names(patti))
      if(!is.na(idx)) {
        result <- c(result, list(recursfinder(x[[1]], patti[[idx]])))
        names(result) <- names(patti[idx])
      }
      else
        print("no matching found in adegpar")
    }
    return(result)
  } ## end recurs finder
  
  nonames <- function(userlist, pattili) {
    ## return the right values list
    sep <- sapply(userlist, strsplit, split = ".", fixed = TRUE) ## a list
    values <- userlist
    val <- .getlist(keys = sep, values = values)
    return(sapply(val, FUN = recursfinder, patti = pattili)) ## get result
  }
  
  value <- list()
  assignement <- FALSE
  if(try(is.list(...), silent = TRUE) == TRUE)
    argu <- as.list(...) ## ... is still a list
  else
    argu <- list(...) ## tranforms in list
  ## choose option
  padegr <- get("padegraphic", envir = .ADEgEnv)
  
  ## switching case: recursive
  switchingcase <- function(userlist, patternlist) {
    if(!length(userlist)) ## empty case 0
      return(list(result = patternlist, assigni = list()))
    else {
      lnames <- names(userlist)
      if(is.null(lnames)) { ## no values, case 0 bis or 1
        res <- nonames(userlist, patternlist)
        return(list(result = res, assigni = list()))
      } else {
        result <- list()
        assigni <- list() ## initialization
        for(i in 1:length(lnames)) {
          if(identical(lnames[i], "")) {
            ## no names, meaning value is the key cas 2/1
            result <- c(result,nonames(userlist[i], patternlist))
          } else {
            ## we have names so value to assign, or sublist
            sep <- sapply(lnames[i], strsplit, split = ".", fixed = TRUE) ## a list
            ## get a list of list with right keys (splitting *.* keys)
            val <- .getlist(keys = sep, values = userlist[i])[[1L]]
            
            idx <- pmatch(names(val), names(patternlist))
            if(!is.na(idx)) {
              ## match with patternlist
              if(is.list(val[[1]])) {
                ## sublist val from user list
                ok <- switchingcase(val[[1]], patternlist = patternlist[[idx]])
                if(length(ok$result)) {
                  result <- c(result, list(ok$result))
                  ## level behind
                  names(result)[length(result)] <- names(patternlist[idx])
                }
                if(length(ok$assigni)) {
                  assigni <- c(assigni, list(ok$assigni))
                  names(assigni)[length(assigni)] <- names(patternlist[idx])
                }
              } else {
                ## if not a list, then a value to assign
                if(is.list(patternlist[[idx]]))
                  stop(paste("be careful, intent to replace in adegraphics parameters: ", names(patternlist[idx]), " by a single value element", sep = ""), call. = FALSE)
                assigni <- c(assigni, list(userlist[[i]]))
                names(assigni)[length(assigni)] <- names(patternlist[idx])
              }
            }
          }
        }
      }
    }
    return(list(result = result, assigni = assigni))
  } ## end switching case
  if(!length(argu)) ## ... empty
    return(padegr) ## case 0
  else {
    ## adegpar called with arguments
    switchi <- switchingcase(argu, padegr)
  }
  if(length(switchi$assign)) {
    padegr <- modifyList(padegr, switchi$assign, keep.null = TRUE)
    assign("padegraphic", padegr, envir = .ADEgEnv)
    return(invisible(padegr)) ## must be improve : avoid two calls to padegraphic
  }
  return(switchi$result)  
}
