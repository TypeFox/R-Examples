######################
### global methods ###
######################
#' print \code{\link{dimVar-class}} objects
#'
#' print \code{\link{dimVar-class}} objects in a resonable way
#'
#' @aliases print,dimVar-method
#' @rdname print-method
#' @param x An object of class \code{\link{dimVar-class}}
#' @param ... currently not used
#' @export
#' @docType methods
setMethod(f="print", signature= "dimVar",
  definition=function (x, ...) {
    cat('dimVariable "',g_varname(x),'" has', g_nr_levels(x),'levels!\n')
    cat("the level-structure is given by:\n")
    cO <- g_original_codes(x)
    cD <- g_default_codes(x)
    lev <- g_levels(x)
    for ( i in 1:length(cO) ) {
      cat(paste(rep("-", 3*(lev[i]-1)), collapse=""),"> ", cO[i], "(",cD[i],")\n",sep="")
    }
    cat("--------------------------------\n")
    if ( g_has_dups(x) ) {
      cat("Note: redundant levels (",g_dups(x),") have been removed!\n")
    }
  }
)

#######################################
### methods only for class 'dimVar' ###
#######################################
#' @aliases get.dimVar,dimVar,character-method
#' @rdname get.dimVar-method
setMethod(f="get.dimVar", signature=c("dimVar", "character"),
  definition=function(object, type) {
    if ( !type %in% c("varName", "codesOriginal", "codesDefault",
        "codesMinimal", "levels", "structure", "dims", "dups",
        "dupsUp", "hasDuplicates", "nrLevels", "minimalCodesDefault") ) {
      stop("get.dimVar:: argument 'type' is not valid!\n")
    }
    if ( type == "varName" ) {
      return(g_varname(object))
    }
    if ( type == "codesOriginal" ) {
      return(g_original_codes(object))
    }
    if ( type == "codesDefault" ) {
      return(g_default_codes(object))
    }
    if ( type == "codesMinimal" ) {
      return(g_minimal_codes(object))
    }
    if ( type == "levels" ) {
      return(g_levels(object))
    }
    if ( type == "structure" ) {
      return(g_structure(object))
    }
    if ( type == "dims" ) {
      return(g_dims(object))
    }
    if ( type == "dups" ) {
      return(g_dups(object))
    }
    if ( type == "dupsUp" ) {
      return(g_dups_up(object))
    }
    if ( type == "hasDuplicates" ) {
      return(g_has_dups(object))
    }
    if ( type == "nrLevels" ) {
      return(g_nr_levels(object))
    }
    if ( type == "minimalCodesDefault" ) {
      return(g_minimal_default_codes(object))
    }
  }
)

#' @aliases calc.dimVar,dimVar,character,character-method
#' @rdname calc.dimVar-method
setMethod(f="calc.dimVar", signature=c("dimVar", "character", "character"),
  definition=function(object, type, input) {
    if ( !type %in% c("hasDefaultCodes", "matchCodeOrig", "matchCodeDefault",
        "standardize", "requiredMinimalCodes") ) {
      stop("calc.dimVar:: check argument 'type'!\n")
    }

    if ( type == "hasDefaultCodes" ) {
      return(c_has_default_codes(object, input))
    }

    if ( type == "matchCodeOrig" ) {
      return(c_match_orig_codes(object, input))
    }

    if ( type == "matchCodeDefault" ) {
      return(c_match_default_codes(object, input))
    }
    if ( type == "standardize" ) {
      return(c_standardize(object, input))
    }
    if ( type == "requiredMinimalCodes" ) {
      if ( length(input) != 1 ) {
        stop("requiredMinimalCodes:: length of argument 'code' must equal 1!\n")
      }
      if ( as.numeric(input) == 0 ) {
        dI <- g_default_codes(object)
        if ( length(dI) == 1 ) {
          out <- dI
        } else {
          out <- g_default_codes(object)[g_minimal_codes(object)==TRUE]
        }
      } else {
        isMinimal <- g_minimal_codes(object)[which(g_default_codes(object)==input)]
        if ( isMinimal ) {
          out <- input
        } else {
          dims <- g_dims(object)
          out <- NULL
          out <- dims[[max(which(!is.na(sapply(dims, function(x) { match(input, x )}))))]][-1]

          minCodes <- sapply(out, function(x) { g_minimal_codes(object)[which(g_default_codes(object)==x)] })

          if ( any(minCodes == FALSE) ) {
            runInd <- TRUE
            while(runInd) {
              checkInd <- which(minCodes==FALSE)
              removeVars <- names(minCodes[checkInd])
              for ( i in checkInd ) {
                new <- dims[[max(which(!is.na(sapply(dims, function(x) { match(out[i], x )}))))]][-1]
                out <- c(out, new)
              }
              out <- setdiff(out, removeVars)

              out <- unique(out)
              minCodes <- sapply(out, function(x) {
                g_minimal_codes(object)[which(g_default_codes(object)==x)]
              })
              if ( all(minCodes == TRUE) ) {
                runInd <- FALSE
              }
            }
          }
        }
      }
      out <- c_match_orig_codes(object, input=out)
      # add possible dups (recoding errors in rawData?? #
      dupsUp <- g_dups_up(object)
      if ( !is.null(dupsUp) ) {
        dups <- g_dups(object)
        ind <- which(dupsUp %in% out )
        if ( length(ind) > 0 ) {
          out <- c(out, dups[ind])
        }
      }
      names(out) <- NULL
    }
    return(out)
  }
)

#' @aliases init.dimVar,list-method
#' @rdname init.dimVar-method
setMethod(f='init.dimVar', signature=c('list'),
  definition=function(input) {
    vName <- input$vName
    input <- input$input

    calcInfo <- function(inputList) {
      genLevel <- function(strIDs, dimStructure) {
        strID <- level <- NULL
        cs <- c(0,cumsum(dimStructure))
        dt <- data.table(strID=strIDs, level=1)
        for ( i in 1:(length(cs)-1)) {
          cmd <- paste0("dt[substr(strID,",cs[i]+1,",",cs[i+1],")!='",paste0(rep(0,dimStructure[i]), collapse=""),"'")
          cmd <- paste0(cmd, ",level:=",i,"]")
          #cat(cmd,"\n")
          eval(parse(text=cmd))
        }
        dt[as.numeric(strID)==0, level:=1]
        if ( any(is.na(dt))) {
          stop("error while calculating the levels!\n")
        }
        return(as.integer(dt$level))
      }

      # define variables
      removeInd <- NULL

      ### calculate the levels and the the number of levels
      nrLevels <- length(unique(inputList$levels))

      ### calculate necessary digits to represent this hierarchy
      if ( nrLevels == 1 ) {
        nrDigits <- 1
      } else {
        nrDigits <- c(1, nchar(as.character(length(which(inputList$levels==2))))) # level 2
        if( nrLevels > 2 ) {
          for( i in 2:(nrLevels-1) ) {
            ss <- inputList$levels[which(inputList$levels %in% c(i, i+1))]
            fac <- rep(NA, length(ss))
            ind <- which(ss==i)
            for ( j in 1:(length(ind)) ) {
              if ( j != length(ind) )
                fac[ind[j]:(ind[j+1]-1)] <- j
              else
                fac[ind[j]:length(fac)] <- j
            }
            spl <- split(ss, fac)
            nrDigits <- c(nrDigits, nchar(as.character((max(unlist(lapply(spl, function(x) { length(x)})))-1))))
          }
        }
      }

      ### calculate standard-codes
      # calculate position of levels in standard codec (substrInd)
      if ( nrLevels == 1 ) {
        codes <- "0"
      } else {
        cs <- cumsum(nrDigits)
        substrInd <- list()
        substrInd[[1]] <- c(1,1)
        for( j in 2:nrLevels ) {
          whichDigits <- (cs[(j-1)]+1):cs[j]
          if( length(whichDigits) == 1 )
            whichDigits <- c(whichDigits,whichDigits)
          substrInd[[j]] <- whichDigits
        }
        codes <- rep(paste(rep("0", sum(nrDigits)), collapse=""), length(inputList$levels))
        # calc the standard codes
        for( i in 2:length(inputList$levels) ) {
          actLevel <- inputList$levels[i]
          charsActLevel <- nrDigits[actLevel]

          if( inputList$levels[i] >= inputList$levels[i-1] ) {
            oldInd <- i-1
            codes[i] <- codes[i-1]
            oldVal <- as.integer(substr(codes[oldInd],substrInd[[actLevel]][1],substrInd[[actLevel]][length(substrInd[[actLevel]])]))
            substr(codes[i], substrInd[[actLevel]][1], substrInd[[actLevel]][length(substrInd[[actLevel]])]) <- sprintf(paste("%0",charsActLevel,"d",sep=""),oldVal+1)
          } else if( inputList$levels[i] < inputList$levels[i-1] ) {
            # go back as far as necessary
            candidate <- which(inputList$levels==actLevel)
            oldInd <- candidate[max(which(candidate < i))]
            codes[i] <- codes[oldInd]
            oldVal <- as.integer(substr(codes[oldInd],substrInd[[actLevel]][1],substrInd[[actLevel]][length(substrInd[[actLevel]])]))
            substr(codes[i], substrInd[[actLevel]][1], substrInd[[actLevel]][length(substrInd[[actLevel]])]) <- sprintf(paste("%0",charsActLevel,"d",sep=""),oldVal+1)
          }
        }
        # calculate if a given level is neccessary
        for( i in 1:(length(inputList$levels)-1) ) {
          if( inputList$levels[i+1] > inputList$levels[i] )
            removeInd <- c(removeInd, i)
        }
      }

      codesMinimal <- rep(TRUE, length(inputList$levels))
      if( length(removeInd) > 0 )
        codesMinimal[removeInd] <- FALSE

      ### calculate additional information
      minInd <- codes[codesMinimal==TRUE]

      # calculate all possible characteristics (=sub|totals) for each dimensional variable
      newDims <- NULL
      for ( j in (length(nrDigits)-1):1 ) {
        spl <- split(minInd, substr(minInd, 1, sum(nrDigits[1:j])))
        spl <- lapply(spl, function(x) { as.character(unique(x))[1] } )
        for( z in 1:length(spl) ) {
          # calculating the upper limit and update the levels
          upperHier <- spl[[z]]
          from <- sum(nrDigits[1:j]) + 1
          to <- nchar(upperHier)
          substr(upperHier, from, to) <- paste(rep("0", (to - from + 1)), collapse="")
          if( !upperHier %in% as.character(codes) )
            newDims <- append(newDims, upperHier)
        }
      }

      # combine existing and possible new characteristics
      allDims <- unique(unlist(c(codes, newDims)))
      levels <- genLevel(allDims, nrDigits)
      dimensions <- list()
      if ( nrLevels == 1 ) {
        dimensions[[1]] <- codes
      } else {
        dat <- data.frame(dims=allDims, lev=levels)
        z <- 1
        for( i in 1:(max(dat$lev)-1) ) {
          tmp <- dat[which(dat$lev==i),]
          # top-Level
          if( nrow(tmp) == 1 ) {
            dimensions[[z]] <- as.character(dat[which(dat$lev %in% c(i, i+1)),"dims"])
            z <- z + 1
          }
          # splitting is necessary
          else {
            for( j in 1:nrow(tmp) ) {
              aktDim <- as.character(tmp[j,"dims"])
              aktLev <- tmp[j,"lev"]
              erg <- c(aktDim, as.character(dat[which(dat$lev==(aktLev+1) & substr(dat$dims,1,sum(nrDigits[1:aktLev])) == substr(aktDim,1,sum(nrDigits[1:aktLev]))), "dims"]))
              if( length(erg) > 1 ) {
                dimensions[[z]] <- erg
                z <- z + 1
              }
            }
          }
        }
        dimensions <- lapply(dimensions, sort)
      }

      # recalculate the neccessary (TRUE) and non-neccessary levels based on all possible levels (out$allDims)
      notUsed <- sort(unique(unlist(lapply(dimensions, function(x) x[1]))))
      codesMinimal <- rep(TRUE, length(allDims))
      codesMinimal[allDims %in% notUsed] <- FALSE

      out <- list(
          codesOrig=inputList$codes,
          codesDefault=codes,
          levelsOrig=inputList$levels,
          levelStructure=nrDigits,
          dimensions=dimensions,
          codesMinimal=codesMinimal)
      out
    }

    if ( is.data.frame(input) || is.matrix(input) ) {
      if ( ncol(input) > 2 ) {
        stop('input must only have 2 columns!\n')
      }
      if ( nchar(as.character(input[[1]][1])) != 1 ) {
        stop('"@" must be listed in first row and first column in input!\n')
      }
    } else {
      if ( !file.exists(input) ) {
        stop('check the path of input!\n')
      }
      input <- read.table(input, sep=";", dec=".", colClasses="character")
      if ( ncol(input) > 2 ) {
        stop('input must only have 2 columns!\n')
      }
      if ( nchar(as.character(input[[1]][1])) != 1 ) {
        stop('"@" must be listed in first row and first column in input!\n')
      }
      input <- as.data.table(input)
    }
    if ( sum(input[[1]]=="@") != 1 ) {
      stop(paste0('Error in input for dimension "',vName,'". There can be only one overall total coded with "@"!\n'))
    }
    inputList <- list()
    inputList$levels <- nchar(as.character(input[[1]]))
    inputList$codes <- as.character(input[[2]])

    # get complete level-structure
    infoComplete <- calcInfo(inputList)

    # search for duplicates
    dimLen <- sapply(infoComplete$dimensions, length)
    dups <- dupsUp <- NULL
    removeInd <- NULL
    if ( any(dimLen == 2) ) {
      index <- which(dimLen==2)
      for ( i in 1:length(index)) {
        indexInOrig1 <- match(infoComplete$dimensions[[index[i]]][2], infoComplete$codesDefault)
        levDiff <- setdiff(which(infoComplete$levels < infoComplete$levels[indexInOrig1]), 1:indexInOrig1)
        if ( length(levDiff) > 0 )
          indexInOrig2 <- min(levDiff)
        else
          indexInOrig2 <- length(inputList$levels)+1

        # move one level up
        if ( indexInOrig2 - indexInOrig1 > 1 ) {
          ind <- (indexInOrig1+1):(indexInOrig2-1)
          inputList$levels[ind] <- inputList$levels[ind]-1
        }
        # add info
        dups <- c(dups, infoComplete$codesOrig[indexInOrig1])
        dupsUp <- c(dupsUp, infoComplete$codesOrig[indexInOrig1-1])
        removeInd <- c(removeInd, indexInOrig1)
      }
      inputList$levels <- inputList$levels[-removeInd]
      inputList$codes <- inputList$codes[-removeInd]
    }

    info <- calcInfo(inputList)
    info$dups <- dups
    info$dupsUp <- dupsUp

    dimVar <- new("dimVar",
      codesOriginal=info$codesOrig,
      codesDefault=info$codesDefault,
      codesMinimal=info$codesMinimal,
      vName=vName,
      levels=info$levelsOrig,
      structure=info$levelStructure,
      dims=info$dimensions,
      dups=dups,
      dupsUp=dupsUp
    )
    return(dimVar)
  }
)

setMethod(f="g_varname", signature=c("dimVar"), definition=function(object) {
  return(object@vName)
})
setMethod(f="g_original_codes", signature=c("dimVar"), definition=function(object) {
  return(object@codesOriginal)
})
setMethod(f="g_default_codes", signature=c("dimVar"), definition=function(object) {
  return(object@codesDefault)
})
setMethod(f="g_minimal_codes", signature=c("dimVar"), definition=function(object) {
  return(object@codesMinimal)
})
setMethod(f="g_levels", signature=c("dimVar"), definition=function(object) {
  return(object@levels)
})
setMethod(f="g_structure", signature=c("dimVar"), definition=function(object) {
  return(object@structure)
})
setMethod(f="g_dims", signature=c("dimVar"), definition=function(object) {
  return(object@dims)
})
setMethod(f="g_dups", signature=c("dimVar"), definition=function(object) {
  return(object@dups)
})
setMethod(f="g_dups_up", signature=c("dimVar"), definition=function(object) {
  return(object@dupsUp)
})
setMethod(f="g_has_dups", signature=c("dimVar"), definition=function(object) {
  return(!is.null(g_dups(object)))
})
setMethod(f="g_nr_levels", signature=c("dimVar"), definition=function(object) {
  return (length(g_structure(object)))
})
setMethod(f="g_minimal_default_codes", signature=c("dimVar"),
definition=function(object) {
  codes_default <- g_default_codes(object)
  if (length(codes_default)==1) {
    return(codes_default)
  }
  return(g_default_codes(object)[g_minimal_codes(object)==TRUE])
})
setMethod(f="c_has_default_codes", signature=c("dimVar", "character"), definition=function(object, input) {
  out <- FALSE
  if ( !any(is.na(match(input, g_default_codes(object)))) ) {
    out <- TRUE
  }
  return(out)
})
setMethod(f="c_match_orig_codes", signature=c("dimVar", "character"), definition=function(object, input) {
  d <- g_default_codes(object)
  o <- g_original_codes(object)
  o[match(input, d)]
})

setMethod(f="c_match_default_codes", signature=c("dimVar", "character"), definition=function(object, input) {
  d <- g_default_codes(object)
  o <- g_original_codes(object)
  d[match(input, o)]
})

setMethod(f="c_standardize", signature=c("dimVar", "character"), definition=function(object, input) {
  ind <- which(is.na(match(input, g_original_codes(object))))
  if ( length(ind) > 0 ) {
    matchInd <- match(input[ind], g_dups(object))
    if ( any(is.na(matchInd)) ) {
      stop("c_standardize: elements of 'codesOrig' not listed in 'codesOriginal' or 'dups'!\n")
    }
    input[ind] <- g_dups_up(object)[matchInd]
  }
  out <- c_match_default_codes(object, input)

  if ( any(is.na(out)) ) {
    stop("c_standardize:: matching not successful!\n")
  }
  return(out)
})
setMethod(f="c_required_minimal_codes", signature=c("dimVar", "character"), definition=function(object, input) {
  if ( length(input) != 1 ) {
    stop("c_required_minimal_codes:: length of argument 'code' must equal 1!\n")
  }
  if ( as.numeric(input) == 0 ) {
    dI <- g_default_codes(object)
    if ( length(dI) == 1 ) {
      out <- dI
    } else {
      out <- g_default_codes(object)[g_minimal_codes(object)==TRUE]
    }
  } else {
    isMinimal <- g_minimal_codes(object)[which(g_default_codes(object)==input)]
    if ( isMinimal ) {
      out <- input
    } else {
      dims <- g_dims(object)
      out <- NULL
      out <- dims[[max(which(!is.na(sapply(dims, function(x) { match(input, x )}))))]][-1]

      minCodes <- sapply(out, function(x) { g_minimal_codes(object)[which(g_default_codes(object)==x)] })

      if ( any(minCodes == FALSE) ) {
        runInd <- TRUE
        while(runInd) {
          checkInd <- which(minCodes==FALSE)
          removeVars <- names(minCodes[checkInd])
          for ( i in checkInd ) {
            new <- dims[[max(which(!is.na(sapply(dims, function(x) { match(out[i], x )}))))]][-1]
            out <- c(out, new)
          }
          out <- setdiff(out, removeVars)
          out <- unique(out)
          minCodes <- sapply(out, function(x) {
            g_minimal_codes(object)[which(g_default_codes(object)==x)]
          })
          if ( all(minCodes == TRUE) ) {
            runInd <- FALSE
          }
        }
      }
    }
  }
  out <- c_match_orig_codes(object, input=out)
  # add possible dups (recoding errors in rawData?
  dupsUp <- g_dups_up(object)
  if ( !is.null(dupsUp) ) {
    dups <- g_dups(object)
    ind <- which(dupsUp %in% out )
    if ( length(ind) > 0 ) {
      out <- c(out, dups[ind])
    }
  }
  names(out) <- NULL
  return(out)
})

