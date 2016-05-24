"combined.split" <- function(x)
{ terms <- strsplit(x, "[& ]")[[1]]
  return(terms)
}

"findDecomposition" <- function(icomb, term, source, object, ncomb)
#Looks for the combination of term confounded with source in object 
#that consists of ncomb decompositions. Starts with decompistion icomb
{ while (!(term %in% names(object[[icomb]][[source]])))
  { icomb <- icomb + 1
    if (icomb > ncomb)
      stop(paste("Could not find",term,"when confounded with",source," "))
  }
  return(icomb)
}

"findLastwithSource" <- function(icomb, source, object, ncomb)
  #Looks for the last combination that contains source
{ while (!(source %in% names(object[[icomb]])))
  { icomb <- icomb + 1
    if (icomb > ncomb)
      stop(paste("Could not find",source," "))
  }
  #Check for later sources
  repeat
  { if (!(source %in% names(object[[icomb+1]])))
      break
    icomb <- icomb + 1
  }
  return(icomb)
}

"makeCombinedSource" <- function(terms)
{ source <- paste(terms, collapse="&")
  return(source)     
}

"replace.proj" <- function(object)
  #Replace projectors with their degrees of freedom
{ if (!inherits(object, "p2canon"))
  stop("object must be of class p2canon as produced by proj.2canon")
  
  #Loop through p2canon object
  Q1labels <- names(object)
  efficiencies <- vector(mode = "list", length = 0)
  for (i in Q1labels)
  { object[[i]][["Q1res"]] <- degfree(object[[i]][["Q1res"]])
    Q2labels <- names(object[[i]])[-1]
    if (length(Q2labels) > 0)
    { for (j in Q2labels)
      object[[i]][[j]][["Qproj"]] <- degfree(object[[i]][[j]][["Qproj"]])
    }
  }
  return(object)
}

"projs.canon" <- function(formulae, orthogonalize = "differencing", 
                          which.criteria = c("aefficiency", "eefficiency", "order"), 
                          omit.projectors = c("p2canon", "combined"), data = NULL, ...)
{ #examine the relationships between the sets of mutually orthogonal projection matrices for the supplied formulae
  
  #Check arguments and intitialize
  if (!is.list(formulae))
    stop("formulae must be a list")
  else
  { ntiers <- length(formulae)
    ncomb <- ntiers - 1 
    if (!all(unlist(lapply(formulae, inherits, what="formula"))))
      stop("formulae must contain a list of formulae")
  }
  options <- c("differencing", "eigenmethods")
  if (length(orthogonalize) == 1)
    orthogonalize <- rep(orthogonalize, ntiers)
  else
    if (length(orthogonalize) != ntiers)
    { warning("Length of orthogonalize is not equal to 1 or the number of formulae - only using first value")
      orthogonalize <- rep(orthogonalize[1], ntiers)
    }
  orthog <- options[unlist(lapply(orthogonalize, check.arg.values, options=options))]
  
  #check which.criteria arguments
  criteria <- c("aefficiency", "eefficiency", "mefficiency", "sefficiency", "order", "dforthog")
  options <- c(criteria, "none", "all")
  kcriteria <- options[unlist(lapply(which.criteria, check.arg.values, 
                                     options=options))]
  if ("all" %in% kcriteria)
    kcriteria <- criteria
  anycriteria <- !("none" %in% kcriteria)
  
  options <- c("p2canon", "combined", "none")
  omit.proj <- options[unlist(lapply(omit.projectors, check.arg.values, 
                                     options=options))]
  if ("none" %in% omit.proj)
    omit.proj <- "none"
  
  if (is.null(data) | !is.data.frame(data))
    stop("Must supply a data.frame for data")
  Q <- vector(mode = "list", length=ntiers)
  
  #Get set of projectors for the first formula
  if (ntiers == 1)
    Q[[1]] <- projs.structure(formulae[[1]], which.criteria = kcriteria, orthogonalize = orthog, 
                              data=data, ...)
  else 
    Q[[1]] <- projs.structure(formulae[[1]], which.criteria = kcriteria, orthogonalize = orthog[1], 
                              data=data, ...)
  
  #Loop over remaining  formulae
  CombinedSets <- vector(mode="list", length=ntiers)
  CombinedSets[[ntiers]] <- Q[[1]]
  if (ncomb > 0)
  { for (k in 1:ncomb)
    { ktier <- k + 1
      Q[[ktier]] <- projs.structure(formulae[[ktier]], which.criteria = kcriteria, 
                                    orthogonalize = orthog[ktier], data=data, ...)
      CombinedSets[[k]] <- projs.2canon(CombinedSets[[ntiers]], Q[[ktier]])
      CombinedSets[[ntiers]] <- projs.combine.p2canon(CombinedSets[[k]])
      if ("p2canon" %in% omit.proj)
        CombinedSets[[k]] <-  replace.proj(CombinedSets[[k]])
    }
    if ("combined" %in% omit.proj)
    { comb.labels <- names(CombinedSets[[ntiers]])
      for (i in comb.labels)
        CombinedSets[[ntiers]][[i]] <- degfree(CombinedSets[[ntiers]][[i]])
    }
  }
  class(CombinedSets) <- "pcanon"
  return(CombinedSets)
}

"summary.pcanon" <- function(object, which.criteria = c("aefficiency", "eefficiency", "order"), ...)
  #Routine to output a summary of the projector analysis in a pcanon object
{ if (!inherits(object, "pcanon"))
  stop("object must be of class pcanon as produced by projs.pcanon")
  #check which.criteria arguments
  criteria <- c("aefficiency", "eefficiency", "mefficiency", "sefficiency", "order", "dforthog")
  options <- c(criteria, "none", "all")
  kcriteria <- options[unlist(lapply(which.criteria, check.arg.values, 
                                     options=options))]
  if ("all" %in% kcriteria)
    kcriteria <- criteria
  anycriteria <- !("none" %in% kcriteria)
  
  ntiers <- length(object)
  ncomb <- ntiers - 1
  #Check if have projector or df in pcanon object
  have.proj <- FALSE
  if (ntiers == 1)
  { if (inherits(object[[1]], "projector"))
    have.proj <- TRUE
  }  
  else 
    if (inherits(object[[1]][[1]][["Q1res"]], "projector"))
      have.proj <- TRUE
  
  #Form data frame for summary table
  #Initialize
  sources <- names(object[[ntiers]])
  nlines <- length(sources)
  terms <- lapply(sources, combined.split)
  kcomb <- rep(1, nlines)
  ktiers <- lapply(terms, length)
  nc <- ntiers*2
  if (ntiers == 1)
    srcdf.names <- c("Source", "df")
  else
    srcdf.names <- as.vector(outer(c("Source", "df"), as.character(1:ntiers), paste, sep=""))
  orthogonaldesign <- TRUE
  if (anycriteria & ntiers > 1)
  { nc <- nc + length(kcriteria)
  res.criteria <- vector(mode = "list", length = length(kcriteria))
  names(res.criteria) <- kcriteria
  res.criteria[kcriteria] <- NA
  summary <- data.frame(matrix(nrow = nlines, ncol=nc))
  colnames(summary) <- c(srcdf.names, kcriteria)
  }
  else
  {   summary <- data.frame(matrix(nrow = nlines, ncol=nc))
  colnames(summary) <- srcdf.names
  }
  #Loop over sources
  for (kl in 1:nlines)
    if (ntiers ==1)
    { #Get first tier source and df
      summary[[1]][kl] <- terms[[kl]][1]
      summary[[2]][kl] <- degfree(object[[1]][[terms[[kl]][1]]])
    } else
    { #If next term is a residual, check it is for the previous term and output analogously
      if (terms[[kl]][ktiers[[kl]]] == "Residual")
      { #Check same as last term
        if (any(terms[[kl]][1:(ktiers[[kl]]-1)] != terms[[kl-1]][1:(ktiers[[kl]]-1)]))
          stop("Residual term and previous source are not confounded with the same source")
        label <- makeCombinedSource(terms[[kl]][1:(ktiers[[kl]]-1)])
        kcomb[kl] <- findLastwithSource(kcomb[kl], label, object, ncomb)
        summary[kl, 1:(kcomb[kl]*2)] <- summary[(kl-1), 1:(kcomb[kl]*2)]
        summary[kl, (kcomb[kl]*2+1)] <- "Residual"
        if (have.proj)
          summary[[kcomb[kl]*2+2]][kl] <- degfree(object[[kcomb[kl]]][[label]][["Q1res"]])
        else
          summary[[kcomb[kl]*2+2]][kl] <- object[[kcomb[kl]]][[label]][["Q1res"]]
        if (ktiers[[kl]]  == 1)
          stop("Have a Residual in the first formula - include a factor for the units")
        else
        { lastNonResidual <- max(c(1:ktiers[[kl]])["Residual" != terms[[kl]]])
        if (lastNonResidual  > 1)
        { pcomb <- 1
        label <- makeCombinedSource(terms[[kl]][1:(lastNonResidual-1)])
        term2 <- terms[[kl]][lastNonResidual]
        pcomb <- findDecomposition(pcomb, term2, label, object, ncomb) 
        if (abs(1 - unlist(object[[pcomb]][[label]][[term2]][["adjusted"]]["aefficiency"])) > 1e-04)
          orthogonaldesign <- FALSE
        if (anycriteria)
          summary[kl, kcriteria] <- unlist(object[[pcomb]][[label]][[term2]][["adjusted"]][kcriteria])
        }
        }
        #      if (anycriteria)
        #      { if (ktiers[[kl]]  == 1)
        #          stop("Have a Residual in the first formula - include a factor for the units")
        #        else
        #        { lastNonResidual <- max(c(1:ktiers[[kl]])["Residual" != terms[[kl]]])
        #          if (lastNonResidual  > 1)
        #          { pcomb <- 1
        #            label <- makeCombinedSource(terms[[kl]][1:(lastNonResidual-1)])
        #            term2 <- terms[[kl]][lastNonResidual]
        #            pcomb <- findDecomposition(pcomb, term2, label, object, ncomb) 
        #            summary[kl, kcriteria] <- unlist(object[[pcomb]][[label]][[term2]][["adjusted"]][kcriteria])
        #          }
        #        }
        #      }
      }
      else  #Not a Residual
      { #Get first tier source and df
        summary[[1]][kl] <- terms[[kl]][1]
        kconf.sources <- names(object[[1]][[terms[[kl]][1]]])[-1]
        if (have.proj)
        { totdf <- sum(unlist(lapply(kconf.sources, 
                                     function(ksrc, obj){ degfree(obj[[ksrc]][["Qproj"]])}, 
                                     obj = object[[1]][[terms[[kl]][1]]])))
        totdf <- totdf + degfree(object[[1]][[terms[[kl]][1]]][["Q1res"]])
        }
        else
        { totdf <- sum(unlist(lapply(kconf.sources, 
                                     function(ksrc, obj){ obj[[ksrc]][["Qproj"]]}, 
                                     obj = object[[1]][[terms[[kl]][1]]])))
        totdf <- totdf + object[[1]][[terms[[kl]][1]]][["Q1res"]]
        }
        summary[[2]][kl] <- totdf
        label <- summary[[1]][kl]
        #Add remaining tiers
        if (ktiers[[kl]] > 1)
        { for (kterm in 2:ktiers[[kl]])
        { term2 <- terms[[kl]][kterm]
        if (term2 == "Residual")
        { if (have.proj)
          summary[[kcomb[kl]*2+2]][kl] <- degfree(object[[kcomb[kl]]][[label]][["Q1res"]])
        else
          summary[[kcomb[kl]*2+2]][kl] <- object[[kcomb[kl]]][[label]][["Q1res"]]
        }
        else
        { kcomb[kl] <- findDecomposition(kcomb[kl], term2, label, object, ncomb)
        if (have.proj)
          summary[[kcomb[kl]*2+2]][kl] <- degfree(object[[kcomb[kl]]][[label]][[term2]][["Qproj"]])
        else
          summary[[kcomb[kl]*2+2]][kl] <- object[[kcomb[kl]]][[label]][[term2]][["Qproj"]]
        }
        summary[[kcomb[kl]*2+1]][kl] <- term2
        if (kterm != ktiers[[kl]])
        { label <- paste(label,term2,sep="&")
        kcomb[kl] <- kcomb[kl] + 1
        }
        }
          if (abs(1 - unlist(object[[kcomb[kl]]][[label]][[term2]][["adjusted"]]["aefficiency"])) > 1e-04)
            orthogonaldesign <- FALSE
        }
        if (anycriteria & ktiers[[kl]] > 1)
          summary[kl, kcriteria] <- unlist(object[[kcomb[kl]]][[label]][[term2]][["adjusted"]][kcriteria])
      }
    }
  class(summary) <- c("summary.pcanon", "data.frame")
  if (ntiers == 1)
    attr(summary, which = "title") <- 
    "\n\nSummary table of the decomposition\n\n"
  else
    attr(summary, which = "title") <- 
    "\n\nSummary table of the decomposition (based on adjusted quantities)\n\n"
  attr(summary, which = "ntiers") <- ntiers
  attr(summary, which = "orthogonal") <- orthogonaldesign
  return(summary)
}

print.summary.pcanon <- function(x, ...)
{ if (!inherits(x, "summary.pcanon"))
  stop("Must supply an object of class summary.p2canon")
  cat(attr(x, which="title"))
  ntiers <- attr(x, which="ntiers")
  y <- x
  nlines <- nrow(y)
  if (ntiers ==1)
  { dffw = max(3, floor(log10(max(y$df))) + 1)
    y$df <- formatC(y$df, format="f", digits=0, width=dffw)
  } else
    dffw = max(3, floor(log10(max(y$df1))) + 1)
  repeats <- vector(mode="list", length=ntiers)
  if (ntiers > 1)
  { for (ktier in 1:ntiers)
    { src.name <- paste("Source", ktier, sep="")
      repeats[[ktier]] <- c(FALSE, y[2:nlines, src.name] == y[1:(nlines-1),src.name])
      repeats[[ktier]][is.na(repeats[[ktier]])] <- FALSE
      if (ktier > 1)
        repeats[[ktier]] <- (repeats[[ktier]] & repeats[[(ktier-1)]])
      y[repeats[[ktier]], src.name] <- "  "
      df.name <- paste("df", ktier, sep="")
      y[[df.name]] <- formatC(y[[df.name]], format="f", digits=0, width=dffw)
      y[repeats[[ktier]], df.name] <- "  "
      y[[df.name]] <- gsub("NA", "  ", y[[df.name]])
    }
    if ("aefficiency" %in% names(y))
    { y$aefficiency <- formatC(y$aefficiency, format="f", digits=4, width=11)
      y$aefficiency <- gsub("NA", "  ", y$aefficiency)
    }
    if ("mefficiency" %in% names(y))
    { y$mefficiency <- formatC(y$mefficiency, format="f", digits=4, width=11)
      y$mefficiency <- gsub("NA", "  ", y$mefficiency)
    }
    if ("eefficiency" %in% names(y))
    { y$eefficiency <- formatC(y$eefficiency, format="f", digits=4, width=11)
      y$eefficiency <- gsub("NA", "  ", y$eefficiency)
    }
    if ("sefficiency" %in% names(y))
    { y$sefficiency <- formatC(y$sefficiency, format="f", digits=4, width=11)
      y$sefficiency <- gsub("NA", "  ", y$sefficiency)
    }
    if ("order" %in% names(y))
    { y$order <- formatC(y$order, format="f", digits=0, width=5)
      y$order <- gsub("NA", "  ", y$order)
    }
    if ("dforthog" %in% names(y))
    { y$dforthog <- formatC(y$dforthog, format="f", digits=0, width=8)
      y$dforthog <- gsub("NA", "  ", y$dforthog)
    }
  }
  print.data.frame(y, na.print="  ", right=FALSE, row.names=FALSE)
  if (!attr(x, which="orthogonal"))
    cat("\nThe design is not orthogonal\n\n")
  invisible(x)
}

"efficiencies.pcanon" <- function(object, which = "adjusted")
  #function to extract the efficiency factors from a pcanon object
{ if (!inherits(object, "pcanon"))
    stop("must supply an object of class pcanon as produced by projs.pcanon")
  options <- c("adjusted", "pairwise")
  opt <- options[check.arg.values(which, options)]
  
  #Get efficiencies
  ntiers <- length(object)
  efficiencies <- vector(mode="list", length=(ntiers-1))
  for (k in 1:(ntiers-1))
     efficiencies[[k]] <- efficiencies.p2canon(object[[k]], which=opt)
  return(efficiencies)
}  

