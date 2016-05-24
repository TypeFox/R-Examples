"fac.getinTerm" <- function(term)
  #function to return the set of factors/variables in a term separated by ':"
{ unlist(strsplit(term, ":", fixed=TRUE))}

"projs.jandw" <- function(R, Q, which.criteria = c("aefficiency","eefficiency","order"))
#A function to use J&W to orthogonalise the set of Q to R and previous Q
{ if (!is.list(Q))
    stop("The matrices to orthogonalize to must be in a list")

  #check which.criteria arguments
  criteria <- c("aefficiency", "eefficiency", "mefficiency", "sefficiency", "order", "dforthog")
  options <- c(criteria, "none", "all")
  kcriteria <- options[unlist(lapply(which.criteria, check.arg.values, 
                                     options=options))]
  if ("all" %in% kcriteria)
    kcriteria <- criteria
  anycriteria <- !("none" %in% kcriteria)
  #set up efficiency summary
  if (anycriteria)
  { nc <- 2 + length(kcriteria)
    summary <- data.frame(matrix(nrow = 0, ncol=nc))
    colnames(summary) <- c("Source", "df", kcriteria)
    eff.crit <- vector(mode="list",length=length(kcriteria))
    names(eff.crit) <- kcriteria
  }

  #do orthogonalization
  terms <- names(Q)
  for (i in 1:length(terms))
  { decomp <- proj2.combine(R, Q[[terms[i]]])
    Q[[terms[i]]] <- decomp$Qconf
    R <- decomp$Qres
    df <- degfree(Q[[terms[i]]])
    if (df == 0)
    { warning(paste(terms[[i]],"is aliased with previous terms in the formula", sep=" "))
      if (anycriteria)
      { eff.crit <- 0
        summary <- rbind(summary, 
                         data.frame(c(list(Source = terms[[i]], df = df), eff.crit), 
                                    stringsAsFactors = FALSE))
      }
    } else
    { keff.crit <- efficiency.criteria(decomp$efficiencies)
      if ((df - keff.crit[["dforthog"]]) != 0)
      { warning(paste(terms[[i]],"is partially aliased with previous terms in the formula", sep=" "))
        if (anycriteria)
          summary <- rbind(summary, 
                           data.frame(c(list(Source = terms[[i]], df = df), keff.crit[kcriteria]), 
                                      stringsAsFactors = FALSE))
      }
    }
  }
  
  #Print out the efficiency criteria if which.criteria is set
  if (anycriteria & nrow(summary) > 0)
  { cat("\nTable of efficiency criteria for aliasing between terms within a structure\n\n")
    print(summary)
  }
  return(Q)
}
  
"projs.structure" <- function(formula, orthogonalize = "differencing", 
                              which.criteria = c("aefficiency","eefficiency","order"), 
                              data = NULL, ...)
{ #generate a set of mutually orthogonal projection matrices, one for each term in the formula
  options <- c("differencing", "eigenmethods")
  opt <- options[check.arg.values(orthogonalize, options)]
  #check which.criteria arguments
  criteria <- c("aefficiency", "eefficiency", "mefficiency", "sefficiency", "order", "dforthog")
  options <- c(criteria, "none", "all")
  kcriteria <- options[unlist(lapply(which.criteria, check.arg.values, 
                                     options=options))]
  if ("all" %in% kcriteria)
    kcriteria <- criteria
  anycriteria <- !("none" %in% kcriteria)
  
  #Initialize
  if (is.null(data) | !is.data.frame(data))
    stop("Must supply a data.frame for data")
  n <- nrow(data)
  Q.G = projector(matrix(1, nrow=n, ncol=n)/n)
  fac.modl <- model.frame(formula, data=data)
  
  #get terms and form mean operators for each term
  terms <- attr(terms(formula, ...), which="term.labels")
  Q <- vector("list", length=length(terms))
  names(Q) <- terms
  for (k in 1:length(terms))
  { Q[[terms[k]]] <- model.matrix(as.formula(paste("~ ",terms[k])), data=fac.modl)
    Q[[terms[k]]] <- Q[[terms[k]]] %*% ginv(t(Q[[terms[k]]]) %*% Q[[terms[k]]]) %*% t(Q[[terms[k]]])
    Q[[terms[k]]] <- projector(Q[[terms[k]]] - Q.G)
  }
  
  #form projection matrices of the structure
  if (length(terms) > 1)
  { if (orthogonalize == "differencing") #by difference
    { orthogonal <- TRUE
      fac.mat <- attr(terms(formula), which="factors")
      for (i in 2:length(terms))
      { Q.work <- Q[[terms[i]]]
        for (j in 1:length(terms))
        { if (i != j)
          { if (all((fac.mat[,j] != 0) == (fac.mat[,j] & fac.mat[,i])))
               Q.work <- Q.work - Q[[terms[j]]]
          }
        }
        Q.work <- projector(Q.work)
        if (degfree(Q.work) == 0)
          warning(paste(terms[[i]],"is aliased with previous terms in the formula", sep=" "))
        else #Check that this term is orthogonal to previous projectors
        { i1 <- i - 1 
          if (i1 > 0)
            for (j in 1:i1)
              if (!is.allzero(Q.work %*% Q[[terms[j]]]))
              { warning(paste("** Projection matrices for ",terms[i], " and ", terms[j], 
                              " are not orthogonal", sep=""))
              }
        }
        Q[[terms[i]]] <- Q.work
      }
    }
    else  #by recursive orthogonalization
    { R <- projector(diag(1, nrow = n, ncol = n) - Q[[terms[1]]] - Q.G)
      Q[terms[2:length(terms)]] <- projs.jandw(R, Q[terms[2:length(terms)]],
                                               which.criteria = kcriteria)
    }
  }
  return(Q)
}    

"efficiency.criteria" <- function(efficiencies)
{ daeTolerance <- get("daeTolerance", envir=daeEnv)
  criteria <- vector(mode="list", length = 6)
  names(criteria) <- c('aefficiency','mefficiency','sefficiency','eefficiency','order',"dforthog")
  df.orthog <- sum(abs(1-efficiencies) <  daeTolerance[["eigen.tol"]])
  eff.unique <- remove.repeats(efficiencies, daeTolerance[["eigen.tol"]])
  K <- length(eff.unique)
  if (K == 1)
  { if (eff.unique == 0)
    { criteria["aefficiency"] <- 0
      criteria["mefficiency"] <- 0
      criteria["sefficiency"] <- 0
      criteria["eefficiency"] <- 0
      criteria["order"] <- 0
      criteria["dforthog"] <- 0
  }
    else
    { criteria["aefficiency"] <- eff.unique
      criteria["mefficiency"] <- eff.unique
      criteria["sefficiency"] <- 0
      criteria["eefficiency"] <- eff.unique
      criteria["order"] <- K
      criteria["dforthog"] <- df.orthog
    }
  }
  else
  { criteria["aefficiency"] <- harmonic.mean(efficiencies)
    criteria["mefficiency"] <- mean(efficiencies)
    criteria["sefficiency"] <- var(efficiencies)
    criteria["eefficiency"] <- min(efficiencies)
    criteria["order"] <- K
    criteria["dforthog"] <- df.orthog
  }
  return(criteria)
}

print.summary.p2canon <- function(x, ...)
{ if (!inherits(x, "summary.p2canon"))
    stop("Must supply an object of class summary.p2canon")
  cat(attr(x, which="title"))
  y <- x
  nlines <- nrow(y)
  repeats <- c(FALSE, y[2:nlines,"Source"] == y[1:(nlines-1),"Source"])
  y[repeats, "Source"] <- "  "
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
  print.data.frame(y, na.print="  ", right=FALSE, row.names=FALSE)
  if (!attr(x, which="orthogonal"))
    cat("\nThe design is not orthogonal\n\n")
  invisible(x)
}

"proj2.sweep" <- function(Q1, Q2, Eff.Q1.Q2)
{ #A procedure to compute the Residual operator for P remove Q when P and Q are nonorthogonal.
  #  Corresponding projection operator for Q in P is also obtained.
  #This version, which allows the supply of efficiency factors, is not being used
  n <- nrow(Q1)
  if (n != nrow(Q2))
    stop("Matrices not conformable.")
  isproj <- is.projector(Q1) & is.projector(Q2)

  if (length(Eff.Q1.Q2) == 1 & Eff.Q1.Q2[1]==0) #check matrices are orthogonal
  { Qconf <- projector(matrix(0, nrow = n, ncol = n))
    Qres <- Q1
    Eff.Q1.Q2 <- 0
    warning("Matrices are orthogonal.")
  }
  else
  { daeTolerance <- get("daeTolerance", envir=daeEnv)
    EffUnique.Q1.Q2 <-remove.repeats(Eff.Q1.Q2, daeTolerance[["eigen.tol"]])
    K <- length(EffUnique.Q1.Q2)
    #check for all confounded (i.e. eff = 1)
    if (K == 1 & EffUnique.Q1.Q2[1] == 1 & length(Eff.Q1.Q2) == degfree(Q2))
    { Qconf <- projector(Q2)
      Qres <- projector(Q1 - Q2)
    }
    else      #compute projection operators for partially confounded case
    { I <- diag(1, nrow = n, ncol = n)
      Qres <- I
      Q121 <- Q1 %*% Q2 %*% Q1
      for(i in 1:K)
        Qres <- Qres %*% (Q1 - (Q121/EffUnique.Q1.Q2[i]))
      Qres <- projector(Qres)
      Qconf <- projector(Q1 - Qres)
    }
  }
  list(Qconf = Qconf, Qres = Qres)
}


"projs.2canon" <- function(Q1, Q2)
  #Function to do an eigenanalysis of the relationship between two sets of projection matrices
{ if (!is.list(Q1) | !is.list(Q2))
    stop("Both Q1 and Q2 must be lists")
  daeTolerance <- get("daeTolerance", envir=daeEnv)
  
  #Get sizes and set up labels
  nQ1 <- length(Q1)
  nQ2 <- length(Q2)
  n <- nrow(Q1[1])
  Q1labels <- names(Q1)
  if (is.null(Q1labels))
    Q1labels <- as.character(1:nQ1)
  Q2labels <- names(Q2)
  if (is.null(Q2labels))
    Q2labels <- as.character(1:nQ2)
  criteria <- c('aefficiency','mefficiency','sefficiency','eefficiency','order',"dforthog")
  
  #Perform the analysis
  kQ1Q2 <- 0
  multieffic <- FALSE
  results <- vector(mode = "list", length = 0)
  for (i in Q1labels)
  { results[[i]][["Q1res"]] <- Q1[[i]]
    rdf <- degfree(Q1[[i]])
    for (j in Q2labels)
    { if (rdf >0) #only do this if there are df left in Q1
      { #Get unadjusted criteria and store only if confounded
        Q1Q2.eff <- suppressWarnings(proj2.efficiency(Q1[[i]], Q2[[j]]))
        if (Q1Q2.eff[1] > 0)
        { 
          #Get adjusted efficiencies and matrices
          adj.Q1Q2 <- suppressWarnings(proj2.combine(results[[i]][["Q1res"]], Q2[[j]]))
          if (degfree(adj.Q1Q2$Qconf) > 0)
          { 
            #Check for adjusted orthogonality
            Qfitlab <- names(results[[i]])[-1]
            if (length(Qfitlab) > 0)
            { for (k in Qfitlab)
              { Qjik <- Q2[[k]] %*% Q1[[i]] %*% Q2[[j]]
                if(!is.allzero(Qjik))
                  warning(paste(j,"and",k,"are partially aliased in",i, sep=" "))
              }
            }  
            
            results[[i]][[j]] <- vector(mode = "list", length = 0)
            #store pairwise efficiencies
            results[[i]][[j]][["pairwise"]][["efficiencies"]] <- Q1Q2.eff
            results[[i]][[j]][["pairwise"]][criteria] <- efficiency.criteria(Q1Q2.eff)
            
            if (degfree(adj.Q1Q2$Qres) == 0)
            { adj.Q1Q2$Qres <- matrix(0, nrow=nrow(adj.Q1Q2$Qres), ncol=ncol(adj.Q1Q2$Qres))
              adj.Q1Q2$Qres <- projector(adj.Q1Q2$Qres)
            }
            results[[i]][["Q1res"]] <- adj.Q1Q2$Qres
          
            #store adjusted efficiencies
            results[[i]][[j]][["adjusted"]][["efficiencies"]] <- adj.Q1Q2$efficiencies
            results[[i]][[j]][["adjusted"]][criteria] <- efficiency.criteria(adj.Q1Q2$efficiencies)
            if (results[[i]][[j]][["adjusted"]][["aefficiency"]] - 
                  results[[i]][[j]][["adjusted"]][["eefficiency"]] > daeTolerance[["eigen.tol"]])
              multieffic <- TRUE
            
            #Store adjusted projector
            results[[i]][[j]][["Qproj"]] <- adj.Q1Q2$Qconf
            rdf <- degfree(adj.Q1Q2$Qres)
          }
          else
            warning(paste(j,"is aliased with previous terms in",i, sep=" "))
        }  
      }
    }
  }
  class(results) <- "p2canon"
  return(results)
}

"summary.p2canon" <- function(object, which.criteria = c("aefficiency", "eefficiency", "order"), ...)
  #Routine to output a summary of the projector analysis
{ if (!inherits(object, "p2canon"))
    stop("object must be of class p2canon as produced by projs.2canon")
  #check which.criteria arguments
  criteria <- c("aefficiency", "eefficiency", "mefficiency", "sefficiency", "order", "dforthog")
  options <- c(criteria, "none", "all")
  kcriteria <- options[unlist(lapply(which.criteria, check.arg.values, 
                                     options=options))]
  if ("all" %in% kcriteria)
    kcriteria <- criteria
  anycriteria <- !("none" %in% kcriteria)
  
  #Form data frame for summary table
  orthogonaldesign <- TRUE
  nc <- 3
  if (anycriteria)
  { nc <- nc + length(kcriteria)
    res.criteria <- vector(mode = "list", length = length(kcriteria))
    names(res.criteria) <- kcriteria
    res.criteria[kcriteria] <- NA
    summary <- data.frame(matrix(nrow = 0, ncol=nc))
    colnames(summary) <- c("Source", "Confounded.source", "df", kcriteria)
  }
  else
  {   summary <- data.frame(matrix(nrow = 0, ncol=nc))
      colnames(summary) <- c("Source", "Confounded.source", "df")
  }  
  Q1labels <- names(object)
  for (i in Q1labels)
  { Q2labels <- names(object[[i]])[-1]
    nconf.terms <- 0
    if (length(Q2labels) > 0)
    { for (j in Q2labels)
      { kdf <- degfree(object[[i]][[j]]$Qproj)
        if (kdf > 0)
        { nconf.terms <- nconf.terms + 1
          if (abs(1 - unlist(object[[i]][[j]][["adjusted"]]["aefficiency"])) > 1e-04)
            orthogonaldesign <- FALSE
          if (anycriteria)
            summary <- rbind(summary,
                             data.frame(Source = i,
                                        Confounded.source = j, 
                                        df = kdf, 
                                        object[[i]][[j]][["adjusted"]][kcriteria], 
                                        stringsAsFactors = FALSE))
          else
            summary <- rbind(summary,
                             data.frame(Source = i,
                                        Confounded.source = j, 
                                        df = kdf, 
                                        stringsAsFactors = FALSE))
        }
      }
    }
    kdf <- degfree(object[[i]]$Q1res)
    if (kdf > 0)
      if (nconf.terms > 0)
      { if (anycriteria)
        summary <- rbind(summary,
                         data.frame(Source = i,
                                    Confounded.source = "Residual", 
                                    df = kdf, 
                                    res.criteria[kcriteria], 
                                    stringsAsFactors = FALSE))
        else
          summary <- rbind(summary,
                           data.frame(Source = i,
                                      Confounded.source = "Residual", 
                                      df = kdf, 
                                      stringsAsFactors = FALSE))
      }  
      else
      { if (anycriteria)
        summary <- rbind(summary,
                         data.frame(Source = i,
                                    Confounded.source = "   ", 
                                    df = kdf, 
                                    res.criteria[kcriteria], 
                                    stringsAsFactors = FALSE))
        else
          summary <- rbind(summary,
                           data.frame(Source = i,
                                      Confounded.source = "   ", 
                                      df = kdf, 
                                      stringsAsFactors = FALSE))
      }
  }
  summary <- as.data.frame(summary)
  class(summary) <- c("summary.p2canon", "data.frame")
  attr(summary, which = "title") <- 
    "\n\nSummary table of the decomposition (based on adjusted quantities)\n\n"
  attr(summary, which = "orthogonal") <- orthogonaldesign
  return(summary)
}


"efficiencies.p2canon" <- function(object, which = "adjusted")
#function to extract the efficiency factors from a p2canon object 
{ if (!inherits(object, "p2canon"))
    stop("object must be of class p2canon as produced by projs.2canon")
  options <- c("adjusted", "pairwise")
  opt <- options[check.arg.values(which, options)]
  
  #Get efficiencies
  Q1labels <- names(object)
  efficiencies <- vector(mode = "list", length = 0)
  for (i in Q1labels)
  { Q2labels <- names(object[[i]])[-1]
    if (length(Q2labels) > 0)
    { efficiencies[[i]] <- vector(mode = "list", length = 0)
      for (j in Q2labels)
        efficiencies[[i]][[j]] <- object[[i]][[j]][[opt]][["efficiencies"]]
    }  
  }
  return(efficiencies)
}  

"projs.combine.p2canon" <- function(object)
#function to extract, from a p2canon object, the projectors that give the combined decomposition  
{ 
  #Initialize list
  Q1combineQ2 <- vector(mode = "list", length = 0)
  
  #Loop through p2canon object
  Q1labels <- names(object)
  efficiencies <- vector(mode = "list", length = 0)
  for (i in Q1labels)
  { Q2labels <- names(object[[i]])[-1]
    if (length(Q2labels) > 0)
    { for (j in Q2labels)
      { Q1Q2label <- paste(Q1labels[[match(i, Q1labels)]], Q2labels[[match(j, Q2labels)]], sep="&")
        Q1combineQ2[[Q1Q2label]] <-object[[i]][[j]][["Qproj"]]
      }
      #Get the residual if any
      if (degfree(object[[i]][["Q1res"]]) > 0)
      { Q1Q2label <- paste(Q1labels[[match(i, Q1labels)]],"Residual", sep="&")
        Q1combineQ2[[Q1Q2label]] <- object[[i]][["Q1res"]]
      }
    }
    else
      Q1combineQ2[[i]] <- object[[i]][["Q1res"]]
  }
  return(Q1combineQ2)
}  
