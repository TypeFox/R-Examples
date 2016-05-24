`preparePredictor.fnc` <-
function(pred, model, m, ylabel, fun, val, xlabel, ranefs, ...) {
   


  ###############################################################################
  # we first figure out whether we are dealing with a polynomial predictor, or
  # perhaps a restricted cubic spline
  ###############################################################################

  X = model@pp$X

  polynomial = FALSE
  namesplit = strsplit(pred, ", ")[[1]]
  a = regexpr("poly\\(", namesplit[1])
  if ((a==1) & (attr(a, "match.length")==5)) {
    polynomial = TRUE
    degree = degreesOrKnots.fnc(pred)
  }

  rcspline = FALSE
  namesplit = strsplit(pred, ", ")[[1]]
  a = regexpr("rcs\\(", namesplit[1])
  if ((a==1) & (attr(a, "match.length")==4)) {
    rcspline = TRUE
    knots = degreesOrKnots.fnc(pred)
  }

  if ((!polynomial) & (!rcspline)) {
    pred2 = paste("rcs\\(", pred, sep="")
    if (length(grep(pred2, colnames(X))) > 0) {
      rcspline = TRUE
    } else {
      pred2 = paste("poly\\(", pred, sep="")
      if (length(grep(pred2, colnames(X))) > 0) {
        polynomial = TRUE
      }
    }
  }

  isfactor = FALSE
  ###############################################################################
  #  we first handle the case that a predictor is not polynomial nor spline (if)
  #  and later (else) handle the polynomial case
  ###############################################################################

  fixefs = fixef(model) 
  if (!is.na(ranefs[[1]])) {
     nm = as.vector(ranefs[[4]])
     if (nm %in% names(fixefs)) {
       blup = ranef(model)[[ranefs[[1]]]][ranefs[[2]],ranefs[[3]]]
       fixefs[nm] = fixefs[nm]+blup
       fixefs = as.numeric(fixefs)
     } else 
       stop(paste(nm, "is not a valid predictor name, check 'fixef(model)'\n", sep=" "))
  }

  if ((pred %in% colnames(model@frame)) & polynomial==FALSE & rcspline==FALSE) {
    if (is.numeric(model@frame[,pred])) {
      # oke, so this is a numeric predictor
      if (pred %in% colnames(X)) {
        m[,pred] = seq(min(X[,pred]), max(X[,pred]), length=nrow(m))
        # adjust m so that interactions are now properly represented in the columns of m
        m = implementInteractions.fnc(m)
        # calculate predicted values
        vals = m %*% as.numeric(fixefs)
        # if necessary apply transformation to expected values
        vals = transforming.fnc(vals, fun)
        dfr = data.frame(X=m[,pred], Y = vals)
        dfr$Predictor = rep(xlabel, nrow(dfr))
        dfr$Type = rep(isfactor, nrow(dfr))
        if (is.na(val)) {
          dfr$Interaction = rep(NA, nrow(dfr))
        } else {
          dfr$Interaction = rep(val, nrow(dfr))
        }
      } else {
        stop(paste(pred, " is not plotted (not a fixed effect predictor)\n"))
      }
    } else {
      if (is.logical(model@frame[,pred])) model@frame[,pred] = factor(model@frame[,pred])
      if (is.factor(model@frame[,pred])) {
        #-------------------------------------------------------------------
        # so now we are dealing with a factor as predictor, and we need
        # to reconstruct the factor names
        #-------------------------------------------------------------------
        isfactor=TRUE
        factnames = paste(pred, levels(model@frame[,pred])[-1],sep="")
        m = m[1:(length(factnames)+1),]
        for (i in 1:length(factnames)) {
          m[i+1, factnames[i]] = 1
        }
        # and then implement the proper interactions in the model matrix
        m = implementInteractions.fnc(m)
        # calculate the expected values
        vals = m %*% as.numeric(fixefs)
        # and transform expected values if so desired
        vals = transforming.fnc(vals,fun)
        x = 1:nrow(m)
        dfr = data.frame(X=x, Y = vals)
        dfr$Predictor = rep(xlabel, nrow(dfr))
        dfr$Type = rep(isfactor, nrow(dfr))
        if (is.na(val)) {
          dfr$Interaction = rep(FALSE, nrow(dfr))
        } else {
          dfr$Interaction = rep(TRUE, nrow(dfr))
        }
        dfr$Levels = levels(model@frame[,pred])
      } else {
        cat("warning: I don't know how to handle ", pred, "\n")
      }
    }
  }  # of if pred in colnames
  else {

    ##########################################################################
    # some preprocessing
    ##########################################################################
 
    if (!(pred %in% colnames(X))) {
      # this is the case in which the predictor is not in the list of predictors
      # and has not been detected as polynomial or spline; in other words, the user is
      # specifying the predictor by its name, e.g., "PrevRT" while in the model
      # formula it is specified as "poly(PrevRT, degree)" or "rcs(prevRT, knots)".  
      # So we reconstruct the name for the polynomial terms or for the spline terms, 
      # and also extract the degree of the polynomial c.q. the number of knots
      pos = grep(pred, colnames(X), fixed=TRUE) 
      degree = 1
      knots = 1
      if (length(pos) > 0) {
        name=colnames(X)[pos][1]
        namesplit = strsplit(name, ", ")[[1]]
        # check whether polynomial
        a = regexpr("poly", namesplit[1])
        if ((a==1) & (attr(a, "match.length")==4)) {
          polynomial = TRUE
          degree = as.numeric(namesplit[2])
          xlabel = parsePredName.fnc(pred)[[1]]
          name = pred
        }
        if (!polynomial) {
          # check whether restricted cubic spline
          a = regexpr("rcs", namesplit[1])
          if ((a==1) & (attr(a, "match.length")==3)) {
            rcspline = TRUE
            #knots = as.numeric(substr(namesplit[2], 1, nchar(namesplit[2])-1))
            aa = parsePredName.fnc(name)
            knots = aa[[2]]
            xlabel = aa[[1]]
          }
        }
      } 
    }  # of if pred in colnames
    else {  # so we know this is a term in the model frame, 
      # so we are dealing with something like "poly(PrevRT, 2, raw = TRUE)"
      namesplit = strsplit(pred, ", ")[[1]]
      name = pred
      arg2 = as.numeric(substr(namesplit[2], 1, nchar(namesplit[2])-1))
      cat("DIT ZOU DOOD STUK CODE MOETEN ZIJN\n")
    }
         
    if (polynomial|rcspline) {
      if (is.na(xlabel)) {
        xlabel = pred
      }
      if (polynomial) {   # we install the appropriate values in the columns in m
        hasPoly = FALSE
        if (length(grep("^poly\\(", pred))>0) {
          #xlabel = strsplit(strsplit(name, " ")[[1]][1],"[^a-zA-Z]")[[1]][2]
          #degree = as.numeric(strsplit(name, ",")[[1]][2])
          vec = paste(name, "1", sep="")
          hasPoly = TRUE
        } else {
          xlabel = pred
          vec = paste("poly(", name, ", ", degree, ", raw = TRUE)1", sep="")
        }
        name1 = vec  # name of first column, need this for MCMC intervals
        m[,vec] = seq(min(X[,vec]), max(X[,vec]), length=nrow(m))
        for (i in 2:degree) {
          if (hasPoly) {
            vec = c(vec, paste(name, as.character(i), sep=""))
          } else {
            vec = c(vec, paste("poly(", name, ", ", degree, ", raw = TRUE)", as.character(i), sep=""))
          }
          m[,vec[i]] = m[,vec[i-1]]*m[,vec[1]]
        }
      } else {  # restricted cubic spline
        if (length(grep("^rcs\\(", pred))>0) {
          nms = unlist(parsePredName.fnc(pred))
          basename = nms[1]
          knots = as.numeric(nms[2])
          name1 = paste("rcs(", basename, ", ", knots, ")", basename, sep="")
          xlabel = basename
        } else {
          knots = getKnots.fnc(colnames(X), pred)
          name1 = paste("rcs(", pred, ", ", knots, ")", pred, sep="")
        }
        vec = rep(name1, knots-1)
        vec[2] = paste(vec[1], "'", sep="")
        if (knots > 3) {
          for (i in 3:(knots-1)) {
            vec[i] = paste(vec[i-1], "'", sep="")
          }
        }
        mtmp = unique(X[,vec])
        if (nrow(mtmp) <= nrow(m)) {
          m = m[1:nrow(mtmp),]
          m[,vec] = mtmp
        } else {
          vecIndices = c(1, sort(sample(2:(nrow(mtmp)-1), nrow(m)-2)), nrow(mtmp))
          m[,vec] = mtmp[vecIndices,]
        }
        m = m[order(m[, vec[1]]),]
      }
      # make sure the consequences for interactions are taken into account
      m = implementInteractions.fnc(m)
      # calculate expected values
      vals = m %*% as.numeric(fixefs)
      # and transform the expected values if necessary
      vals = transforming.fnc(vals, fun)
      dfr = data.frame(X=m[,vec[1]], Y = vals)
      dfr$Predictor = rep(xlabel, nrow(dfr))
      dfr$Type = rep(isfactor, nrow(dfr))
      if (is.na(val)) {
        dfr$Interaction = rep(FALSE, nrow(dfr))
      } else {
        dfr$Interaction = rep(val, nrow(dfr))
      }

    } # if spline or polynomial
    else {
      stop(paste("unknown function used in", pred, "\n"))
    }
  }
  return(dfr)
}

