# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.
xp.bootgam <- function (object,
                        n = NULL,                   # number of bootstrap iterations
                        id = "ID",                 # column name of id
                        oid = "OID",               # create a new column with the original ID data
                        seed = NULL,               # random seed
                        parnam = xvardef("parms", object),
                        covnams = xvardef("covariates", object),
                        wts.col = NULL,
                        ask.for.input = TRUE,
                        overwrite = TRUE,
                        ...) {
  ask.for.n <- function(...) {
    cat("\nEnter maximum number of bootstrap replicates to perform in bootGAM (0 to exit): ")
    ans <- as.numeric(readline())
    if ((!is.numeric(ans)) || (is.na(ans))) {
      cat("\nPlease specify a number.\n")
      ans <- Recall(...)
    }
    return(as.numeric(ans))
  } 
  ask.for.seed <- function(...) {
    cat("\nSpecify a seed number for the bootstrap (0 to exit): ")
    ans <- as.numeric(readline())
    if ((!is.numeric(ans)) || (is.na(ans))) {
      cat("\nPlease specify a number.\n")
      ans <- Recall(...)
    }
    if (ans == 0) {
      return(NULL)
    }
    return(as.numeric(ans))
  } 
  ask.for.par <- function(...) {
    cat("\nEnter name of parameter to use for this bootGAM (0 to exit): ")
    ans <- readline()
    if (ans == 0) {
      return(NULL)
    }
    if (length(ans) > 1) {
      cat("\nYou have specified more than one parameter.\n")
      cat("The bootGAM can be run on only one parameter at a time.\n")
      ans <- Recall(...)
    }
    else {
      ans.exists <- check.vars(ans, object)
      if (is.null(ans.exists)) {
        ans <- Recall(...)
      }
    }
    return(ans)
  }
  get.par <- function(nams, get.input = FALSE, ...) {
    ans <- NULL
    if (length(nams) == 0) {
      cat("\nNo parameter is defined for this bootGAM\n")
      if (get.input) {
                ans <- ask.for.par()
              }
      else {
        cat("\nType '?xp.bootgam' for more information.\n")
      }
    }
    if (length(nams) > 1) {
      cat("\nThere is more than one parameter defined\n")
      cat("for this bootGAM. The parameters defined are:\n\n")
      cat(nams, fill = 60)
      cat("\nThe bootGAM can be run on only one parameter at a time.\n")
      if (get.input) {
        ans <- ask.for.par()
      }
      else {
        cat("\nType '?xp.bootgam' for more information.\n")
      }
    }
    if (length(nams) == 1) {
      ans <- nams
    }
    return(ans)
  }
  ask.for.wts <- function(...) {
    cat("\nWeight column to use (0 to exit, NULL for no weight):")
    ans <- readline()
    if (ans == "NULL") 
      return("NULL")
    if (ans == 0) 
      return(NULL)
    if (length(ans) > 1) {
      cat("\nYou have specified more than one weight.\n")
      cat("Only one weight is allowed.\n")
      ans <- Recall(...)
    }
    else {
      if (is.na(pmatch(ans, names(object@Data.firstonly)))) {
        cat(paste("\n", "-----------Variable(s) not defined!-------------\n", 
                  ans, "is not defined in the current database\n", 
                  "and must be defined for this command to work!\n", 
                  "------------------------------------------------\n"))
        ans <- Recall(...)
      }
      return(ans)
    }
  }
  get.wts <- function(nams, get.input = FALSE, ...) {
    ans <- NULL
    if (length(nams) == 0) {
      cat("\nNo weights are defined for this bootGAM\n")
      if (get.input) {
        ans <- ask.for.wts()
      }
      else {
        cat("\nType '?xp.bootgam' and '?xpose.bootgam' for more information.\n")
      }
    }
    if (length(nams) > 1) {
      cat("\nPlease specify a the weights for the parameter.\n")
      cat("The weights come from columns in the data contained\n")
      cat("in the Data.firstonly section of the xpose data object.\n")
      cat("These values usualy come from the .phi file of a NONMEM run.\n")
      cat("Possible weight values (column names) are:\n\n")
      cat(nams, fill = 60)
      cat("\nOnly one weight can be specified.\n")
      if (get.input) {
        ans <- ask.for.wts()
      }
      else {
        cat("\nType '?xp.bootgam' and '?xpose.bootgam' for more information.\n")
      }
    }
    if (length(nams) == 1) {
      ans <- nams
    }
    return(ans)
  }
  pars <- get.par(parnam, get.input = ask.for.input, ...)
  if (is.null(pars)) {
    return(invisible())
  }
  if (object@Prefs@Bootgam.prefs$n) {
    n <- object@Prefs@Bootgam.prefs$n
  } else {
    if (ask.for.input) {
      n <- ask.for.n()
      if (is.null(n)) {
        return(invisible())
      }
    }
  }
  if (object@Prefs@Gam.prefs$wts & ask.for.input) {
    wts <- get.wts(names(object@Data.firstonly), get.input = ask.for.input, 
                   ...)
    if (is.null(wts)) {
      return(invisible())
    }
    if (wts == "NULL") {
      wts <- NULL
    }
    wts.col <- wts
  }
  if (exists(paste("bootgam.xpose.", pars, ".", object@Runno, sep = ""), 
             where = 1) & !overwrite) {
    if (ask.for.input) {
      cat("\nThere is already a bootgam object associated with the current\n")
      cat("run number and parameter. It will be overwritten if you proceed.\n")
      cat("Proceed? n(y): ")
      ans <- readline()
      if (ans != "y") {
        return()
      }
    }
    else {
      cat("\nThere is already a bootgam object associated with the current\n")
      cat("run number and parameter. It will NOT be overwritten.\n")
      return()
    }
  }

  # Invoke actual bootGAM
  if (is.null(seed)) {
    seed <- ask.for.seed()
    if (is.null(seed)) {
      return(invisible)
    }
  }

  bootgam.obj <- xpose.bootgam (object,                        # Xpose database
                                n = n,                         # number of replicates
                                id = object@Prefs@Xvardef$id,  # column label of id
                                oid = "OID",                   # create a new column with the original ID data
                                seed = seed,                   # random seed
                                parnam = pars,
                                covnams = covnams,
                                wts.col = wts.col,
                                ...)


  cat ("\nBootstrap completed\n")
  if (bootgam.obj$algo == "fluct.ratio") {
    cat (paste("Fluctuation ratio: ", bootgam.obj$fluct.ratio.last, "\n\n", sep = ""))
  } else {
    cat (paste("Lowest absolute joint inclusion frequency: ", bootgam.obj$ljif.last, "\n\n", sep = ""))  
  }
  if (length(bootgam.obj$failed[bootgam.obj$failed == 1]) > 0) {
    cat ("Note: the gam failed for the following bootstrap replicate(s), since not all\ncovariates had enough levels (>1) in the dataset: \n")
    cat (paste ( seq(1:n)[bootgam.obj$failed == 1] ) )
    cat ("\nThese bootstrap replicates will not be used in summaries or diagnostic plots.\n")
    cat ("\n\n")
  }
  c1 <- call("assign",pos = 1, paste("bootgam.xpose.", pars, ".", object@Runno, 
           sep = ""), bootgam.obj, immediate = T)
  eval(c1)
  if (exists("current.bootgam", where = 1)) {
    remove(pos = 1, "current.bootgam")
  }
  c2 <- call("assign",pos = 1, "current.bootgam", bootgam.obj, immediate = T)
  eval(c2)
  return(invisible(bootgam.obj))
  
}

data.long <- function (data) { # same as melt() in reshape: create a dataframe in the long format based on a matrix/data.frame (for xyplot)
  all <- c()
  for (i in seq(along=colnames(data))) {
    all <- data.frame (rbind (all, cbind(var = colnames(data)[i], row = 1:length(data[,1]), value = data[[colnames(data)[i]]] )))
  }
  all[,2] <- as.numeric(as.character(all[,2]))
  all[,3] <- as.numeric(as.character(all[,3]))
  return (all)
}

xpose.bootgam <- function (object,
                           n = n,                         # number of bootstrap iterations
                           id = object@Prefs@Xvardef$id,  # column name of id
                           oid = "OID",                   # create a new column with the original ID data
                           seed = NULL,                   # random seed
                           parnam = xvardef("parms", object)[1],
                           covnams = xvardef("covariates", object),
                           conv.value = object@Prefs@Bootgam.prefs$conv.value, 
                           check.interval = as.numeric(object@Prefs@Bootgam.prefs$check.interval),   
                           start.check = as.numeric(object@Prefs@Bootgam.prefs$start.check),   
                           algo = object@Prefs@Bootgam.prefs$algo,   
                           start.mod = object@Prefs@Bootgam.prefs$start.mod,   
                           liif = as.numeric(object@Prefs@Bootgam.prefs$liif),   
                           ljif.conv = as.numeric(object@Prefs@Bootgam.prefs$ljif.conv),
                           excluded.ids = as.numeric(object@Prefs@Bootgam.prefs$excluded.ids),
                           ...) {
  ids   <- unique (object@Data[[id]])
  if (!is.null(seed)) {
    set.seed (seed)
  }

  ## create template object 
  bootgam.obj <- list("results.raw" = list(),
                      "results.tab" = data.frame (matrix(0, nrow = n, ncol=length(covnams),
                        dimnames = list(NULL, covnams))),
                      "incl.freq" = data.frame (matrix(0, nrow = n, ncol=length(covnams),
                        dimnames = list(NULL, covnams))),
                      "oid"       = data.frame (matrix(0, nrow = n, ncol=length(ids))),
                      "seed"      = seed,
                      "runno"     = object@Runno,
                      "failed"    = rep(0, n),
                      "parnam"    = parnam,
                      "covnams"   = covnams,
                      "n" = n,
                      "conv.value" = conv.value,
                      "check.interval" = check.interval,
                      "start.check"= start.check,
                      "algo"      = algo,
                      "start.mod" = start.mod,
                      "liif"      = liif, 
                      "ljif.conv" = ljif.conv,
                      "excluded.ids" = excluded.ids,
                      "fluct.ratio.last" = NULL,
                      "ljif.last" = NULL
                      )

  ## Define convergence criteria
  get.crit1 <- function(cov.table, check.interval, i, conv) {
    #cov.table <- current.bootgam$incl.freq
    n.del <- sum(apply (cov.table[1:i,], 1, sum) == 0)
    i <- i - n.del
    cov.table <- cov.table[apply (cov.table, 1, sum) > 0,] # unselect failed gams
    i.start <- i - check.interval
    if (i.start < 1) {i.start <- 1}
    min.int <- apply (cov.table[i.start:i,], 2, min)
    max.int <- apply (cov.table[i.start:i,], 2, max)
    prob <- apply(cov.table[i.start:i, ], 2, mean)
    del <- min.int != 0  # if min.int is 0, then division by zero! remove those
    min.int <- min.int[del]
    max.int <- max.int[del]
    prob <- prob[del]
    crit <- sum((max.int/min.int) * prob)/sum(prob)
    cnv <- 0
    if(!is.na(crit)) { 
      if(is.numeric(crit) && (crit < conv)) {
        cnv <- 1
      }
    } 
    retrn <- c(crit, cnv)
    return(retrn)
  }
  get.crit2 <- function(ind.mat, i, liif = 0.2, ljif.conv = 25, failed) {
    ind.mat <- ind.mat[failed == 0,]
    n.del <- sum(failed[1:i])
    i <- i - n.del
    min.inc <- min(apply(ind.mat[1:i,], 2, mean))
    f.low   <- i * min.inc * liif
    cnv <- 0
    if(f.low >= ljif.conv) {
      cnv <- 1
    }
    retrn <- c (f.low,cnv)
    return(retrn)
  }
  
  ## Initialize bootstrap
  cat (paste("\nStarting bootstrap (max ", n, " replicates, seed = ", seed,").\n", sep=""))
  cat ("Bootgam progress:\n")
  pb <- txtProgressBar(min=0, max=n, initial=1, style=3)
  next.check <- as.numeric(start.check)
  i <- 1
  
  ## filter out excluded individuals (if exist in dataset)
  for (i in seq(along=excluded.ids)) {
    if (excluded.ids[i] %in% object@Data[[id]]) {
      object@Data <- object@Data[object@Data[[id]]!=excluded.ids[i],]
    }
  }
  
  ## The actual bootstrap loop
  i <- 1
  while (i <= n) { 
    setTxtProgressBar(pb, i)
    bs.object <- object
    bs.sample <- xpose.bootgam.drawsample (object@Data)
    if (test.covariate.levels (bs.sample$Data, covnams) == 0) { # test number of levels
      bs.object@Data <- bs.sample$Data           # Bootstrapped dataset
      for (k in seq(along = bs.sample$id)) {
        bootgam.obj$oid[i,bs.sample$id[k]] <- bootgam.obj$oid[i,bs.sample$id[k]] + 1                   # save original ID numbers selected in the bootstrap
      }
      capture.output (                           # don't show all print output from the gam
                      gam.results <- xpose.gam (bs.object,
                                                start.mod = start.mod,
                                                parnam = bootgam.obj$parnam)
                      )
      bootgam.obj$results.raw[[i]] <- gam.results
      covs <- xpose.bootgam.extract.covnames(gam.results)
      bootgam.obj$results.tab[i, match (covs, covnams)] <- 1
      tmp <- bootgam.obj$results.tab[bootgam.obj$failed==0,]
      bootgam.obj$incl.freq[i,] <- apply (tmp[1:(i-sum(bootgam.obj$failed)),], 2, mean)
    } else {
      bootgam.obj$results.tab[i, ] <- -99 
      bootgam.obj$failed[i] <- 1
    }

    ## Convergence criteria
    if(i == next.check) { # If it is time to check the convergence
      next.check <- as.numeric(next.check) + as.numeric(check.interval)
      ## The criteria depends on the algo
      if(algo == "fluct.ratio") {
	crit.vec <- get.crit1(bootgam.obj$incl.freq, check.interval, i, conv.value)
        bootgam.obj$fluct.ratio.last <- crit.vec[1]
      } else {
	crit.vec <- get.crit2(bootgam.obj$oid, i, liif, ljif.conv, bootgam.obj$failed)
        bootgam.obj$ljif.last <- crit.vec[1]
      }
      cat (paste(" Conv. crit :", round(crit.vec[1],3), "\n"))	
      ## Convergence!!
      if(crit.vec[2] == 1) {
	  bootgam.obj$results.tab <- bootgam.obj$results.tab[1:i, ]
	  bootgam.obj$incl.freq <- bootgam.obj$incl.freq[1:i, ]
          bootgam.obj$failed <- bootgam.obj$failed[1:i]
          n <- 0 # stop the bootstrap
	}
    }
    i <- i+1
  }
  close(pb)
  return (bootgam.obj)
}

test.covariate.levels <- function (data, covariates) {
  # all covariates must at least have more than 1 level to be able to run a gam, otherwise it will crash
  flag <- 0
  for (i in 1:length(covariates)) {
    if (length(unique(data[[covariates[i]]])) < 2) {
      flag <- 1
    }
  }
  return (flag)
}

xpose.bootgam.extract.covnames <- function (gam.object) {
  covs <- names(gam.object[1]$coefficients)[-1]  # first element is intercept    
  covs <- gsub ("ns\\(", "", covs)
  covs <- gsub ("\\,.*", "", covs)
  covs <- gsub ("[0-9]*", "", covs)
  return (unique (covs))
}

xpose.bootgam.drawsample <- function (data,                      # data.frame() 
                                      id = "ID",                 # column name of id
                                      oid = "OID"               # create a new column with the original ID data
                                      ) {
  errors_encountered <- c()
  
  ## get vector with individuals in dataset
  ids   <- unique (data[[id]])
  n_ids <- length(ids)
  if (n_ids <= 1) {
    errors_encountered <- c("Please check arguments supplied to bootstrap function.")
  }

  if (is.null(errors_encountered)) {# only draw bootstrap sample when dataset seems valid
    
    ## initialize bootstrap step
    id_draw <- sample (n_ids, n_ids, replace = T)
    bs_draw <- ids[id_draw]
    if (length(oid) > 0) { # store original id number as new column
      data[[oid]] <- data[[id]] 
    }

    bs_data <- c()
    ## Loop over the ids that were drawn and give a new number
    for (i in 1:n_ids) {
      tmp_data <- data[data[[id]] == bs_draw[i],]
      tmp_data[[id]] <- i  # give a new id number
      bs_data <- rbind (bs_data, tmp_data)
    }
    bs.sample <- list ("Data" = bs_data, "oid" = bs_draw, "id" = id_draw)
    return (bs.sample)
  }
}
