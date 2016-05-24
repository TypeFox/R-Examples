generateValues <- function(dataSample, dataPop, params) {
  if ( !nrow(dataSample) ) {
    return(character())
  }

  meth <- params$method
  cur.var <- params$cur.var
  excludeLevels <- params$excludeLevels
  hasNewLevels <- params$hasNewLevels
  newLevels <- params$newLevels
  w <- params$w
  formula.cmd <- params$formula.cmd
  eps <- params$eps
  limit <- params$limit[[cur.var]]
  censor <- params$censor[[cur.var]]
  if(nrow(dataSample[!duplicated(dataSample[,cur.var,with=FALSE])])==1){
    invisible(unlist(head(dataSample[,cur.var,with=FALSE],1)))
  }else{
    # temporarily recode response vector
    dataSample[,cur.var:=cleanFactor(.SD),.SDcols=cur.var,with=FALSE]
    levelsResponse <- levels(unlist(dataSample[,cur.var,with=FALSE]))

    #indices for unique occurence
    indGrid <- split(1:nrow(dataPop), dataPop, drop=TRUE)

    #get only the first obs with unique combinations
    grid <- dataPop[sapply(indGrid, function(i) i[1])]

    # in sample, observations with NAs have been removed to fit the
    # model, hence population can have additional levels
    # these need to be removed since those probabilities cannot
    # be predicted from the model
    if ( excludeLevels ) {
      exclude <- mapply(function(pop, new) pop %in% new,
          pop=grid[, hasNewLevels, drop=FALSE,with=FALSE],
          new=newLevels[hasNewLevels])
      if ( is.null(dim(exclude)) ) {
        exclude <- which(any(exclude))
      } else {
        exclude <- which(apply(exclude, 1, any))
      }
    } else {
      exclude <- integer()
    }
    # fit multinomial model
    # command needs to be constructed as string
    # this is actually a pretty ugly way of fitting the model
    mod <- eval(parse(text=formula.cmd))  # fitted model

    # predict probabilities
    if ( length(exclude) == 0 ) {
      newdata <- copy(grid)
    } else {
      newdata <- copy(grid[-exclude])
    }
    ind <- match(colnames(newdata), colnames(dataSample))
    for ( i in 1:length(ind) ) {
      if (is.factor(unlist(newdata[,i,with=FALSE]))) {
        newdata[,colnames(newdata)[i]:=factor(as.character(unlist(newdata[,colnames(newdata)[i],with=FALSE])),levels(dataSample[[ind[i]]])),with=FALSE]
      }
    }

    if ( meth %in% "multinom" ) {
      probs <- predict(mod, newdata=newdata, type="probs")
    }else if ( meth %in% c("ctree","cforest") ) {
      probs <- predict(mod, newdata=data.table(newdata), type="prob")
      probs <- do.call("rbind",probs)
    }
    #if ( meth %in% "naivebayes" ) {
    #  probs <- predict(mod, newdata=newdata, type="raw")
    #}
    # TODO: fix error if level-sets are not equal!
    #if ( meth %in% "ctree" ) {
    #  probs <- do.call("rbind", predict(mod, newdata=newdata, type="prob"))
    #}
    # set too small probabilities to exactly 0
    if ( !is.null(eps) ) {
      probs[probs < eps] <- 0
    }

    # ensure code works for missing levels of response
    ind <- as.integer(which(table(dataSample[[cur.var]]) > 0))
    if( length(ind) > 2 && (nrow(grid)-length(exclude)) == 1 ) {
      probs <- t(probs)
    }

    # account for structural zeros
    if ( (!is.null(limit) || !is.null(censor)) && !is.null(dim(probs)) ) {
      if(length(exclude) == 0) {
        probs <- adjustProbs(probs, grid, names(indGrid), limit[[i]], censor[[i]])
      } else {
        probs <- adjustProbs(probs, grid[-exclude, , drop=FALSE], names(indGrid)[-exclude], limit[[i]], censor[[i]])
      }
    }
    # local function for sampling from probabilities
    if ( length(ind) == 1 ) {
      resample <- function(k, n, p) rep.int(1, n[k])
    } else if ( length(ind) == 2 ) {
      resample <- function(k, n, p) spSample(n[k], c(1-p[k],p[k]))
    } else {
      resample <- function(k, n, p) spSample(n[k], p[k,])
    }
    # generate realizations for each combination

    if ( length(exclude) == 0 ) {
      ncomb <- as.integer(sapply(indGrid, length))
      sim <- lapply(1:length(ncomb), resample, ncomb, probs)
    } else {
      ncomb <- as.integer(sapply(indGrid[-exclude], length))
      sim <- as.list(rep.int(NA, length(indGrid)))
      sim[-exclude] <- lapply(1:length(ncomb), resample, ncomb, probs)
    }
    sim <- unsplit(sim, dataPop, drop=TRUE)
    invisible(levelsResponse[ind][sim])
  }
}

# simulation of variables using random draws from the observed
# conditional distributions of their multivariate realizations
generateValues_distribution <- function(dataSample, dataPop, params) {
  grid <- params$grid
  additional <- params$additional
  basic <- params$basic
  w <- params$w

  if( !nrow(dataSample) ) {
    return(character())
  }

  # population data
  splitS <- split(1:nrow(dataSample), dataSample[, basic, with=F], drop=TRUE)
  pSplit <- lapply(splitS, function(i) {
        tmp <- tableWt(dataSample[i, additional, with=F], dataSample[[w]][i])
        tmp <- as.data.frame(tmp)
        p <- ncol(tmp)
        tmp[, p]/sum(tmp[, p])
      })
  splitP <- split(1:nrow(dataPop), dataPop[, basic, with=F])
  NSplit <- sapply(splitP, length)
  # in sample, observations with NAs have been removed to fit the
  # model, hence population can have additional levels
  whichP <- which(names(splitP) %in% names(splitS))
  # generate realizations for each combination
  sim <- as.list(rep.int(NA, length(splitP)))
  sim[whichP] <- mapply(spSample, NSplit[whichP], pSplit, SIMPLIFY=FALSE)
  sim <- unsplit(sim, dataPop[, basic, with=F])
  sim <- grid[sim,,drop=FALSE]
  rownames(sim) <- rownames(dataPop)
  sim
}






#' Simulate categorical variables of population data
#'
#' Simulate categorical variables of population data. The household structure
#' of the population data needs to be simulated beforehand.
#'
#' The number of cpus are selected automatically in the following manner. The
#' number of cpus is equal the number of strata. However, if the number of cpus
#' is less than the number of strata, the number of cpus - 1 is used by
#' default. This should be the best strategy, but the user can also overwrite
#' this decision.
#'
#' @name simCategorical
#' @param simPopObj a \code{simPopObj} containing population and household
#' survey data as well as optionally margins in standardized format.
#' @param additional a character vector specifying additional categorical
#' variables available in the sample object of \code{simPopObj} that should be
#' simulated for the population data.
#' @param method a character string specifying the method to be used for
#' simulating the additional categorical variables. Accepted values are
#' \code{"multinom"} (estimation of the conditional probabilities using
#' multinomial log-linear models and random draws from the resulting
#' distributions) or \code{"distribution"} (random draws from the observed
#' conditional distributions of their multivariate realizations).
#' \code{"ctree"}  for using Classification trees
#' \code{"cforest"}  for using random forest
#' @param limit if \code{method} is \code{"multinom"}, this can be used to
#' account for structural zeros. If only one additional variable is requested,
#' a named list of lists should be supplied. The names of the list components
#' specify the predictor variables for which to limit the possible outcomes of
#' the response. For each predictor, a list containing the possible outcomes of
#' the response for each category of the predictor can be supplied. The
#' probabilities of other outcomes conditional on combinations that contain the
#' specified categories of the supplied predictors are set to 0. If more than
#' one additional variable is requested, such a list of lists can be supplied
#' for each variable as a component of yet another list, with the component
#' names specifying the respective variables.
#' @param censor if \code{method} is \code{"multinom"}, this can be used to
#' account for structural zeros. If only one additional variable is requested,
#' a named list of lists or \code{data.frame}s should be supplied. The names of
#' the list components specify the categories that should be censored. For each
#' of these categories, a list or \code{data.frame} containing levels of the
#' predictor variables can be supplied. The probability of the specified
#' categories is set to 0 for the respective predictor levels. If more than one
#' additional variable is requested, such a list of lists or \code{data.frame}s
#' can be supplied for each variable as a component of yet another list, with
#' the component names specifying the respective variables.
#' @param maxit,MaxNWts control parameters to be passed to
#' \code{\link[nnet]{multinom}} and \code{\link[nnet]{nnet}}. See the help file
#' for \code{\link[nnet]{nnet}}.
#' @param eps a small positive numeric value, or \code{NULL} (the default). In
#' the former case and if \code{method} is \code{"multinom"}, estimated
#' probabilities smaller than this are assumed to result from structural zeros
#' and are set to exactly 0.
#' @param nr_cpus if specified, an integer number defining the number of cpus
#' that should be used for parallel processing.
#' @param regModel allows to specify the variables or model that is used when
#' simulating additional categorical variables. The following choices are
#' available if different from NULL.  \itemize{ \item'basic'only the basic
#' household variables (generated with \code{\link{simStructure}}) are used.
#' \item'available'all available variables (that are common in the sample and
#' the synthetic population such as previously generated varaibles) excluding
#' id-variables, strata variables and household sizes are used for the
#' modelling. This parameter should be used with care because all factors are
#' automatically used as factors internally.  \item formula-objectUsers may also
#' specify a specifiy formula (class 'formula') that will be used. Checks are
#' performed that all required variables are available.  } If method
#' 'distribution' is used, it is only possible to specify a vector of length
#' one containing one of the choices described above.  If parameter 'regModel'
#' is NULL, only basic household variables are used in any case.
#' @param seed optional; an integer value to be used as the seed of the random
#' number generator, or an integer vector containing the state of the random
#' number generator to be restored.
#' @param verbose set to TRUE if additional print output should be shown.
#' @param by defining which variable to use as split up variable of the estimation. Defaults to the strata variable.
#' @return An object of class \code{\linkS4class{simPopObj}} containing survey
#' data as well as the simulated population data including the categorical
#' variables specified by argument \code{additional}.
#' @note The basic household structure needs to be simulated beforehand with
#' the function \code{\link{simStructure}}.
#' @author Bernhard Meindl and Andreas Alfons and Stefan Kraft
#' @seealso \code{\link{simStructure}}, \code{\link{simRelation}},
#' \code{\link{simContinuous}}, \code{\link{simComponents}}
#' @export
#' @keywords datagen
#' @examples
#' data(eusilcS) # load sample data
#' \dontrun{
#' ## approx. 20 seconds computation time
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
#' ## in the following, nr_cpus are selected automatically
#' simPop <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
#' simPop <- simCategorical(simPop, additional=c("pl030", "pb220a"), method="multinom", nr_cpus=1)
#' simPop
#' }
simCategorical <- function(simPopObj, additional,
    method=c("multinom", "distribution","ctree","cforest"),
    limit=NULL, censor=NULL, maxit=500, MaxNWts=1500,
    eps=NULL, nr_cpus=NULL, regModel=NULL, seed=1,
    verbose=FALSE,by="strata") {

  x <- newAdditionalVarible <- NULL

  method <- match.arg(method)
  dataP <- popObj(simPopObj)
  dataS <- sampleObj(simPopObj)
  # required because we do not want to change existing populaton by reference
  # thus a copy of the dataset is taken. additional variables will be added to
  # this dataset
  data_pop_o <- copy(popData(simPopObj))
  data_pop <- popData(simPopObj)
  data_sample <- sampleData(simPopObj)
  basic <- simPopObj@basicHHvars
  if(verbose){
    cat("Dimension of the population:\n")
    print(dim(data_pop))
    cat("Dimension of the sample:\n")
    print(dim(data_sample))
  }
  if ( any(additional %in% colnames(data_pop)) ) {
    stop("variables already exist in the population!\n")
  }

  if ( (length(regModel)==1|class(regModel)=="formula") & length(additional)>1 ) {
    if(class(regModel)=="formula"){
      regModelL <- list()
      for(i in seq_along(additional)){
        regModelL[[i]] <- regModel
      }
      regModel <- regModelL
    }else if ( regModel %in% c("available","basic") ) {
      regModel <- rep(regModel, length(additional))
    }
  }
  if (!is.null(regModel) ) {
    if ( class(regModel)=="formula" ) {
      regModel <- list(regModel)
    }
  }
  if ( method=="distribution" ) {
    if ( is.null(regModel) ) {
      regModel <- "basic"
    } else {
      if ( length(regModel)!=1 ) {
        stop("For method 'distribution' parameter regModel must bei either NULL, a formula or
                'basic' or 'available'!\n")
      }
    }
  } else {
    if ( is.null(regModel) ) {
      regModel <- rep("basic", length(additional))
    }
  }
  # parameters for parallel computing
  if(by=="strata"){
    curStrata <- dataS@strata
  }else{
    curStrata <- by
  }
  if(!curStrata%in%colnames(data_sample)){
    stop(curStrata," is defined as by variable, but not in the sample data set.")
  }
  if(!curStrata%in%colnames(data_pop)){
    stop(curStrata," is defined as by variable, but not in the population data set.")
  }
  nr_strata <- length(levels(data_sample[[curStrata]]))
  data_sample[,curStrata,with=FALSE]
  pp <- parallelParameters(nr_cpus=nr_cpus, nr_strata=nr_strata)
  parallel <- pp$parallel
  nr_cores <- pp$nr_cores
  have_win <- pp$have_win; rm(pp)

  ##### initializations
  if ( !missing(seed) ) {
    set.seed(seed)  # set seed of random number generator
  }


  # check arguments to account for structural zeros
  if ( length(additional) == 1 ) {
    if ( !(length(limit) == 1 && isTRUE(names(limit) == additional)) ) {
      limit <- list(limit)
      names(limit) <- additional
    }
    if ( !(length(censor) == 1 && isTRUE(names(censor) == additional)) ) {
      censor <- list(censor)
      names(censor) <- additional
    }
  }

  # list indStrata contains the indices of dataP split by strata
  N <- nrow(data_pop)
  indStrata <- split(1:N, data_pop[[curStrata]])

  ##### simulation
  if ( method == "distribution" ) {
    regInput <- regressionInput(simPopObj, additional=additional[1], regModel=regModel[1])
    predNames <- setdiff(regInput[[1]]$predNames, c(dataS@hhsize, curStrata))

    # observations with missings are excluded from simulation
    # fix #31?
    exclude <- getExclude(data_sample[,c(additional,predNames),with=F])
    if ( length(exclude) > 0 ) {
      data_sample <- data_sample[-exclude,]
    }
    data_sample <- checkFactor(data_sample, c(curStrata, predNames, additional))
    data_pop <- checkFactor(data_pop, c(curStrata, predNames))

    params <- list()
    params$grid <- expand.grid(lapply(data_sample[,additional, with=F], levels))
    params$additional <- additional
    params$basic <- predNames
    if(verbose) cat("Variables used for method 'distribution':\n"); print(params$basic)
    params$w <- dataS@weight

    if ( parallel ) {
      # windows
      if ( have_win ) {
        cl <- makePSOCKcluster(nr_cores)
        registerDoParallel(cl,cores=nr_cores)
        values <- foreach(x=levels(data_sample[[curStrata]]), .options.snow=list(preschedule=FALSE)) %dopar% {
          generateValues_distribution(
              dataSample=data_sample[data_sample[[curStrata]] == x,],
              dataPop=data_pop[indStrata[[x]], params$basic, with=F], params
          )
        }
        stopCluster(cl)
      }
      # linux/max
      if ( !have_win ) {
        values <- mclapply(levels(data_sample[[curStrata]]), function(x) {
              generateValues_distribution(
                  dataSample=data_sample[data_sample[[curStrata]] == x,],
                  dataPop=data_pop[indStrata[[x]], params$basic, with=F], params)
            }, mc.cores=nr_cores)
      }
    } else {
      values <- lapply(levels(data_sample[[curStrata]]), function(x) {
            generateValues_distribution(
                dataSample=data_sample[data_sample[[curStrata]] == x,c(additional,params$basic),with=F],
                dataPop=data_pop[indStrata[[x]], params$basic, with=F], params)
          })
    }
    values <- do.call("rbind", values)
    values <- values[unlist(indStrata),,drop=F]

    ## add new categorical variables to data set and return
    for ( i in additional ) {
      data_pop_o[,newAdditionalVarible:= values[,i]]
      setnames(data_pop_o,"newAdditionalVarible",i)
    }
    simPopObj@pop@data <- data_pop_o
    return(invisible(simPopObj))
  }

  # any other method
  counter <- 0
  for ( i in additional ) {
    counter <- counter+1
    if(verbose) cat(paste0("Simulating variable '",i,"'.\n"))
    if(length(regModel)>1){
      curRegModel <- regModel[counter]
    }else{
      curRegModel <- regModel
    }
    regInput <- regressionInput(simPopObj, additional=additional[counter], regModel=curRegModel)
    predNames <- setdiff(regInput[[1]]$predNames, c(dataS@hhsize, curStrata))

    # observations with missings are excluded from simulation
    exclude <- getExclude(data_sample[,c(additional,predNames),with=F])
    if ( length(exclude) > 0 ) {
      sampWork <- data_sample[-exclude,]
    } else {
      sampWork <- data_sample
    }

    # variables are coerced to factors
    sampWork <- checkFactor(sampWork, unique(c(curStrata, predNames, additional)))
    data_pop <- checkFactor(data_pop_o, unique(c(curStrata, predNames)))

    # components of multinomial model are specified
    levelsResponse <- levels(sampWork[[i]])

    # simulation of variables using a sequence of multinomial models
    if ( method == "multinom" ) {
      formula.cmd <- paste(i, "~", paste(predNames, collapse = " + "))
      formula.cmd <- paste0("suppressWarnings(multinom(", formula.cmd,
          ", weights=", dataS@weight, ", data=dataSample, trace=FALSE",
          ", maxit=",maxit, ", MaxNWts=", MaxNWts,"))")
      if(verbose) cat("we are running the following multinom-model:\n")
      if(verbose) cat(strwrap(cat(gsub("))",")",gsub("suppressWarnings[(]","",formula.cmd)),"\n"), 76), sep = "\n")
    }else if ( method == "ctree" ) {
      # simulation via recursive partitioning and regression trees
      formula.cmd <- paste(i, "~", paste(predNames, collapse = " + "))
      formula.cmd <- paste("suppressWarnings(ctree(", formula.cmd, ", weights=as.integer(dataSample$", dataS@weight, "), data=dataSample))", sep="")
      if(verbose) cat("we are running recursive partitioning:\n")
      if(verbose) cat(strwrap(cat(gsub("))",")",gsub("suppressWarnings[(]","",formula.cmd)),"\n"), 76), sep = "\n")
    }else if ( method == "cforest" ) {
      # simulation via recursive partitioning and regression trees
      formula.cmd <- paste(i, "~", paste(predNames, collapse = " + "))
      formula.cmd <- paste("suppressWarnings(cforest(", formula.cmd, ", weights=as.integer(dataSample$", dataS@weight, "), data=dataSample))", sep="")
      if(verbose) cat("we are running recursive partitioning:\n")
      if(verbose) cat(strwrap(cat(gsub("))",")",gsub("suppressWarnings[(]","",formula.cmd)),"\n"), 76), sep = "\n")
    }
    #if ( method == "naivebayes" ) {
    #  formula.cmd <- paste(i, "~", paste(predNames, collapse = " + "))
    #  formula.cmd <- paste("naiveBayes(", formula.cmd, ", data=dataSample, usekernel=TRUE)", sep="")
    #}

    # check if population data contains factor levels that do not exist
    # in the sample
    newLevels <- lapply(predNames, function(nam) {
          levelsS <- levels(sampWork[[nam]])
          levelsP <- levels(data_pop[[nam]])
          levelsP[!(levelsP %in% levelsS)]
        })
    hasNewLevels <- sapply(newLevels, length) > 0
    excludeLevels <- any(hasNewLevels)

    # generate values of new variable
    params <- list()
    params$method <- method
    params$cur.var <- i
    params$excludeLevels <- excludeLevels
    params$hasNewLevels <- hasNewLevels
    params$newLevels <- newLevels
    params$w <- dataS@weight
    params$formula.cmd <- formula.cmd
    params$eps <- eps
    params$limit <- limit
    params$censor <- censor
    params$levelsResponse <- levelsResponse

    # windows
    if ( parallel ) {
      if ( have_win ) {
        cl <- makePSOCKcluster(nr_cores)
        registerDoParallel(cl,cores=nr_cores)
        values <- foreach(x=levels(data_sample[[curStrata]]), .options.snow=list(preschedule=FALSE)) %dopar% {
          generateValues(
              dataSample=sampWork[sampWork[[curStrata]] == x,],
              dataPop=data_pop[indStrata[[x]], predNames, with=F], params
          )
        }
        stopCluster(cl)
      }
      # linux/mac
      if ( !have_win) {
        values <- mclapply(levels(data_sample[[curStrata]]), function(x) {
              generateValues(
                  dataSample=sampWork[sampWork[[curStrata]] == x,],
                  dataPop=data_pop[indStrata[[x]], predNames, with=F], params
              )
            }, mc.cores=nr_cores)
      }
    } else {
      values <- lapply(levels(data_sample[[curStrata]]), function(x) {
            generateValues(
                dataSample=sampWork[sampWork[[curStrata]] == x,],
                dataPop=data_pop[indStrata[[x]], predNames, with=F], params
            )
          })
#      print(str(values))
    }
    values <- factor(unsplit(values, data_pop[[curStrata]]), levels=levelsResponse)
    ## add new categorical variable to data set
#    print(str(values))
#    print(length(values))
#    print(length(data_pop_o[[i]]))
#    print(i)
    data_pop_o[[i]] <- values
    simPopObj@pop@data <- data_pop_o
  }
  invisible(simPopObj)
}
