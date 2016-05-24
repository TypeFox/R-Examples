#' Simulate the household structure of population data
#' 
#' Simulate basic categorical variables that define the household structure
#' (typically variables such as household ID, age and gender) of population
#' data by resampling from survey data.
#' 
#' @name simStructure
#' @param dataS an object of class \code{dataObj} containing household survey
#' data that is usually generated with \code{\link{specifyInput}}.
#' @param method a character string specifying the method to be used for
#' simulating the household sizes.  Accepted values are \code{"direct"}
#' (estimation of the population totals for each combination of stratum and
#' household size using the Horvitz-Thompson estimator), \code{"multinom"}
#' (estimation of the conditional probabilities within the strata using a
#' multinomial log-linear model and random draws from the resulting
#' distributions), or \code{"distribution"} (random draws from the observed
#' conditional distributions within the strata).
#' @param basicHHvars a character vector specifying important variables for the
#' household structure that need to be available in \code{dataS}. Typically
#' variables such as age or sex may be used.
#' @param seed optional; an integer value to be used as the seed of the random
#' number generator, or an integer vector containing the state of the random
#' number generator to be restored.
#' @return An object of class \code{simPopObj} containing the simulated
#' population household structure as well as the underlying sample that was
#' provided as input.
#' @note The function \code{\link{sample}} is used, which gives results
#' incompatible with those from < 2.2.0 and produces a warning the first time
#' this happens in a session.
#' @author Bernhard Meindl and Andreas Alfons
#' @seealso \code{\link{simCategorical}}, \code{\link{simContinuous}},
#' \code{\link{simComponents}}, \code{\link{simEUSILC}}
#' @keywords datagen
#' @export
#' @examples
#' 
#' data(eusilcS)
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
#' eusilcP <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
#' class(eusilcP)
#' eusilcP
#' 
simStructure <- function(dataS, method=c("direct", "multinom", "distribution"), basicHHvars, seed=1) {
  if ( !class(dataS) == "dataObj" ) {
    stop("Error. Please provide the input sample in the required format.\n
      It must be an object of class 'dataObj' that can be created using function specifyInput()!\n")
  }
  if ( dataS@ispopulation ) {
    cat("the structure is created from a population! \n")
    #stop("dataS must contain sample information!\n")
  }

  if ( is.null(basicHHvars) ) {
    stop("please provide valid variables defining the household structure!\n")
  }

  if ( !all(basicHHvars %in% colnames(dataS@data)) ) {
    stop("not all variables listed in argument 'basicHHvars' is available in the sample input data!\n")
  }

  ##### initializations
  if ( !missing(seed) ) {
    set.seed(seed)
  }

  method <- match.arg(method)

  # order sample (needed later on)
  if ( is.null(dataS@pid) ) {
    setkeyv(dataS@data, dataS@hhid)
  } else {
    setkeyv(dataS@data, dataS@pid)
  }

  # extract variables
  hid <- dataS@data[[dataS@hhid]]
  if ( dataS@ispopulation ) {
    w <- rep(1,length(hid))
  }else{
    w <- dataS@data[[dataS@weight]]  
  }
  
  hsize <- dataS@data[[dataS@hhsize]]
  if ( !is.null(dataS@strata) ) {
    strata <- factor(dataS@data[[dataS@strata]])
  } else {
    strata <- factor(rep(1, nrow(dataS@data)))
  }

  ##### set up household structure

  # generate variables on household level (indicated by H)
  # these variables have the same value for all household members, hence
  # we take the first occurence (faster than aggregation with 'unique')
  hfirst <- !duplicated(hid)
  dataH <- data.frame(hsize=as.factor(hsize)[hfirst], strata=strata[hfirst])
  wH <- w[hfirst]
  households <- tableWt(dataH, wH)  # expected number of households

  if ( method != "direct" ) {
    ### simulation of household size
    ls <- colnames(households)  # strata levels
    NH <- apply(households, 2, sum)  # number of households by stratum
    if ( method == "multinom" ) {
      empty <- which(households == 0)
      # TODO: allow setting some more arguments
      mod <- suppressWarnings(multinom(hsize~strata, weights=wH, data=dataH, trace=FALSE))
      newdata <- data.frame(strata=ls)
      rownames(newdata) <- ls
      probs <- as.matrix(predict(mod, newdata=newdata, type="probs"))
      # # experimental: boosting factor for small probabilities
      # fraction <- nrow(dataH)/nrow(dataPH)
      # boost <- 1/sqrt(fraction)
      # probs[probs < fraction] <- pmin(fraction, boost*probs[probs < fraction])
      hsizePH <- unlist(lapply(ls, function(l) spSample(NH[l], probs[l,])))
    } else if ( method == "distribution" ) {
      hsizePH <- unlist(lapply(ls, function(l) spSample(NH[l], households[, l])))
    }
    dataPH <- data.frame(hsize=as.factor(hsizePH), strata=factor(rep(ls, times=NH), levels=ls, ordered=is.ordered(strata)))
    households <- tableWt(dataPH)  # recompute number of households
  }

  ### simulation of household structure

  # all combinations of strata and household size (that do occur
  # in the sample) are stored in data frame 'grid'
  grid <- expand.grid(hsize=rownames(households), strata=colnames(households), stringsAsFactors=FALSE)
  if ( method != "multinom" ) {
    grid <- grid[households != 0,]  # remove those that do not occur
  }
  ncomb <- nrow(grid)

  # indices of households by strata and household size, with combinations
  # that do not occur in the sample left out (drop = TRUE)
  split <- split(1:nrow(dataH), dataH, drop = (method != "multinom"))

  # to be on the safe side, names are reconstructed from 'grid'
  # and the list 'split' is sorted according to those names
  nam <- apply(as.matrix(grid), 1, paste, collapse=".")
  split <- split[nam]
  if ( method == "multinom" ) {
    # combinations with strata that do not occur in the sample borrow from
    # indices of all households in the sample with the same household size
    donor <- grid[empty, 1]
    split[empty] <- split(1:nrow(dataH), dataH$hsize)[donor]
  }
  # for each stratum, draw from original sample
  numbers <- lapply(1:ncomb, function(i) {
    n <- households[grid[i, 1], grid[i, 2]]
    w <- wH[split[[i]]]
    p <- w / sum(w)  # probability weights
    spSample(n, p)
  })

  ### generation of the household structure for the population

  # the sampled household numbers in list 'numbers' are transformed to
  # indices of data frame 'dataS' and stored in vector 'indices'
  hidH <- hid[hfirst]
  indices <- lapply(1:ncomb, function(i) {
    pn <- which(hid %in% hidH[split[[i]]])
    pnM  <- matrix(pn, nrow=as.numeric(grid[i, 1]))
    c(pnM[, numbers[[i]]])
  })
  indices <- unlist(indices)

  # new household IDs are created for sampling from population
  if ( method == "direct" ) {
    # household IDs are constructed from the cells in the table of
    # population households
    upper <- cumsum(households[households != 0])
    lower <- c(1, upper[-ncomb] + 1)
    hidNew <- lapply(1:ncomb, function(i) {
      rep(lower[i]:upper[i], each=as.numeric(grid[i,1]))
    })
  } else {
    # household IDs are constructed from the population household data
    hidNew <- split(1:nrow(dataPH), dataPH, drop=(method == "distribution"))
    hidNew <- lapply(1:ncomb, function(i) {
      rep(hidNew[[i]], each=as.numeric(grid[i,1]))
    })
  }
  hidNew <- unlist(hidNew)

  ##### return simulated population household structure
  # it is much faster to generate replications for each variable
  # individually than for the data frame as a whole

  # build command
  if ( !is.null(basicHHvars) ) {
    if ( !all(basicHHvars %in% colnames(dataS@data)) ) {
      stop("please specify valid variable names for argument 'basicHHvars'!\n")
    }
    expr <- paste(basicHHvars, "=dataS@data[[\"", basicHHvars,"\"]][indices]", sep="")
    expr <- paste(",", expr, collapse="")
  } else {
    expr <- ""
  }

  cmd <- paste0("data.table(", dataS@hhid, "=hidNew, ", dataS@hhsize, "=hsize[indices]")
  if ( is.null(dataS@strata) ) {
    cmd <- paste0(cmd,")")
  } else {
    if ( dataS@strata %in% basicHHvars ) {
        cmd <- paste0(cmd, expr,")")
    } else {
      if ( method == "direct" ) {
        cmd <- paste0(cmd, ", ",dataS@strata,"=strata[indices]", expr,")")
      } else {
        cmd <- paste0(cmd, ", ",dataS@strata,"=dataPH$strata[hidNew]", expr, ")")
      }      
    }
  }
  # evaluate command and return result
  dataP <- eval(parse(text=cmd))
  setkeyv(dataP, dataS@hhid)

  sizes <- dataP[,.N, by=key(dataP)]
  pid <- paste(dataP[[dataS@hhid]], ".",unlist(sapply(sizes[["N"]], function(x) seq(1,x))), sep="")
  dataP[[dataS@pid]] <- pid
  dataP$weight <- 1

  pop <- new("dataObj", data=dataP, hhid=dataS@hhid, hhsize=dataS@hhsize, pid=dataS@pid, strata=dataS@strata, weight="weight", ispopulation=TRUE)
  out <- new("simPopObj", sample=dataS, table=NULL, pop=pop)
  out@basicHHvars <- basicHHvars
  invisible(out)
}
