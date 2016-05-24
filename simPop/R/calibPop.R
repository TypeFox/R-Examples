# required for parallelisation
calcFinalWeights <- function(data0, totals0, params) {
  # TODO: speed-improvements
  # each row of output contains indices of a specific margin
  indices_by_constraint <- function(data0, totals0, parameter) {
    # recode to character and set NA to "NA"
    myfun <- function(x) {
      x <- as.character(x)
      ii <- which(is.na(x))
      if ( length(ii) > 0 ) {
        x[ii] <- "NA"
      }
      x
    }

    out <- matrix(NA, nrow=nrow(totals0), ncol=nrow(data0))
    totals0 <- as.data.frame(totals0)
    totals0[,parameter] <- apply(totals0[,parameter], 2, myfun)
    for ( x in parameter ) {
      data0[[x]] <- myfun(data0[[x]])
    }
    for (i in 1:nrow(totals0) ) {
      ex <- paste("out[i,] <- as.integer(", sep="")
      for ( k in seq_along(parameter) ) {
        if ( k > 1 ) {
          ex <- paste(ex, "&", sep="")
        }
        ex <- paste(ex, 'data0[,"',parameter[k],'", with=F] =="',totals0[i, parameter[k]],'"', sep="")
      }
      ex <- paste(ex, ")", sep="")
      eval(parse(text=ex))
    }
    out
  }

  inp <- indices_by_constraint(data0, totals0, params$parameter)
  current_totals <- as.numeric(totals0$N)
  weights <- as.integer(data0[[params$weight]])

  hh_info <- list()
  hh_info$hh_ids <- as.integer(data0[[params$hhid]])
  hh_info$hh_head <- rep(0L, nrow(data0))
  index <- which(sapply(data0[[params$pid]], function(x) {
    unlist(strsplit(x, "[.]"))[2] }) == "1")
  hh_info$hh_head[index] <- 1L
  hh_info$hh_size <- as.integer(data0[[params$hhsize]])
  hh_info$median_hhsize <- median(hh_info$hh_size[hh_info$hh_head==1], na.rm=TRUE)

  w <- .Call("simPop_calibPop_work", inp=inp, totals=current_totals,
    weights=weights, hh_info=hh_info, params=params, package="simPop")
  invisible(w)
}



#' Calibration of 0/1 weights by Simulated Annealing
#'
#' A Simulated Annealing Algorithm for calibration of synthetic population data
#' available in a \code{\linkS4class{simPopObj}}-object. The aims is to find,
#' given a population, a combination of different households which optimally
#' satisfy, in the sense of an acceptable error, a given table of specific
#' known marginals. The known marginals are also already available in slot
#' 'table' of the input object 'inp'.
#'
#' Calibrates data using simulated annealing. The algorithm searches for a
#' (near) optimal combination of different households, by swaping housholds at
#' random in each iteration of each temperature level. During the algorithm as
#' well as for the output the optimal (or so far best) combination will be
#' indicated by a logical vector containg only 0s (not inculded) and 1s
#' (included in optimal selection). The objective function for simulated
#' annealing is defined by the sum of absolute differences between target
#' marginals and synthetic marginals (=marginals of synthetic dataset). The sum
#' of target marginals can at most be as large as the sum of target marginals.
#' For every factor-level in \dQuote{split}, data must at least contain as many
#' entries of this kind as target marginals.
#'
#' Possible donors are automatically generated within the procedure.
#'
#' The number of cpus are selected automatically in the following manner. The
#' number of cpus is equal the number of strata. However, if the number of cpus
#' is less than the number of strata, the number of cpus - 1 is used by
#' default.  This should be the best strategy, but the user can also overwrite
#' this decision.
#'
#' @name calibPop
#' @docType data
#' @param inp an object of class \code{\linkS4class{simPopObj}} with slot
#' 'table' being non-null! (see \code{\link{addKnownMargins}}.
#' @param split given strata in which the problem will be split. Has to
#' correspond to a column population data (slot 'pop' of input argument 'inp')
#' . For example \code{split = c("region")}, problem will be split for
#' different regions. Parallel computing is performed automatically, if
#' possible.
#' @param temp starting temperatur for simulated annealing algorithm
#' @param eps.factor a factor (between 0 and 1) specifying the acceptance
#' error. For example eps.factor = 0.05 results in an acceptance error for the
#' objective function of 0.05*sum(totals)
#' @param maxiter maximum iterations during a temperature step.
#' @param temp.cooldown a factor (between 0 and 1) specifying the rate at which
#' temperature will be reduced in each step.
#' @param factor.cooldown a factor (between 0 and 1) specifying the rate at
#' which the number of permutations of housholds, in each iteration, will be
#' reduced in each step.
#' @param min.temp minimal temperature at which the algorithm will stop.
#' @param nr_cpus if specified, an integer number defining the number of cpus
#' that should be used for parallel processing.
#' @param verbose boolean variable; if TRUE some additional verbose output is
#' provided, however only if \code{split} is NULL. Otherwise the computation is
#' performed in parallel and no useful output can be provided.
#' @return Returns an object of class \code{\linkS4class{simPopObj}} with an
#' updated population listed in slot 'pop'.
#' @author Bernhard Meindl, Johannes Gussenbauer and Matthias Templ
#' @keywords datasets
#' @export
#' @examples
#' data(eusilcS) # load sample data
#' data(eusilcP) # population data
#' \dontrun{
#' ## approx. 20 seconds computation time
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
#' simPop <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
#' simPop <- simCategorical(simPop, additional=c("pl030", "pb220a"), method="multinom", nr_cpus=1)
#'
#' # add margins
#' margins <- as.data.frame(
#'   xtabs(rep(1, nrow(eusilcP)) ~ eusilcP$region + eusilcP$gender + eusilcP$citizenship))
#' colnames(margins) <- c("db040", "rb090", "pb220a", "freq")
#' simPop <- addKnownMargins(simPop, margins)
#' }
#'
#' # apply simulated annealing
#' \dontrun{
#' ## long computation time
#' simPop_adj <- calibPop(simPop, split="db040", temp=1, eps.factor=0.1)
#' }
calibPop <- function(inp, split, temp = 1, eps.factor = 0.05, maxiter=200,
  temp.cooldown = 0.9, factor.cooldown = 0.85, min.temp = 10^-3, nr_cpus=NULL, verbose=FALSE) {

  if ( class(inp) != "simPopObj" ) {
    stop("argument 'inp' must be of class 'simPopObj'!\n")
  }

  if ( length(split) > 1 ) {
    split <- split[1]
    warning("only first variable will be used to divide the population into strata")
  }

  x <- NULL
  data <- popData(inp)
  hid <- popObj(inp)@hhid
  pid <- popObj(inp)@pid
  hhsize <- popObj(inp)@hhsize
  totals <- tableObj(inp)
  parameter <- colnames(totals)[-ncol(totals)]

  if ( !is.null(split) ) {
    verbose <- FALSE
    if ( !split %in% colnames(data) ) {
      stop("variable specified in argument 'split' must be a column in synthetic population (slot 'data' of argument 'inp')!\n")
    }
    if ( !split %in% parameter ) {
      stop("variable specified in argument 'split' must be a column in slot 'table' of argument 'inp'!\n")
    }
  } else {
    data$tmpsplit <- 1
    totals$tmpsplit <- 1
    split <- "tmpsplit"
  }

  params <- list()
  params$temp <- as.numeric(temp)[1]
  params$eps_factor = as.numeric(eps.factor)[1]
  params$maxiter = as.integer(maxiter)[1]
  params$temp_cooldown = as.numeric(temp.cooldown)[1]
  params$factor_cooldown = as.numeric(factor.cooldown)[1]
  params$min_temp = as.numeric(min.temp)[1]
  params$verbose <- ifelse(verbose, 1L, 0L)
  params$parameter <- parameter
  params$weight <- popObj(inp)@weight
  params$hhid <- hid
  params$pid <- pid
  params$hhsize <- hhsize

  # generate donors
  data2 <- sampHH(data, sizefactor=2, hid=hid, strata=split, hsize=hhsize)
  data2[[params$weight]] <- 0
  data2 <- data2[, which(!grepl(hid, colnames(data2))), with=FALSE]
  cn <- colnames(data2)
  cn[length(cn)] <- hid
  setnames(data2, cn)
  data2 <- data2[,match(colnames(data), cn), with=FALSE]
  data <- rbind(data, data2)

  # parameters for parallel computing
  nr_strata <- length(unique(data[[split]]))
  pp <- parallelParameters(nr_cpus=nr_cpus, nr_strata=nr_strata)
  parallel <- pp$parallel
  nr_cores <- pp$nr_cores
  have_win <- pp$have_win; rm(pp)

  ii <- match(parameter, colnames(data))
  data[,ii] <- data[,lapply(.SD,as.factor),.SDcols=parameter]
  rm(ii)

  ## split problem by "split"-factor
  setkeyv(data, split)
  setkeyv(totals, split)
  split.number <- unique(data[,split,with=F])

  if ( parallel ) {
    # windows
    if ( have_win ) {
      cl <- makePSOCKcluster(nr_cores)
      registerDoParallel(cl,cores=nr_cores)
      final_weights <- foreach(x=1:nrow(split.number), .options.snow=list(preschedule=TRUE)) %dopar% {
        calcFinalWeights(
          data0=data[split.number[x]],
          totals0=totals[which(totals[,split,with=F]==as.character(split.number[x][[split]])),],
          params=params
        )
      }
      stopCluster(cl)
    }else if ( !have_win ) {# linux/mac
      final_weights <- mclapply(1:nrow(split.number), function(x) {
        calcFinalWeights(
          data0=data[split.number[x]],
          totals0=totals[which(totals[,split,with=F]==as.character(split.number[x][[split]])),],
          params=params)
      },mc.cores=nr_cores)
    }
  } else {
    final_weights <- lapply(1:nrow(split.number), function(x) {
      calcFinalWeights(
        data0=data[split.number[x]],
        totals0=totals[which(totals[,split,with=F]==as.character(split.number[x][[split]])),],
        params=params)
    })
  }
  # return dataset with new weights
  data$new.weights <- as.integer(unlist(final_weights))
  data <- data[data$new.weights==1,]
  data[[inp@pop@weight]] <- data$new.weights
  data$new.weights <- NULL
  inp@pop@data <- data
  invisible(inp)
}
NULL

