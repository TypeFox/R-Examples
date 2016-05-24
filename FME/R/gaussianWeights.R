## This function adds a gaussian weighing function to a dataset in long format. By adding weights in a column called err. At the same time the dataset is
## extended to replicate each observation as many times as needed (see details).
##
## Parameters:
## obs         : data set in long format as is typically used in modCost
## x           : name of the independent variable in the dataset
## xmodel      : unique times at which model output is produced
## spread      : standard deviation used to calculate the weights from a normal distribution (stdev = standard deviation of that normal distribution)
## ordering    : list of variables in obs used for grouping during scale factor calculation (cf. tapply)
## weight      : scaling factor of modCost function (sd, mean, or none)

gaussianWeights <- function(obs,                            # data set in long format as is typically used in modCost (e.g. Meso1 in long format)
                            x = x,                          # name of the independent variable (e.g. "day" in Meso1)
                            y = y,
                            xmodel,                         # unique times at which model output is produced
                            spread,                         # standard deviation used to calculate the weights from a normal distribution (stdev = standard deviation of that normal distribution)
                            weight      = "none",           # scaling factor of modCost function (sd, mean, or none)
                            aggregation = x,
                            ordering                        # list of variables in obs used for grouping during scale factor calculation (cf. tapply)
                           )
{

  varlist <- seq_along(names(obs))
  names(varlist) <- names(obs)

  xname <- as.character(substitute(x))
  yname <- as.character(substitute(y))
  aggname <- as.character(substitute(aggregation))
  if(missing(ordering)) orderingnames <- c(aggname,xname)
  else orderingnames <- as.character(substitute(ordering))
  if(!(aggname %in% orderingnames)) orderingnames <- c(aggname,orderingnames)
  if(!(xname %in% orderingnames))  orderingnames <- c(orderingnames,xname)

  # Exception/error handling
  if(!(xname %in% names(obs)))
    stop("Independent variable in observational dataset not known")
  if(!(yname %in% names(obs)))
    stop("Dependent variable in observational dataset not known")
  if(!(aggname %in% names(obs)))
    stop("Aggregation factor unknown")
  if(length(union(orderingnames,names(varlist))) > length(varlist))
    stop("At least one ordering variable not defined in observations")

  xind <- varlist[names(varlist) == substitute(x)]
  yind <- varlist[names(varlist) == substitute(y)]
  aggind <- varlist[names(varlist) == substitute(aggregation)]


  if(missing(xmodel)) {
    warning("No independent variable for output specified: using independent variable of input instead")
    xmodel <- unique(unlist(obs[,xind]))
  }

  if(!(weight %in% c("mean","sd","none")))
    stop("Unknown weighing factor")

  # Selection of indices of relevant model times
  obsTimesSelect <- lapply(as.list(unique(obs[,xind])),
                           FUN = function(t) which(abs(xmodel - t) == 0)
                           )
  obsRangeSelect <- lapply(as.list(unique(obs[,xind])),
                           FUN = function(t) which(abs(xmodel - t) <= 3*spread)
                           )

  # Collection of actual model times
  obsTimes   <- lapply(obsTimesSelect, function(t) xmodel[t])
  obsRanges  <- lapply(obsRangeSelect, function(t) xmodel[t])
  obsWeights <- obsRanges # just to copy data structure

  # Calculation of inverse weights
  for(i in 1:length(obsTimes)) {
     obsWeights[[i]] <- dnorm(obsWeights[[i]],mean=obsTimes[[i]],sd=spread)
     obsWeights[[i]] <- obsWeights[[i]]/sum(obsWeights[[i]]) # rescaling to sum up to 1
     obsTimes[[i]] <- rep(obsTimes[[i]],times=length(obsWeights[[i]])) # replication of original observation times
  }

  # Construction of error data frame
  errors <- data.frame(time = unlist(obsTimes),
                       err = 1/unlist(obsWeights),
                       xmodel=unlist(obsRanges)
                      )
  # Merging with observations to obtain an extended observational dataset
  fullObs <- merge(obs,errors,by.x=xind,by.y="time",all.y=T)
  fullObs[,xname] <- fullObs[,"xmodel"]   # observational times are replaced by model times
  fullObs[,"xmodel"] <- NULL                             # the latter column is redundant, hence removed

  if(weight != "none") {
    aggregationFunction <- tapply(fullObs[,yname],INDEX=fullObs[,rev(aggname)],weight,na.rm=T)
    factorLevels <- dimnames(aggregationFunction)
    names(factorLevels) <- rev(aggname)
    aggregationFactor <- data.frame(expand.grid(factorLevels),
                                    scaling  = abs(array(aggregationFunction))
                                   )
    fullObs <- merge(fullObs,aggregationFactor)
    fullObs$err <- fullObs$err * fullObs$scaling
    fullObs$scaling <- NULL
  }

  fullObs <- eval(parse(text=paste("fullObs[order(",
                                   paste("fullObs$",orderingnames,sep="",collapse=","),
                                   "),]",
                                   collapse=""
                                  )
                       )
                 )

  firstColumns <- c(aggname[1],xname,yname,"err")
  lastColumns  <- setdiff(names(fullObs),firstColumns)
  fullObs <- data.frame(fullObs[,firstColumns],fullObs[,lastColumns])
  names(fullObs) <- c(firstColumns,lastColumns)

  return(fullObs)
}
