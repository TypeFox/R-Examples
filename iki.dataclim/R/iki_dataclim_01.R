## make zoo (S3) useable for S4
setOldClass("zoo")

setClass("dataclim",
         slots=c(basePeriod="numeric", flagged="list", flaggedData="zoo",
             data="zoo",
             monthlyAvg="zoo", annualAvg="zoo",
             monthlyMin="zoo", annualMin="zoo",
             monthlyMax="zoo", annualMax="zoo",
             convertFlaggedToNA="logical")
         )

.qcDailyPrec <- function(x){
  if (is.zoo(x) && length(x)>0){
    ## check precip in 0...300 mm
    inds <- which(findInterval(x, c(0.0, 300)) != 1)

    ## check repetiveness
    rle <- rle(coredata(x))
    ## ... > 1mm
    rep1 <- which(rle$value > 1 & rle$length >= 10)
    if (length(rep1) > 0){
      for (i in rep1){
        start <- 1
        if (i>1)
            start <- sum(rle$length[1:(i - 1)]) + 1
        inds <- c(inds, start:sum(rle$length[1:i]))
      }
    }

    ## ... > 5mm
    rep5 <- which(rle$value > 5 & rle$length >= 5)
    if (length(rep5) > 0){
      for (i in rep5){
        start <- 1
        if (i>1)
            start <- sum(rle$length[1:(i - 1)]) + 1
        inds <- c(inds, start:sum(rle$length[1:i]))
      }
    }

    ## ... <= 1mm
    repdry <- which(rle$value <= 1 &
                    rle$length >= 14 * sd(rle$length[rle$value <= 1]))
    if (length(repdry) > 0){
      for (i in repdry){
        start <- 1
        if (i>1)
            start <- sum(rle$length[1:(i - 1)]) + 1
        inds <- c(inds, start:sum(rle$length[1:i]))
      }
    }

    return(sort(unique(inds)))
  } else {
    return(NULL)
  }
}

.qcDailyTemp <- function(tmin, tmax, tavg=.5 * (tmin + tmax),
                        basePeriod){
  tempDF <- zoo(data.frame(tmin=tmin, tmax=tmax, tavg=tavg), index(tmin))

  inds <- list()

  ## check individual temperatures
  for (v in names(tempDF)){
    ##check intervals
    inds[[v]] <-  which(!(findInterval(tempDF[, v], c(-90, 60),
                                       rightmost.closed=TRUE) == 1))
    
    ##check runs of constant values
    rle <- rle(coredata(tempDF[, v]))
    rep5 <- which(rle$length >= 5)
    if (length(rep5) > 0){
      for (i in rep5){
        start <- 1
        if (i>1)
            start <- sum(rle$length[1:(i - 1)]) + 1
        inds[[v]] <- c(inds[[v]], start:sum(rle$length[1:i]))
      }
    }
  }

  ## check >< inequalities
  bad.inequalities <- which(tempDF$tmax < tempDF$tavg |
      tempDF$tmax < tempDF$tmin |
          tempDF$tavg < tempDF$tmin)
  for (v in names(tempDF))
      inds[[v]] <- c(inds[[v]], bad.inequalities)

  ## check outliers wrt annual cycle
  ## -- define indices per julian day +- 2 days
  yday.groups <- apply(cbind(month(tempDF), mday(tempDF)), 1, paste,
                       collapse="-")
  yday.inds <- split(1:nrow(tempDF), yday.groups)
  yday.inds.window <- sapply(yday.inds, function(ii) {
    iii <- as.vector(t(outer(ii, -2:2, "+")))
    iii <- iii[iii > 0 & iii <= nrow(tempDF)]
    return(iii)
  })

  ## calc "normal" climatological bandwidth of each julian day and
  ## check each day within this bandwidth
  for (v in names(tempDF)){
    clim.avg <- sapply(yday.inds.window,
                       function(ii) mean(tempDF[ii, v], na.rm=TRUE))
    clim.sd <- sapply(yday.inds.window,
                      function(ii) sd(tempDF[ii, v], na.rm=TRUE))
    clim.upper <- clim.avg + 5 * clim.sd
    clim.lower <- clim.avg - 5 * clim.sd
    for (d in names(yday.inds)){
      d.days <- yday.inds[[d]]
      inds[[v]] <- c(inds[[v]],
                     d.days[findInterval(tempDF[d.days, v],
                                         c(clim.lower[d], clim.upper[d]),
                                         rightmost.closed=TRUE) != 1]
                     )
    } ## loop over yday.inds
  } ## loop over v
  inds <- lapply(inds, function(ii) sort(unique(ii)))
  
  return(inds)
}

setMethod(f="initialize",
          signature=signature(.Object="dataclim"),
          definition=function(.Object, 
              date, tmin, tmax, prec,
              basePeriod, convertFlaggedToNA, ...){

            ## store basePeriod & convertFlaggedToNA
            .Object@basePeriod <- basePeriod
            .Object@convertFlaggedToNA <- convertFlaggedToNA
            
            ## check validity of arguments
            if ( !( class(date) == "Date" && length(date) > 0 &&
                   all( c(length(tmin), length(tmax), length(prec)) ==
                       length(date)) ) ){
              stop("Arguments of different lengths.")
              return(NULL)
            }

            if (!all(basePeriod %in% unique(year(date)))){
              stop("basePeriod outside temporal data range.")
            }

            ## create minimum data.frame
            df <- data.frame(tmin=tmin, tmax=tmax, prec=prec)
            
            ## check additional variables
            more.vars <- list(...)
            if ( (length(more.vars) > 0) &&
                !all(unlist(sapply(more.vars, length)) == length(date)) ) {
              warning("Additional variables  of wrong length, not included.")
            } else if (length(more.vars) > 0) {
              for (v in names(more.vars))
                  df[[v]] <- more.vars[[v]]
            }
            remove(more.vars)

            ## build zoo-object
            df <- zoo(df, date)

            ## quality control
            flagged.inds <- list()
            flagged.inds$prec <- .qcDailyPrec(df$prec)

            flagged.temps <- .qcDailyTemp(tmin=df$tmin, tmax=df$tmax,
                                         basePeriod=basePeriod)
            flagged.inds$tmax <- flagged.temps$tmax
            flagged.inds$tmin <- flagged.temps$tmin
            if ("tavg" %in% names(df))
                flagged.inds$tavg <- flagged.temps$tavg

            .Object@flaggedData <- df[unique(unlist(flagged.inds)), ]
            .Object@flagged <- flagged.inds

            ## if selected, set entries of bad quality to NA
            if (convertFlaggedToNA){
              for (v in names(df))
                  df[, v][flagged.inds[[v]]] <- NA
            }

            ##pad with NAs to have complete years
            complete.year.days <-
                seq(as.Date(paste(year(df)[1],"-01-01",sep="")),
                    as.Date(paste(year(df)[nrow(df)],"-12-31",sep="")),
                    by="days")
            complete.df <-
                zoo(matrix(NA, ncol=ncol(df),
                           nrow=length(complete.year.days)),
                    complete.year.days)
            complete.df[index(complete.df) %in% index(df),] <- df
            names(complete.df) <- names(df)
            df <- complete.df

            .Object@data <- df

            ## calculate annual and monthly averages/sums
            monthly <- aggregate(df, as.yearmon(time(df)),
                                 .monthlyAgg, "mean", na.rm=TRUE)
            monthly$prec <- monthly$prec * days_in_month(index(monthly))
            annual <- aggregate(monthly, year(monthly),
                                .annualAgg, "mean", na.rm=TRUE)
            annual$prec <- annual$prec * 12

            .Object@monthlyAvg <- monthly
            .Object@annualAvg <- annual

            ## calculate annual and monthly minimum
            monthly <- aggregate(df, as.yearmon(time(df)),
                                 .monthlyAgg, "min", na.rm=TRUE)
            monthly$prec <- monthly$prec * days_in_month(index(monthly))
            annual <- aggregate(monthly, year(monthly),
                                .annualAgg, "min", na.rm=TRUE)
            annual$prec <- annual$prec * 12

            .Object@monthlyMin <- monthly
            .Object@annualMin <- annual

            ## calculate annual and monthly maximum
            monthly <- aggregate(df, as.yearmon(time(df)),
                                 .monthlyAgg, "max", na.rm=TRUE)
            monthly$prec <- monthly$prec * days_in_month(index(monthly))
            annual <- aggregate(monthly, year(monthly),
                                .annualAgg, "max", na.rm=TRUE)
            annual$prec <- annual$prec * 12

            .Object@monthlyMax <- monthly
            .Object@annualMax <- annual

            ## call inspector
            validObject(.Object)
            return(.Object)
          } )

.monthlyAgg <- function(x, FUN, ...){
  FUN <- match.fun(FUN)
  if (length(x) < 25 || sum(is.na(x)) > 3){
    return(NA)
  } else {
    return(FUN(x, ...))
  }
}

.annualAgg <- function(x, FUN, ...){
  FUN <- match.fun(FUN)
  if (length(x) != 12 || sum(is.na(x)) > 0){
    return(NA)
  } else {
    return(FUN(x, ...))
  }
}

      

createDataclim <- function(date=NULL, tmin=NULL, tmax=NULL, prec=NULL,
                        basePeriod=1961:1990, convertFlaggedToNA=TRUE, ...){
  if ( !( length(date) > 0 &&
         all( c(length(tmin), length(tmax), length(prec)) %in%
             length(date)) )  &&
         all( c(is.numeric(tmin), is.numeric(tmax), is.numeric(prec)) )  ){
    print("Arguments of different lengths.")
    return(NULL)
  }
  
  myDataclim <- new("dataclim", date=date, tmin=tmin, tmax=tmax, prec=prec,
                    basePeriod=basePeriod,
                    convertFlaggedToNA=convertFlaggedToNA,
                    ...)

  return(myDataclim)
}

.myCoef <- function(lm){
  return(c(trend=unname(coef(lm)["y"]), pval=summary(lm)$coef["y", 4]))
}

.colApplyZoo <- function(X, FUN, ...){
  if (class(X) != "zoo"){
    print("Wrong input to .colApplyZoo.")
    return(NULL)
  }
  
  FUN <- match.fun(FUN)
  dl <- length(dim(X))
  if (!dl || dl != 2) 
      stop("In .colApplyZoo: dim(X) must have length 2")
  output <- sapply(1:ncol(X), function(j) FUN(X[, j], ...))
  if (is.null(dim(output))){
    names(output) <- names(X)
  } else {
    colnames(output) <- names(X)
  }
  
  return(output)
}

## show method
setMethod(
    f="show",
    signature="dataclim",
    definition=function(object){
      cat("*** Class dataclim, method Show *** \n")
      cat("* Slot data =\n")
      print(object@data[1:5, ])
      cat("[...] Truncated\n\n")
      
      cat("* Slot annualAvg =\n")
      print(object@annualAvg[1:5, ])
      cat("[...] Truncated\n\n")

      cat("* Slot monthlyAvg =\n")
      print(object@monthlyAvg[1:5, ])
      cat("[...] Truncated\n\n")

      cat("* Slot flagged =\n")
      print(lapply(object@flagged, function(l) l[1:(min(length(l), 5))]))

      cat("* Slot basePeriod =\n")
      print(object@basePeriod)

      cat("\n\n* Slot convertFlaggedToNA =")
      print(object@convertFlaggedToNA)

      cat("\n* Further slots: annualMax, annualMin, monthlyMax, monthlyMin\n")

      cat("\n* Access slots by slot(dataclim, slotname), see ?slot for help.\n")
      cat("\n******* End Show (dataclim) ******* \n")
}
)

## summary method
setMethod(
    f="summary",
    signature="dataclim",
    definition=function(object, basePeriod=NULL){
      if (is.null(basePeriod))
          basePeriod <- object@basePeriod

      ## overall statistics
      means <- .colApplyZoo(object@data[year(object@data) %in% basePeriod, ],
                     mean, na.rm=TRUE)

      sds <- .colApplyZoo(object@data[year(object@data) %in% basePeriod, ],
                     sd, na.rm=TRUE)

      ## annual cycle
      monthNames <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
                      "Sep", "Oct", "Nov", "Dec")
      acAvg <- .colApplyZoo(object@monthlyAvg[year(object@monthlyAvg) %in%
                                              basePeriod, ],
                            aggregate, month, mean, na.rm=TRUE)
      rownames(acAvg) <- monthNames

      acMin <- .colApplyZoo(object@monthlyMin[year(object@monthlyMin) %in%
                                              basePeriod, ],
                            aggregate, month, min, na.rm=TRUE)
      rownames(acMin) <- monthNames

      acMax <- .colApplyZoo(object@monthlyMax[year(object@monthlyMax) %in%
                                              basePeriod, ],
                            aggregate, month, max, na.rm=TRUE)
      rownames(acMax) <- monthNames

      ## trends of the averages
      annualTrends <-
          .colApplyZoo(object@annualAvg[time(object@annualAvg) %in%
                                        basePeriod, ],
                function(x) {.myCoef(lm(v~y, data=data.frame(v=coredata(x),
                                                 y=index(x))))})

      monthlyTrends <- list(trend=matrix(NA, ncol=ncol(object@monthlyAvg),
                                   nrow=12),
                               pval=matrix(NA, ncol=ncol(object@monthlyAvg),
                                   nrow=12))
      for (m in 1:12){
        mm <- .colApplyZoo(object@monthlyAvg[year(object@monthlyAvg) %in%
                                             basePeriod &
                                   month(object@monthlyAvg) == m, ],
                    function(x) {
                      .myCoef(lm(v~y,
                                data=data.frame(v=coredata(x), y=year(x))))})
        monthlyTrends$trend[m, ] <- mm["trend", ]
        monthlyTrends$pval[m, ] <- mm["pval", ]
        colnames(monthlyTrends$trend) <- colnames(monthlyTrends$pval) <-
            names(object@monthlyAvg)
        rownames(monthlyTrends$trend) <- rownames(monthlyTrends$pval) <-
            monthNames
      }
      ## return
      return(list(means=means, stdDevs=sds,
                  annualCycleAvg=acAvg, annualCycleMin=acMin,
                  annualCycleMax=acMax,
                  annualTrends=annualTrends,
                  monthlyTrends=monthlyTrends, basePeriod=basePeriod))
    }
    )

## homogeneity tests from knmi
SNHtest <- function(X){
  if (!is.numeric(X) || length(X) < 20 || sum(is.na(X)) > 0){
    print("Wrong input to SNHtest.")
    return(NULL)
  }
  
  s <- sd(X)
  counts <- 1:length(X)
  
  ## standard normal homogeneity test (SNH)
  z1 <- cumsum(X - mean(X)) / counts / s
  z2 <- (sum(X - mean(X)) - cumsum(X - mean(X))) /
      (rev(counts) - 1) / s
  ## test statistic
  T0 <- max(z1^2 * counts + z2^2 * (rev(counts) - 1), na.rm=TRUE)
  breakpoint <- which.max(z1^2 * counts + z2^2 * (rev(counts) - 1))
  ## compare to critical values
  crit.vals <- data.frame(n=c(20, 30, 40, 50, 70, 100),
                          p1=c(9.56, 10.45, 11.01, 11.38, 11.89, 12.32),
                          p5=c(6.95, 7.65, 8.10, 8.45, 8.80, 9.15))
  ## save result
  my.slot <- min(findInterval(length(X) - 0.1, crit.vals$n) + 1,
                 nrow(crit.vals))

  output <- list(statistic=T0, breakpoint=breakpoint,
              significance=(c(names(crit.vals)[2:3],
                              "NS")[which(c(T0 > crit.vals[my.slot, 2:3],
                                            TRUE))[1]]))
  return(output)
}

BHRtest <- function(X){
  if (!is.numeric(X) || length(X) < 20 || sum(is.na(X)) > 0){
    print("Wrong input to BHRtest.")
    return(NULL)
  }

  s <- sd(X)
  counts <- 1:length(X)
  
  S <- c(0, cumsum(X - mean(X)))
  ## test statistic
  R <- (max(S) - min(S)) / s / sqrt(length(X))
  breakpoint <- which.max(S)
  if (abs(max(S)) < abs(min(S)))
      breakpoint <- which.min(S)
  
  ## compare to critical values
  crit.vals <- data.frame(n=c(20, 30, 40, 50, 70, 100),
                          p1=c(1.60, 1.70, 1.74, 1.78, 1.81, 1.86),
                          p5=c(1.43, 1.50, 1.53, 1.55, 1.59, 1.62))
  ## save result
  my.slot <- min(findInterval(length(X) - 0.1, crit.vals$n) + 1,
                 nrow(crit.vals))

  output <- list(statistic=R, breakpoint=breakpoint,
              significance=(c(names(crit.vals)[2:3],
                              "NS")[which(c(R > crit.vals[my.slot, 2:3],
                                            TRUE))[1]]) )
  return(output)
}

PETtest <- function(X){
  if (!is.numeric(X) || length(X) < 20 || sum(is.na(X)) > 0){
    print("Wrong input to PETtest.")
    return(NULL)
  }

  counts <- 1:length(X)

  X <- 2 * cumsum(rank(X)) - counts * (length(X) + 1)
  ## test statistic and position of timestep of break
  Xe <- max(abs(X))
  breakpoint <- which.max(abs(X))
  
  ## compare to critical values
  crit.vals <- data.frame(n=c(20, 30, 40, 50, 70, 100),
                          p1=c(71, 133, 208, 293, 488, 841),
                          p5=c(57, 107, 167, 235, 393, 677))
  ## save result
  my.slot <- min(findInterval(length(X) - 0.1, crit.vals$n) + 1,
                 nrow(crit.vals))

  output <- list(statistic=Xe, breakpoint=breakpoint,
              significance=(c(names(crit.vals)[2:3],
                              "NS")[which(c(Xe > crit.vals[my.slot, 2:3],
                                            TRUE))[1]]) )
  return(output)
              
}

VONtest <- function(X){
  if (!is.numeric(X) || length(X) < 20 || sum(is.na(X)) > 0){
    print("Wrong input to VONtest.")
    return(NULL)
  }

  ## test statistic
  N <- sum((-diff(X))^2) / sum((X - mean(X))^2)
  ## compare to critical values
  crit.vals <- data.frame(n=c(20, 30, 40, 50, 70, 100),
                          p1=c(1.04, 1.20, 1.29, 1.36, 1.45, 1.54),
                          p5=c(1.30, 1.42, 1.49, 1.54, 1.61, 1.67))
  ## save result
  my.slot <- min(findInterval(length(X) - 0.1, crit.vals$n) + 1,
                 nrow(crit.vals))

  output <- list(statistic=N, breakpoint=NULL,
              significance=(c(names(crit.vals)[2:3],
                              "NS")[which(c(N < crit.vals[my.slot, 2:3],
                                            TRUE))[1]]))
  return(output)
}


setGeneric(name="evalHomogeneity",
           def=function(X){standardGeneric("evalHomogeneity")})
setMethod(
    f="evalHomogeneity",
    signature="dataclim",
    definition=function(X){
      
      ## test variables
      tv.daily <- data.frame(DTR=X@data[, "tmax"] - X@data[, "tmin"])
      tv.daily$vDTR <- c(NA, abs(diff(tv.daily$DTR)))
      tv.daily$RR1 <- X@data[, "prec"] > 1
      tv.daily <- zoo(tv.daily, time(X@data))
      
      tv.annual <- aggregate(tv.daily, year,
                             function(v){if (sum(!is.na(v)) > 330){
                               mean(v, na.rm=TRUE)
                             } else {
                               NA
                             }})
      
      ## treat missing values: at least 70% valid years are required
      ## and at least 20 valid years. gaps within this limit are
      ## filled (linear interpolation).
      for (v in names(tv.annual)){
        if (sum(!is.na(tv.annual[, v])) / nrow(tv.annual) >= .7 &&
            sum(!is.na(tv.annual)[, v]) >= 20){
          tv.annual[, v] <- na.fill(tv.annual[, v], "extend")
        } else {
          print(paste("Not enough valid data in", v))
          return(NULL)
        }
      }
      
      ## now evaluate these series
      breaks <- test.results <-
          data.frame(matrix(ncol=ncol(tv.annual), nrow=4))
      names(test.results) <- names(breaks) <- names(tv.annual)
      rownames(test.results) <- rownames(breaks) <-
          c("SNH", "BHR", "PET", "VON")
      
      ## loop over variables
      for (v in names(tv.annual)){
        X <- as.vector(tv.annual[, v])

        ## standard normal homogeneity test (SNH)
        myTest <- SNHtest(X)
        test.results["SNH", v] <- myTest$significance
        breaks["SNH", v] <- index(tv.annual)[myTest$breakpoint]

        ## Buishand range test (BHR)
        myTest <- BHRtest(X)
        test.results["BHR", v] <- myTest$significance
        breaks["BHR", v] <- index(tv.annual)[myTest$breakpoint]
        
        ## Pettitt test (PET)
        myTest <- PETtest(X)
        test.results["PET", v] <- myTest$significance
        breaks["PET", v] <- index(tv.annual)[myTest$breakpoint]
        
        ## Von Neumann ratio (VON)
        myTest <- VONtest(X)
        test.results["VON", v] <- myTest$significance
        breaks["VON", v] <- NA

      } ## loop over variables
      
      ## aggregate test results
      test.agg <- apply(test.results == "p1", 2, sum)
      temp.score <- max(test.agg[c("DTR", "vDTR")])
      prec.score <- test.agg["RR1"]
      classes <- c("useful", "doubtful", "suspect")
      temp.class <- classes[findInterval(temp.score, c(0, 1.5, 2.5))]
      prec.class <- classes[findInterval(prec.score, c(0, 1.5, 2.5))]
      
      return(list(tests=test.results, breaks=breaks,
                  classes=c(temp=temp.class,
                                          prec=prec.class)))
    }
    )

setMethod(
    f="evalHomogeneity",
    signature="data.frame",
    definition=function(X){
      if (nrow(X) < 20 || sum(is.na(X)) > 0){
        print("Data.frame for evalHomogeneity too small or contains NAs.")
        return(NULL)
      }
      
      ## now evaluate these series
      test.results <- breaks <- data.frame(matrix(ncol=ncol(X), nrow=4))
      names(test.results) <- names(breaks) <- names(X)
      rownames(test.results) <- rownames(breaks) <-
          c("SNH", "BHR", "PET", "VON")
      
      ## loop over variables
      for (v in names(X)){
        myX <- as.numeric(X[, v])

        ## standard normal homogeneity test (SNH)
        myTest <- SNHtest(myX)
        test.results["SNH", v] <- myTest$significance
        breaks["SNH", v] <- myTest$breakpoint

        ## Buishand range test (BHR)
        myTest <- BHRtest(myX)
        test.results["BHR", v] <- myTest$significance
        breaks["BHR", v] <- myTest$breakpoint
        
        ## Pettitt test (PET)
        myTest <- PETtest(myX)
        test.results["PET", v] <- myTest$significance
        breaks["PET", v] <- myTest$breakpoint
        
        ## Von Neumann ratio (VON)
        myTest <- VONtest(myX)
        test.results["VON", v] <- myTest$significance
        breaks["VON", v] <- NA
        
      } ## loop over variables
      
      ## aggregate test results
      test.agg <- apply(test.results == "p1", 2, sum)
      classes <- c("useful", "doubtful",
                   "suspect")[findInterval(test.agg, c(0, 1.5, 2.5))]
      names(classes) <- names(test.agg)

      return(list(tests=test.results, breaks=breaks,
                  classes=classes))
    }
    )


createClimdex <- function(myDataclim, basePeriod=NULL){
  if (!(class(myDataclim) == "dataclim" && (is.null(basePeriod) ||
                 (is.numeric(basePeriod) && length(basePeriod) >= 2)))){
    print("Wrong input to createClimdex.")
    return(NULL)
  }
  if (is.null(basePeriod))
      basePeriod <- myDataclim@basePeriod
  
  pcic.dates <- as.PCICt(as.POSIXlt(time(myDataclim@data)), cal="gregorian")
  my.pcic <- climdexInput.raw(tmax=myDataclim@data$tmax,
                              tmin=myDataclim@data$tmin,
                              prec=myDataclim@data$prec,
                              tmax.dates=pcic.dates,
                              tmin.dates=pcic.dates,
                              prec.dates=pcic.dates,
                              base.range=c(min(basePeriod),
                                  max(basePeriod)))
  return(my.pcic)
}

computeClimdex <- function(myClimdex){
  if (class(myClimdex) != "climdexInput"){
    print("Wrong input to computeClimdex.")
    return(NULL)
  }

  ##annual
  all.annual <- data.frame(
      cdd=climdex.cdd(myClimdex),
      cwd=climdex.cwd(myClimdex),
      fd=climdex.fd(myClimdex),
      id=climdex.id(myClimdex),
      r20mm=climdex.r20mm(myClimdex),
      r99ptot=climdex.r99ptot(myClimdex),
      rx1day=climdex.rx1day(myClimdex,freq="annual"),
      sdii=climdex.sdii(myClimdex),
      tn10p=climdex.tn10p(myClimdex,freq="annual"),
      tnn=climdex.tnn(myClimdex,freq="annual"),
      tr=climdex.tr(myClimdex),
      tx90p=climdex.tx90p(myClimdex,freq="annual"),
      txx=climdex.txx(myClimdex,freq="annual"),
      csdi=climdex.csdi(myClimdex),
      dtr=climdex.dtr(myClimdex,freq="annual"),
      gsl=climdex.gsl(myClimdex),
      prcptot=climdex.prcptot(myClimdex),
      r10mm=climdex.r10mm(myClimdex),
      r95ptot=climdex.r95ptot(myClimdex),
      rnnmm=climdex.rnnmm(myClimdex),
      rx5day=climdex.rx5day(myClimdex,freq="annual"),
      su=climdex.su(myClimdex),
      tn90p=climdex.tn90p(myClimdex,freq="annual"),
      tnx=climdex.tnx(myClimdex,freq="annual"),
      tx10p=climdex.tx10p(myClimdex,freq="annual"),
      txn=climdex.txn(myClimdex,freq="annual"),
      wsdi=climdex.wsdi(myClimdex))
  all.annual <- zoo(all.annual, as.numeric(rownames(all.annual)))
  
  ##monthly
  all.monthly <- data.frame(
      rx1day=climdex.rx1day(myClimdex,freq="monthly"),
      tn10p=climdex.tn10p(myClimdex,freq="monthly"),
      tnn=climdex.tnn(myClimdex,freq="monthly"),
      tx90p=climdex.tx90p(myClimdex,freq="monthly"),
      txx=climdex.txx(myClimdex,freq="monthly"),
      dtr=climdex.dtr(myClimdex,freq="monthly"),
      rx5day=climdex.rx5day(myClimdex,freq="monthly"),
      tn90p=climdex.tn90p(myClimdex,freq="monthly"),
      tnx=climdex.tnx(myClimdex,freq="monthly"),
      tx10p=climdex.tx10p(myClimdex,freq="monthly"),
      txn=climdex.txn(myClimdex,freq="monthly"))
  all.monthly <- zoo(all.monthly, as.yearmon(rownames(all.monthly)))
  
  return(list(annual.climdex=all.annual, monthly.climdex=all.monthly))
}



setGeneric(name="plotWalterLieth",
           def=function(X, station="", alt=NA, per="", margin=c(4, 4, 5, 4),
               pcol="#005ac8", tcol="#e81800", pfcol="#79e6e8",
               sfcol="#09a0d1",
               shem=FALSE, ...){standardGeneric("plotWalterLieth")})


setMethod(
    f="plotWalterLieth",
    signature="matrix",
    definition=function (X, station="", alt=NA, per="",
        margin=c(4, 4, 5, 4), pcol="#005ac8",
        tcol="#e81800", pfcol="#79e6e8", 
        sfcol="#09a0d1", shem=FALSE, p3line=FALSE, ...) {
  #This function is a customized version of the function provided in
  #The climatol-package, allowing for some extra annotation
  #http://www.inside-r.org/packages/cran/climatol/docs/diagwl
  # WALTER H & LIETH H (1960): Klimadiagramm Weltatlas. G. Fischer, Jena.
  #
  # Args:
  #   X: Monthly climatic data for which the diagram will be plotted.
  #   station: Name of the climatological station
  #   alt: Altitude of the climatological station
  #   per:  Period on which the averages have been computed
  #   margin: Margins vector for the plot (to be passed to par).
  #   mlab: Month labels for the X axis:
  #   pcol: Color pen for precipitation.
  #   tcol: Color pen for temperature.
  #   pfcol: Fill color for probable frosts.
  #   sfcol: Fill color for sure frosts.
  #   shem: Set to TRUE for southern hemisphere stations.
  #   p3line: Set to TRUE to draw a suplementary precipitation line referenced
  #           to three times the temperature (as suggested by Bogdan Rosca).
  #   ...: Other graphic parameters
  #
  # Returns:
  #   NA
  #
  ## 

    old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
    par(mar=margin, pty="s", las=1, new=FALSE)
    nr <- nrow(X)
    if (nrow(X) != 4 || sum(is.na(X)) > 0 || !class(X) == "matrix") 
        stop("Need a 4x12 matrix as input, no NAs!\n")

    mlab <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", 
              "O", "N", "D")
    X <- as.matrix(X)
    if (shem) {
        m1 <- X[, 1:6]
        m2 <- X[, 7:12]
        X <- cbind(m2, m1)
        mlab <- c(mlab[7:12], mlab[1:6])
    }
    p <- X[1, ]
    if (nr == 2) 
        tm <- X[2, ]
    else tm <- apply(X[2:3, ], 2, mean)
    pmax <- max(p)
    ymax <- 60
    if (pmax > 300) 
        ymax <- 50 + 10 * floor((pmax + 100)/200)
    ymin <- min(-1.5, min(tm))
    if (ymin < -1.5) {
        ymin = floor(ymin/10) * 10
        labT <- paste(ymin)
        labP <- ""
        if (ymin < -10) {
            for (i in (ymin/10 + 1):-1) {
                labT <- c(labT, i * 10)
                labP <- c(labP, "")
            }
        }
        labT <- c(labT, "0", "10", "20", "30", "40", "50", "")
        labP <- c(labP, "0", "20", "40", "60", "80", "100", "300")
    }
    else {
        labT <- c("0", "10", "20", "30", "40", "50", "")
        labP <- c("0", "20", "40", "60", "80", "100", "300")
    }
    if (ymax > 60) {
        for (i in 6:(ymax/10 - 1)) {
            labT <- c(labT, "")
            labP <- c(labP, 100 * (2 * i - 7))
        }
    }
    plot(0:13 - 0.5, c(tm[12], tm[1:12], tm[1]), xlim = c(0, 
        12), ylim = c(ymin, ymax), type = "n", xaxs = "i", yaxs = "i", 
        xaxp = c(0, 12, 12), xlab = "", ylab = "", xaxt = "n", 
        yaxt = "n", bty = "n")
    lmin <- ymin
    if (lmin == -1.5) 
        lmin = 0
    axis(2, ((lmin/10):(ymax/10)) * 10, labels = labT, col.axis = tcol)
    axis(4, ((lmin/10):(ymax/10)) * 10, labels = labP, col.axis = pcol)
    mtext("C", 2, col = tcol, las = 1, line = 3, adj = 0, at = 55)
    mtext("mm", 4, col = pcol, las = 1, line = 3, adj = 1, at = 55)
    abline(0, 0)
    abline(50, 0)
    if (is.na(alt)) 
        mtext(station, line = 2, adj = 0)
    else mtext(paste(station, " (", alt, " m)", sep = ""), line = 2, 
        adj = 0)
    mtext(per, line = 1, adj = 0)
    mtext(paste("Average Temp. ",round(mean(tm * 10))/10,
                "C\n Annual Precip. ", round(sum(p)), 
                " mm", sep = ""), line = 1, adj = 1)
    x <- 0:13 - 0.5
    p2 <- c(p[12], p[1:12], p[1])
    if (p3line) {
        yl3 <- c(p[12], p[1:12], p[1])/3
        yl3[yl3 > 50] <- 50
    }
    if (pmax <= 100) {
        xl <- x
        yl <- c(p[12], p[1:12], p[1])/2
        n2 <- 14
    }
    else {
        xp <- numeric(30)
        yp <- numeric(30)
        xl <- numeric(25)
        yl <- numeric(25)
        n <- 0
        n2 <- 0
        gr <- FALSE
        if (p2[1] > 100) {
            n <- n + 1
            xp[n] <- x[1]
            yp[n] <- 50
            n <- n + 1
            xp[n] <- x[1]
            yp[n] <- 50 + (p2[1] - 100)/20
            n2 <- n2 + 1
            xl[n2] <- x[1]
            yl[n2] <- 50
            gr <- TRUE
        }
        else {
            n2 <- n2 + 1
            xl[n2] <- x[1]
            yl[n2] <- p2[1]/2
        }
        for (i in 2:14) {
            if (gr) {
                n <- n + 1
                if (p2[i] > 100) {
                  xp[n] <- x[i]
                  yp[n] <- 50 + (p2[i] - 100)/20
                }
                else {
                  xp[n] <- x[i - 1] + (100 - p2[i - 1])/(p2[i] - 
                    p2[i - 1])
                  yp[n] <- 50
                  n2 <- n2 + 1
                  xl[n2] <- xp[n]
                  yl[n2] <- 50
                  n <- n + 1
                  xp[n] <- NA
                  yp[n] <- NA
                  n2 <- n2 + 1
                  xl[n2] <- x[i]
                  yl[n2] <- p2[i]/2
                  gr <- FALSE
                }
            }
            else {
                if (p2[i] > 100) {
                  n <- n + 1
                  xp[n] <- x[i - 1] + (100 - p2[i - 1])/(p2[i] - 
                    p2[i - 1])
                  yp[n] <- 50
                  if (xl[n2] != x[i - 1]) {
                    n2 <- n2 + 1
                    xl[n2] <- x[i - 1]
                    yl[n2] <- p2[i - 1]/2
                  }
                  n2 <- n2 + 1
                  xl[n2] <- xp[n]
                  yl[n2] <- 50
                  n <- n + 1
                  xp[n] <- x[i]
                  yp[n] <- 50 + (p2[i] - 100)/20
                  gr <- TRUE
                }
                else {
                  n2 <- n2 + 1
                  xl[n2] <- x[i]
                  yl[n2] <- p2[i]/2
                }
            }
        }
        if (!is.na(yp[n])) {
            n <- n + 1
            xp[n] <- xp[n - 1]
            yp[n] <- 50
            n2 <- n2 + 1
            xl[n2] <- 12.5
            yl[n2] <- 50
        }
        polygon(xp[1:n], yp[1:n], col = pcol, border = pcol)
    }
    pi <- approx(xl[1:n2], yl[1:n2], n = 66)$y
    ti <- approx(x, c(tm[12], tm[1:12], tm[1]), n = 66)$y
    ti[ti < 0] <- 0
    d <- pi - ti
    xi <- (1:66)/5 - 0.7
    xw <- subset(xi, d > 0)
    y1 <- subset(pi, d > 0)
    y2 <- subset(ti, d > 0)
    if (length(xw) > 0) 
        segments(xw, y1, xw, y2, col = pcol, lty = 1, lwd = 1)
    xw <- subset(xi, d < 0)
    y1 <- subset(pi, d < 0)
    y2 <- subset(ti, d < 0)
    if (length(xw) > 0) 
        segments(xw, y1, xw, y2, col = tcol, lty = 3, lwd = 2)
    for (i in 1:12) if (X[3, i] <= 0) 
        rect(i - 1, -1.5, i, 0, col = sfcol)
    for (i in 1:12) if (X[4, i] <= 0 && X[3, i] > 0) 
        rect(i - 1, -1.5, i, 0, col = pfcol)
    lines(xl[1:n2], yl[1:n2], col = pcol, lwd = 2)
    if (p3line) 
        lines(x, yl3)
    lines(x, c(tm[12], tm[1:12], tm[1]), col = tcol, lwd = 2)
    mtext(paste("Max\n",formatC(max(as.matrix(X[2, ])), digits = 1,
                                     format = "f")), 
        2, las = 1, line = 2, at = 35)
    mtext(paste("Min\n",formatC(min(as.matrix(X[3, ])), digits = 1,
                                format = "f")), 
        2, las = 1, line = 2, at = 15)
    for (i in 0:13) segments(i, 0, i, -1.5)
    mtext(mlab, 1, las = 1, line = 0.5, adj = 0.5, at = x[2:13])
    invisible()
  })

setMethod(
    f="plotWalterLieth",
    signature="dataclim",
    definition=function(X, station="", alt=NA, per="",
        margin=c(4, 4, 5, 4), pcol="#005ac8",
        tcol="#e81800", pfcol="#79e6e8", 
        sfcol="#09a0d1", shem=FALSE, p3line=FALSE, ...){
      summary <- summary(X)

      wlMatrix <- rbind(summary$annualCycleAvg[, "prec"],
                        summary$annualCycleAvg[, "tmax"],
                        summary$annualCycleAvg[, "tmin"],
                        summary$annualCycleMin[, "tmin"])
      per <- paste(min(summary$basePeriod), max(summary$basePeriod), sep=":")
      plotWalterLieth(wlMatrix, station=station, alt=alt, per=per,
                      margin=margin, pcol=pcol, tcol=tcol, pfcol=pfcol,
                      sfcol=sfcol, shem=shem, p3line=FALSE, ...)
    })
