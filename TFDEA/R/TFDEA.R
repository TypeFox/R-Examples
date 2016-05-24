#******************************************************************************
#
# Copyright Tom Shott, ETA, PSU, 2013. 2014
# Use granted under BSD license terms
#
# R TFDEA Package
#
# TFDEA function
#
# NOTE: All efficiencies in program are in 0 - 1 form, output orientatation values are
# converted to 1/eff on calculation.
#
# Note: Man page - comment that default orientation is output, different from DEA
#
# debug option controls how much output the code generates. 0 - none, 3 is a lot
#
#******************************************************************************
#

TFDEA <- function (x, y, dmu_date_rel, date_forecast, rts="vrs", orientation="output",
                   second="min", mode="static", segroc=FALSE, debug=1){

  rts         <- .checkOption(rts,         "rts",          options.rts.l)
  orientation <- .checkOption(orientation, "orientation",  options.orientation.l)
  second      <- .checkOption(second,      "second",       options.second.l)
  mode        <- .checkOption(mode,        "mode",         options.mode.l)
  segroc      <- .checkOption(segroc,      "segroc",       TRUE)
  debug       <- .checkOption(debug,       "debug",        0)

  #
  # Check Input Values for correctness and consistency in size
  #
  x <- .checkData(x, "x")
  y <- .checkData(y, "y")
  if (nrow(x) != nrow(y))
    stop("Number of DMU's in inputs != number of DMU's in outputs")

  .checkDataGood(x, y, debug=debug)

  dmu_date_rel <- .checkVector(dmu_date_rel, "dmu_date_rel")
  if (length(dmu_date_rel) != nrow(x))
    stop("Number of DMU's in date intro != number of DMU's in inputs", call. = FALSE)
  dmu_date_rel <- array(dmu_date_rel, c(nrow(x)))         # Make same orientation other array

  .checkOption(date_forecast, "date_forecast", 0)
  if(date_forecast < min(dmu_date_rel) || date_forecast > max(dmu_date_rel))
    stop("Forecast date must be between min and max intro_date", call. = FALSE)

  # Unique sorted dates for calculations, as.vector to get right orientation
  date.soa.l      <- as.vector(sort(unique(dmu_date_rel[dmu_date_rel <= date_forecast])))
  if (debug >= 2) cat("Date Cur Forecast:", date_forecast, "Unique SOA Dates: ",
                      paste0(date.soa.l, collapse = ", "), "\n")


  #<New Page>
  ###############################################################################################
  #
  # Phase 1 - for each technology SOA set determine which DMU's are efficient at intro
  #
  ###############################################################################################
  nd                <- nrow(x)              # number of units, firms, DMUs
  dmu.names         <- rownames(x)

  # Values At Release (REL) date product first introduced
  dmu.eff.rel       <- array(NA, c(nd),     list(dmu=dmu.names))
  dmu.lambda.rel    <- array(NA, c(nd,nd),  list(dmu=dmu.names, dmu2=dmu.names))

  # Current Values (CUR) - values for largest date less then or equal forecast date
  dmu.eff.cur       <- array(NA, c(nd),     list(dmu=dmu.names))
  dmu.lambda.cur    <- array(NA, c(nd,nd),  list(dmu=dmu.names, dmu2=dmu.names))
  dmu.roc.cur       <- array(NA, c(nd),     list(dmu=dmu.names))
  dmu.sroc.cur      <- array(NA, c(nd),     list(dmu=dmu.names))
  dmu.date.cur      <- array(date_forecast, c(nd),  list(dmu=dmu.names))

  # Loop for each unique date less then or equal forecast date
  for(t in date.soa.l){
    dmu.eff.cur.b   <- isStdEfficient(dmu.eff.cur)        # DMU's from earlier years still eff
    dmu.cur.b       <- (dmu_date_rel == t)                # DMU's from current year
    dmu.soa.b       <- dmu.eff.cur.b | dmu.cur.b          # DMU's to use for eff calc

    if (debug >= 3) {
      cat("\nEvaluate SOA for date=", t, "\n")
      cat("Eff Cur DMUs:", dmu.names[dmu.eff.cur.b],"\n")
      cat("New Cur DMUs:", dmu.names[dmu.cur.b], "\n")
      cat("SOA Cur DMUs:", dmu.names[dmu.soa.b], "\n")
    }

    # Calculate Eff for All DMU's in SOA for each time period
    # WARNING: Use stdeff option to make all eff 0 - 1, stdeff, even for output orientation
    # index options select subset of DMU's
    results <- .dea(x, y, rts, orientation, second=second, z=dmu_date_rel,
                    slack=FALSE, stdeff=TRUE, index.K=dmu.soa.b, index.T=dmu.soa.b)
    dmu.eff.cur     <- results$eff
    dmu.lambda.cur  <- results$lambda

    # Save release eff & lambda for new DMU's released this year
    dmu.eff.rel[dmu.cur.b] <- dmu.eff.cur[dmu.cur.b]
    dmu.lambda.rel[dmu.cur.b, dmu.cur.b] <- dmu.lambda.cur[dmu.cur.b, dmu.cur.b]
  }

  table <- cbind(dmu_date_rel, dmu.eff.rel, dmu.eff.cur)
  colnames(table) <- c("Date", "Eff_Rel", "Eff_Cur")
  if (debug >= 2) {
    print("done Phase 1")
    print(table, digits=7)
  }

  #<New Page>
  ###############################################################################################
  #
  # Phase 2 - Calculate technology rate of change (ROC) for DMU's <= forecast date
  #
  # Calculated and saved eff, lambda in Phase 1
  #
  ###############################################################################################
  # Calc eff for all DMU's <= forecast date, were eff at release, must redo since need all DMU's
  dmu.cur.b <- (dmu_date_rel <= date_forecast) & isStdEfficient(dmu.eff.rel)
  if (debug >= 3) cat("Phase 2 - DMU Forecast Sample Set:", dmu.names[dmu.soa.b], "\n")

  # Calculate Eff for All DMU's in SOA
  results <- .dea(x, y, rts, orientation, second=second, z=dmu_date_rel,
                  slack=FALSE, stdeff=TRUE, index.K=dmu.cur.b, index.T=dmu.cur.b)
  dmu.eff.cur     <- results$eff
  dmu.lambda.cur  <- results$lambda

  if(mode == "dynamic")
    for (k in which(dmu.cur.b))
      dmu.date.cur[k] <- .wMean(dmu_date_rel, dmu.lambda.cur[k,])

  # Calculate ROC values for all DMU's eff at release, not efficient at cur,
  # have dmu.date.cur < date_forecast
  dmu.roc.b  <- (dmu_date_rel < date_forecast) & (dmu.date.cur > dmu_date_rel) &
    isStdEfficient(dmu.eff.rel) & !isStdEfficient(dmu.eff.cur)
  if (debug >= 3) cat("Phase 2 - DMU Forecast ROC Calc Set:", dmu.names[dmu.roc.b], "\n")

  dmu.roc.cur[dmu.roc.b] <- (dmu.eff.rel[dmu.roc.b] / dmu.eff.cur[dmu.roc.b]) ^
    ( 1 / (dmu.date.cur[dmu.roc.b] - dmu_date_rel[dmu.roc.b]) )

  for (k in which(dmu.cur.b & dmu.date.cur < dmu_date_rel)){
    cat("TFDEA Phase 2: DMU k=", k, " effective current ",
        "date < release date due to dynamic ROC and is being dropped\n")
    dmu.roc.cur[k] = NA
  }

  for (k in which(dmu.roc.cur > 10)){
    cat("TFDEA Phase 2: DMU k=", k, " has a numerically unstable ROC and is being dropped\n")
    dmu.roc.cur[k] = NA
  }

  average_roc <- mean(dmu.roc.cur, na.rm=TRUE)

  #
  # Segmented ROC Calc values
  #
  dmu.soa.b <- isStdEfficient(dmu.eff.cur)              # Forecast from DMU eff now (CUR)
  if (debug >= 3) cat("Forecast using SOA DMUs:", dmu.names[dmu.soa.b],"\n")

  dmu.sroc.cur[dmu.soa.b] <- average_roc                # Default ROC is average
  if(segroc){
    for(k in which(dmu.soa.b)){
      w.mean.roc <- .wMean(dmu.roc.cur, dmu.lambda.cur[,k])
      if (is.finite(w.mean.roc)){
        dmu.sroc.cur[k] <- w.mean.roc
      }
    }
  }

  table <- cbind(dmu_date_rel, dmu.eff.rel, dmu.eff.cur, dmu.date.cur, dmu.roc.cur, dmu.sroc.cur)
  colnames(table) <- c("Date", "Eff_Rel", "Eff_Cur", "EDate", "ROC", "S Roc")
  if (debug >= 2) {
    print(c("done Phase 2", "Avg ROC=", average_roc), digits=3)
    print(table, digits=7)
  }

  #<New Page>
  ###############################################################################################
  #
  # Phase 3 - Forecast (FOR) intro date
  #
  # Go through DMU's > forecast, calc forecast date based upon average ROC,
  # calculated stdeff (super efficiency) and cur year
  #
  ###############################################################################################

  # Values at forecast date (FOR) - date in future we are forecasting
  dmu.eff.for       <- array(NA, c(nd),     list(dmu=dmu.names))
  dmu.lambda.for    <- array(NA, c(nd,nd),  list(dmu=dmu.names, dmu2=dmu.names))
  dmu.sroc.for      <- array(NA, c(nd),     list(dmu=dmu.names))
  dmu.date.for      <- array(NA, c(nd),     list(dmu=dmu.names))

  # Calc DMU's to forecast
  dmu.for.b <- (dmu_date_rel > date_forecast)

  # Check that Average ROC is valid, if not valid ROC skip section. Likely no DMU's to use
  if (is.finite(average_roc) && sum(dmu.for.b) > 0){

    # Calc Super Effeciency
    results <- .sdea(x, y, rts=rts, orientation=orientation, second=second, z=dmu_date_rel,
                     stdeff=TRUE, slack=FALSE, index.T=dmu.soa.b, index.K=dmu.for.b)
    dmu.eff.for     <- results$eff
    dmu.lambda.for  <- results$lambda

    for (k in which(dmu.for.b)){

      if (!is.finite(dmu.eff.for[k])) next

      if (dmu.eff.for[k] <= 1){
        cat("TFDEA Phase 3: DMU k=", k, " is not super efficient at forecast and will"
            ,"not be forecast\n")
        next
      }

      if(mode == "dynamic")
        dmu.date.cur[k] <- .wMean(dmu_date_rel, dmu.lambda.for[k,])

      # Calculate Forecast from ROC, super efficiency & cur date
      peers.b <- (is.finite(dmu.lambda.for[k,]) & (dmu.lambda.for[k,] > 0) & (dmu.sroc.cur > 0))
      dmu.sroc.for[k] <- .wMean(dmu.sroc.cur[peers.b], dmu.lambda.for[k, peers.b])
      dmu.date.for[k] <- dmu.date.cur[k] +
        log(dmu.eff.for[k], exp(1)) / log(dmu.sroc.for[k], exp(1))
    }
  }

  table <- cbind(dmu_date_rel, dmu.eff.rel, dmu.eff.cur, dmu.date.cur,
                 dmu.roc.cur, dmu.sroc.cur, dmu.eff.for, dmu.sroc.for, dmu.date.for)
  colnames(table) <- c("Date", "Eff_Rel", "Eff_Cur", "EDate",
                       "ROC", "S Roc", "Eff_For", " S RocF", "Date For")
  if (debug >= 2){
    print(c("done Phase 3", "Avg ROC=", average_roc), digits=3)
    print(table, digits=7)
  }

  #
  # Return results
  # ToDo - Add option to include ux, vy, w
  results <- list(date_soa=date.soa.l,
                  dmu_eff_rel=dmu.eff.rel, dmu_lambda_rel=dmu.lambda.rel,
                  dmu_eff_cur=dmu.eff.cur, dmu_roc=dmu.roc.cur, dmu_lambda_cur=dmu.lambda.cur,
                  dmu_date_cur=dmu.date.cur, dmu_sroc_cur=dmu.sroc.cur,
                  dmu_eff_for=dmu.eff.for, dmu_lambda_for=dmu.lambda.for,
                  dmu_date_for=dmu.date.for, dmu_sroc_for=dmu.sroc.for,
                  roc=average_roc)
  return(results)
}

#
# Weighted mean
#
.wMean <-        function(value, weights){
  valid.b <- is.finite(value) & is.finite(weights)
  result  <- weighted.mean(value[valid.b], weights[valid.b], na.rm=TRUE)

  return(result)
}
