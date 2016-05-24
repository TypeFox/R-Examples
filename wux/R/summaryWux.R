
# ----------------------------------------------------------------
# $Author: thm $
# $Date: 2015-04-07 16:09:44 +0200 (Tue, 07 Apr 2015) $
# $Rev: 339 $
# ----------------------------------------------------------------

summary.wux.df <- function(object, parms = c("perc.delta.precipitation_amount",
                             "delta.air_temperature"),
                           average.over.gcm.runs = FALSE, ...){
  ## Summary method for WUX data.frame (of class "wux.df"). Will be printed
  ## with method "print.summaryWuxdf".
  ##
  ## Args:
  ##   object: WUX data.frame obtained from models2wux()
  ##   parms: String (vector) specifying parameter to be evaluated
  ##   average.over.gcm.runs: Boolean. Should the runs of a GCM be averaged out? This is recommended!
  ##   '...': Of no use yet.
  ##
  ## Returns:
  ##   Object of class "summaryWux", which is a list containing statistcs.
  ## 
  ## structure of output list:
  ## overall.stats  ---\
  ##                    |-n.models
  ##                    |-n.gcms
  ##                    |-n.rcm
  ##                    |-scenarios --\
  ##                                  |-n.gcms
  ##                                  |-n.models.total
  ##                                  |-n.gcm.runs
  ##                                  |-n.rcms
  ##                                  |-n.rcm.runs
  ##                                  |-rcm.gcm.crosstable
  ## parms.stats ------\
  ##                    |-subregs --\
  ##                                 |-parms --\
  ##                                           |-em.scen --\
  ##                                                        |-season --\
  ##                                                                   |-mean
  ##                                                                   |-sd
  ##                                                                   |-coefvar
  ##                                                                   |-min 
  ##                                                                   |-max
  ##                                                                   |-med
  ##                                                                   |-q25
  ##                                                                   |-q75
  ## History:
  ##   2011-05-10 | First code (thm)
  ##   2011-10-13 | adding emission scenario analysis (thm)
  ##   2011-11-07 | dealing with missing values (thm)
  ##   2014-08-21 | now averaging over the runs of GCMs first, omitting missing values printout (thm)
  ##   2014-11-18 | minor bugfix
  ##   2014-11-20 | changed function nam from summary.wux to summaryWux, because I don't want this
  ##                to be a S3 method so far... we will introduce wux objects later (thm)
  ##
  ## TODO: possibility to average over RCM runs (should be easy to implement)


  ## extract some vectors from data.frame
  sea <- object$season
  em.scen <- object$em.scn
  subreg <- object$subreg
  gcm <- object$gcm
  inst <- object$institute
  model.names <- object$acronym
  em.scen.freq <- NULL

  ## initialize output vector
  out.list.entries <- c("overall.stats", "parms.stats")
  output.list <- vector("list", length(out.list.entries))
  names(output.list) <- out.list.entries
  
### OVERALL STATISTICS
  ## init list
  ## ov.stat.entries <- c("n.models", "n.gcms", "gcm.freq", "n.rcms", "rcm.freq",
  ##                      "rcm.gcm.crosstable", "em.scen.freq")
  ov.stat.entries <- c("n.models", "n.gcms", "n.rcms", "em.scen")##levels(em.scen)
  overall.stats <- vector("list", length(ov.stat.entries))
  names(overall.stats) <- ov.stat.entries
  
  em.scen.list <- vector("list", length(levels(factor(em.scen))))
  names(em.scen.list) <- levels(factor(em.scen))

  ov.stat.scenario.entries <- c("n.gcms", "n.models.total"," n.gcm.runs",
                                "n.rcms", "n.rcm.runs", "rcm.gcm.crosstable")
  overall.scenario.stats <- vector("list", length(ov.stat.scenario.entries))
  names(overall.scenario.stats) <- ov.stat.scenario.entries

  overall.stats$em.scen <- lapply(em.scen.list,
                                  function(x) `<-`(x, overall.scenario.stats))
   
  len <- tapply(model.names, list(sea, em.scen, subreg), length)
  ## in all subregions and seasons number of models the same?
  if (all(len[,1,] == len[,1,1][1]))
    overall.stats$n.models <- sum(len[1,,1]) #sum over emission scenarios
  else
    stop("NUMBER OF MODELS FOR SUBREGION-SEASON COMBINATION DIFFER. ADAPT THE SUMMARY FUNCTION (CURRENTLY ITS ASSUMED THAT THE NUMBER OF MODELS IS CONSTANT BETWEEN SEASONS AND SUBREGIONS)") 
  
  ind <- tapply(model.names, list(sea, em.scen, subreg))
  single.gcms <- object[which(ind == 1), ]$gcm
  if (all(is.na(object$rcm))){
    single.rcms <- NA
    overall.stats[["n.rcm.runs"]] <- NA
    overall.stats[["n.rcms"]] <- 0
    overall.stats[["rcm.gcm.crosstable"]] <- NA
  } else {
    single.rcms <- object[which(ind == 1), ]$rcm
    overall.stats[["n.rcm.runs"]] <- summary(single.rcms)
    overall.stats[["n.rcms"]] <- length(summary(single.rcms))
  ##  overall.stats[["rcm.gcm.crosstable"]] <- table(RCM = single.rcms, GCM = single.gcms)
  }
  gcm.freq <- summary(single.gcms)
  overall.stats[["n.gcms"]] <- length(gcm.freq[gcm.freq != 0])

### SCENARIO DEPENDEND INFORMATION
  em.scen.levels <- levels(factor(em.scen))
  ## for each scenario get relevant overall statistics...
  for (scen in levels(factor(em.scen))) {
    object.scen <- object[object[["em.scn"]] == scen, ]

    ind.scen <- tapply(object.scen$acronym,
                       list(object.scen$season, object.scen$subreg))
    single.gcms.scen <- factor(object.scen[which(ind.scen == 1), ]$gcm)

    gcm.freq.scen <- summary(single.gcms.scen)
    gcm.freq.scen <-  gcm.freq.scen[gcm.freq.scen != 0]
    overall.stats$em.scen[[scen]][["n.gcm.runs"]] <- gcm.freq.scen
    overall.stats$em.scen[[scen]][["n.gcms"]] <- length(gcm.freq.scen)
    overall.stats$em.scen[[scen]][["n.models.total"]] <- sum(gcm.freq.scen)
  
    if (all(is.na(object.scen$rcm))){
      single.rcms <- NA
      overall.stats$em.scen[[scen]][["n.rcm.runs"]] <- NA
      overall.stats$em.scen[[scen]][["n.rcms"]] <- NA
      overall.stats$em.scen[[scen]][["rcm.gcm.crosstable"]] <- NA
    } else {
      single.rcms <- object.scen[which(ind.scen == 1), ]$rcm
      overall.stats$em.scen[[scen]][["n.rcm.runs"]] <- summary(single.rcms)
      overall.stats$em.scen[[scen]][["n.rcms"]] <- length(summary(single.rcms))
      ## overall.stats$em.scen[[scen]][["rcm.gcm.crosstable"]] <-
      ##   table(RCM = single.rcms, GCM = single.gcms)
    }
    em.scen.freq <- c(em.scen.freq, sum(gcm.freq.scen))
    names(em.scen.freq)[[length(em.scen.freq)]] <- scen
  }

  overall.stats[["em.scen.freq"]] <- em.scen.freq[em.scen.freq != 0]
  ## append overall stats list to output list
  output.list[["overall.stats"]] <- overall.stats
  
### GET AGGREGATED STATISTICS FOR ALL SEASON/SUBREG/PARM COMBINATIONS
  list.allparms <- vector("list", length(parms))
  names(list.allparms) <- parms

  for (index.parm in seq(along = parms)){
    parm.name <- parms[index.parm]
    parm.val <- object[[parm.name]]
    if (is.null(parm.val)) {
      cat("PARAMETER", parm.name, "NOT IN WUX DATAFRAME\n")
      stop()
    }

    ## should we aggregate over GCM runs? else buisness as usual
    if (average.over.gcm.runs) {
      ## aggregate over same GCMs (average out the runs)
      gcms.avg <- tapply(parm.val, list(sea, em.scen, subreg,gcm),
                         function(x) mean(x, na.rm = TRUE))
      ## statistics for each subreg/season
      n <- apply(gcms.avg, c(1,2,3), function(x)  sum(!is.na(x)))
      ## those missing values make no sense, because NAs are produced by previous tapply averaging
      ## anyway, so I deleted it (thm)
      ## n.na <- apply(gcms.avg, c(1,2,3),
      ##                function(x) sum(is.na(x)))
      mea <- apply(gcms.avg, c(1,2,3), mean, na.rm = TRUE)
      sd <- apply(gcms.avg, c(1,2,3), sd, na.rm = TRUE)
      coef.var <- sd / abs(mea)
      min <- apply(gcms.avg, c(1,2,3), min, na.rm = TRUE)
      max <- apply(gcms.avg, c(1,2,3), max, na.rm = TRUE)
      med <- apply(gcms.avg, c(1,2,3), median, na.rm = TRUE)
      q25 <- apply(gcms.avg, c(1,2,3), quantile, probs = 0.25,
                   na.rm = TRUE)
      q75 <- apply(gcms.avg, c(1,2,3), quantile, probs = 0.75,
                   na.rm = TRUE)
    } else {
      ## statistics for each subreg/season (old version from 2011)
      n <- tapply(parm.val, list(sea, em.scen, subreg),
                  function(x) sum(!is.na(x)))
      n.na <- tapply(parm.val, list(sea, em.scen, subreg),
                     function(x) sum(is.na(x)))
      mea <- tapply(parm.val, list(sea, em.scen, subreg), mean, na.rm = TRUE)
      sd <- tapply(parm.val, list(sea, em.scen, subreg), sd, na.rm = TRUE)
      coef.var <- sd / abs(mea)
      min <- tapply(parm.val, list(sea, em.scen, subreg), min, na.rm = TRUE)
      max <- tapply(parm.val, list(sea, em.scen, subreg), max, na.rm = TRUE)
      med <- tapply(parm.val, list(sea, em.scen, subreg), median, na.rm = TRUE)
      q25 <- tapply(parm.val, list(sea, em.scen, subreg), quantile, probs = 0.25,
                    na.rm = TRUE)
      q75 <- tapply(parm.val, list(sea, em.scen, subreg), quantile, probs = 0.75,
                    na.rm = TRUE)
    }
    list.allparms[[index.parm]] <- list(n = n,
                                        ## missing = n.na,
                                        mean = mea,
                                        sd = sd,
                                        coefvar = coef.var,
                                        min = min, 
                                        max = max,
                                        med = med,
                                        q25 = q25,
                                        q75 = q75
                                        )
  }
 
### generate nested output list (subreg-parms-scenario-season)
  lsub <- length(levels(subreg))
  lsea <- length(levels(sea))
  lemscen <- length(levels(em.scen))

  ## initialize (sub-)lists for metric output (i.e. parameter statistics)
  olist.metric <- vector("list", lsub)
  names(olist.metric) <- levels(subreg)
  ## list of emission scenarios
  list.scenarios <- vector("list", lemscen)
  names(list.scenarios) <- levels(em.scen)
  ## list at seasonal level
  list.seasons <- vector("list", lsea)
  names(list.seasons) <- levels(sea)
  ## list at parameter level
  list.parms <-  vector("list", length(parms))
  names(list.parms) <- parms

  ## extract statistics for subreg-parms-season
  for (i in seq.int(levels(subreg))){
    for (parm.ind in seq(along = parms)){ ##(par in parms){
      ## current parameter name
      par <- parms[parm.ind]
      for (j in seq.int(levels(em.scen))){
        
        for (k in seq.int(levels(sea))){
          ## get value for season j
          list.seasons[[k]] <- sapply(list.allparms[[parm.ind]], '[', k, j, i)
        }
        ## append season list to parm list
        list.scenarios[[j]] <- list.seasons
      }
      list.parms[[parm.ind]] <- list.scenarios    
    }
    ## append parm-season lipst to subreg list
    olist.metric[[i]] <- list.parms
  }
  
  ## append parameter stats list to output list
  output.list[["parms.stats"]] <- olist.metric

  
###  returns olist.metric
  ## assign object to class "summaryWuxdf"
  class(output.list) <- c("summaryWuxdf", "list")
  return(output.list)
}



print.summaryWuxdf <- function(x, ...){
  ## Prints (via "cat") an "summaryWux" object in a nice way
  ##
  ## Args:
  ##   x: "summaryWux" object
  ##   '...': Nothing yet
  ##
  ## Returns:
  ##   Nothing, just printing.
  ## 
  ## History:
  ##   2011-05-10 | First code (thm)
  ##   2011-10-20 | Adapted to scenario cases (thm)

### print overall statistics
  cat("    ----------------------------------------------------------------------\n")
  cat("    ------------------------ OVERALL FREQUENCIES -------------------------\n")
  cat("    ----------------------------------------------------------------------\n")
  overall <- x$overall.stats
  cat("total number of climate models, including runs:",
      overall[["n.models"]], "\n")
  cat("total number of climate models by emission scenarios:\n")
  print(overall[["em.scen.freq"]])

  cat("number of GCMs (independent from scenario):", overall[["n.gcms"]], "\n")
  cat("number of RCMs (independent from scenario or institution):",
      overall[["n.rcms"]], "\n")
  cat("\n")
  ## omit rcm print if there are no RCMs
  if (overall[["n.rcms"]] != 0){
    cat("Number of RCMs used:\n")
    print(overall[["n.rcm.runs"]])
    cat("\n")
    cat("RCM-GCM matrix for all scenarios:\n")
    print(overall[["rcm.gcm.crosstable"]])
    cat("\n")
  }

### overall information by scenario
  cat("    ----------------------------------------------------------------------\n")
  cat("    ----------------------- FREQUENZIES BY SCENARIO ----------------------\n")
  cat("    ----------------------------------------------------------------------\n")
  for (scen in names(overall$em.scen)){
    cat(scen, ":\n", sep="")
    overall.scen <- overall$em.scen[[scen]]
    cat(" ", overall.scen$n.gcms, "GCMs (disregarding runs)\n")
    cat(" ", overall.scen$n.models.total, "models total\n")
    cat("  Number of GCMs used:\n")
    print(overall.scen$n.gcm.runs)
    if (all(!is.na(overall.scen$n.rcm.runs))){
      cat("  Number of RCM runs:\n")
      print(overall.scen$n.rcm.runs)
    }
   if (!is.null(overall.scen$n.rcms)){
     if (!is.na(overall.scen$n.rcms)){
      cat("  Number of RCMs:", overall.scen$n.rcms, "\n")
    }}
    if(all(!is.null(overall.scen$rcm.gcm.crosstable))){
      if(all(!is.na(overall.scen$rcm.gcm.crosstable))){
      cat("  GCM-RCM crosstable:\n")
      print(overall.scen$rcm.gcm.crosstable)
    }}
  }
  
### print statistics for metric parameters
  cat("\n")
  cat("    ----------------------------------------------------------------------\n")
  cat("    ---------------- CLIMATE MODEL STATISTICS BY SUBREGION ---------------\n")
  cat("    ----------------------------------------------------------------------\n")
  metric.list <- x$parms.stats
  ## for each subregion...
  for (i in seq(along = metric.list)){
    subreg.list <- metric.list[[i]]
    cat("\n------------", names(metric.list)[i], "------------\n")
    ## for each parameter
    for (parm in names(subreg.list)){
      parm.list <- subreg.list[[parm]]
      cat(parm, ": \n", sep = "")
      ## for each scenario....
      for (j in seq(along = parm.list)){
        scenario.list <- parm.list[[j]]
        scenario <- names(parm.list)[j]
        cat(" [", scenario,"]\n", sep = "")
        cat("     ", names(parm.list[[1]][[1]]), "\n", sep="\t")
        ## for each season....
        for (k in seq(along = scenario.list)){
          season <- names(scenario.list)[k]
          cat("   ")
          season.stats <- scenario.list[[k]]
          season.stats <- round(unlist(season.stats), 2)
          cat(season, ":\t", sep="")
          cat(season.stats,"\n", sep="\t")
        }
        cat("\n")
     }
    }
  }
  ## returns invisible object
  invisible(x)
}

