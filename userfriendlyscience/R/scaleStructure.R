scaleStructure <- scaleReliability <- function (dat=NULL, items = 'all', digits = 2,
                                                ci = TRUE, interval.type="normal-theory",
                                                conf.level=.95, silent=FALSE,
                                                samples=1000, bootstrapSeed = NULL,
                                                omega.psych = TRUE,
                                                poly = TRUE) {

  ### Make object to store results
  res <- list(input = as.list(environment()),
              intermediate = list(),
              output = list());
  
  ### If no dataframe was specified, load it from an SPSS file
  if (is.null(dat)) {
    dat <- getData(errorMessage=paste0("No dataframe specified, and no valid datafile selected in ",
                                       "the dialog I then showed to allow selection of a dataset.",
                                       "Original error:\n\n[defaultErrorMessage]"),
                   use.value.labels=FALSE);
    res$input$dat.name <- paste0("SPSS file imported from ", attr(dat, "filename"));
  }
  else {
    if (!is.data.frame(dat)) {
      stop("Argument 'dat' must be a dataframe or NULL! Class of ",
           "provided argument: ", class(dat));
    }
    res$input$dat.name <- deparse(substitute(dat));
  }
  
  ### if items contains only 1 element (or less), we
  ### include all items.
  if (length(items) <= 1) {
    ### Remove all cases with missing data (listwise deletion)
    res$input$dat <- na.omit(dat);
  }
  else {
    ### Select relevant items and remove all cases with
    ### missing data (listwise deletion)
    res$input$dat <- na.omit(subset(dat, select=items));
  }
  
  ### If any of the variables in the dataframe are factors, convert
  ### them to numeric vectors:
  if ("factor" %in% unlist(lapply(res$input$dat, 'class'))) {
    res$input$dat <- data.frame(lapply(res$input$dat, 'as.numeric'));
  }

  ### Set number of items and number of observations
  res$input$n.items <- ncol(res$input$dat);
  res$input$n.observations <- nrow(res$input$dat);
  
  if (samples < res$input$n.observations) {
    res$intermediate$samples <- samples <- 1.2*res$input$n.observations;
  }

  ### Also generate a dataframe (useful when
  ### requesting measures for many scales)
  res$output$dat <- data.frame(n.items        = res$input$n.items,
                               n.observations = res$input$n.observations);
                                     
  ### Get correlation matrix
  res$intermediate$cor <- cor(res$input$dat, use="complete.obs");
  
  ### Store number and proportion of positive correlations
  res$intermediate$cor.pos <-
    sum(res$intermediate$cor[lower.tri(res$intermediate$cor)] > 0);
  res$intermediate$cor.total <-
    (res$input$n.items^2 - res$input$n.items) / 2;
  res$intermediate$cor.proPos <-
    res$intermediate$cor.pos / res$intermediate$cor.total;
  
  ### Cronbach's alpha
  res$intermediate$alpha <- alpha(res$input$dat, check.keys=FALSE);
  res$output$cronbach.alpha <- res$intermediate$alpha$total$raw_alpha;
  res$output$dat$cronbach.alpha <- res$output$cronbach.alpha;
  
  ### GLB and Omega can only be computed if the number
  ### of items exceeds two
  
  if (res$input$n.items == 2) {
    ### Otherwise, compute Spearman Brown coefficient
    ### (see Eisinga, te Grotenhuis & Pelzer (2013). The reliability
    ### of a two-item scale: Pearson, Cronbach, or Spearman-Brown?
    ### International journal of public health, 58(4), 637-42.
    ### doi:10.1007/s00038-012-0416-3)
    ### Get r in numeric variable for convenience
    r <- res$intermediate$cor[1,2];
    res$intermediate$spearman.brown <- 1 / (1 + (1 / ((r/(1-r)) + (r/(1-r)))));
    res$output$spearman.brown <- res$intermediate$spearman.brown;
    res$output$dat$spearman.brown <- res$intermediate$spearman.brown;
  }
  else if (res$input$n.items > 2) {
    ### GLB
    suppressWarnings(res$intermediate$glb <- glb(res$input$dat));
    res$output$glb  <- res$intermediate$glb$glb.max;
    res$output$dat$glb  <- res$intermediate$glb$glb.max;
    ### Omega
    res$intermediate$omega <- ci.reliability(res$input$dat, type="omega");
    res$output$omega <- res$intermediate$omega$est;
    res$output$dat$omega <- res$intermediate$omega$est;
    suppressWarnings(res$intermediate$omega.psych <- omega(res$input$dat, plot=FALSE));
    res$output$omega.psych <- res$intermediate$omega.psych$omega.tot;
    res$output$dat$omega.psych.tot <- res$output$omega.psych;
    res$output$dat$omega.psych.h <- res$intermediate$omega.psych$omega_h;

    ### Examine largest number of levels for the items
    res$intermediate$maxLevels <- max(apply(res$input$dat, 2,
                                            function(x){return(nlevels(factor(x)))}));
    res$intermediate$maxRange <- max(res$input$dat, na.rm=TRUE) -
      min(res$input$dat, na.rm=TRUE) + 1;
    
    ### If sufficiently low, also provide ordinal estimates
    if (poly && res$intermediate$maxLevels < 9 && res$intermediate$maxRange < 9) {
      ### Compute polychoric correlation matrix
      res$intermediate$polychoric <-
        suppressWarnings(polychoric(res$input$dat)$rho);
      
      if (!any(is.na(res$intermediate$polychoric))) {
        ### Ordinal omega
        suppressWarnings(res$intermediate$omega.ordinal <-
                           omega(res$input$dat, poly=TRUE, plot=FALSE));
        ### Ordinal alpha
        res$intermediate$alpha.ordinal <- alpha(res$intermediate$polychoric,
                                                check.keys=FALSE);
      }
    }
    
    if (ci) {
      
      if (interval.type %in% c("bi", "perc", "bca")) {
        ### Use bootstrapping for the confidence intervals
        
        ### Check whether seed value is specified
        if (is.null(bootstrapSeed)) {
          bootstrapSeed <- as.numeric(format(Sys.time(), "%Y%m%d"));
          warning("No 'bootstrapSeed' specified. Without a seed for the ",
                  "random number generator, bootstrapping results will not  ",
                  "be exactly replicable. Setting current date as ",
                  "seed (", bootstrapSeed, "). Specify this value using the ",
                  "'bootstrapSeed' parameter to exactly replicate the ",
                  "bootstrapping results.");
        }
        
        ### Set seed
        res$input$bootstrapSeed <- bootstrapSeed;
        set.seed(res$input$bootstrapSeed);
        
        if (!silent) {
          cat("-- STARTING BOOTSTRAPPING TO COMPUTE CONFIDENCE INTERVALS! --\n");
          cat("--    (this might take a while, computing", samples,"samples)    --\n");
        }
      }
      
      ### Compute CI for alpha
      res$intermediate$alpha.ci <- ci.reliability(res$input$dat, type="alpha",
                                                  conf.level = conf.level,
                                                  interval.type=interval.type, B=samples);
      
      ### Compute CI for omega
      res$intermediate$omega.ci <- ci.reliability(res$input$dat, type="omega",
                                                  conf.level = conf.level,
                                                  interval.type=interval.type, B=samples);

      ### Confidence intervals for ordinal estimates
      if (poly && res$intermediate$maxLevels < 9 && res$intermediate$maxRange < 9) {
        
        ### Set interval to wald if bootstrapping was requested
        if (interval.type %in% c("bi", "perc", "bca")) {
          intervalType <- "normal-theory";
        } else {
          intervalType <- interval.type;
        }
        
        ### Compute CI and point estimates for ordinal alpha
        res$intermediate$alpha.ordinal.ci <- ci.reliability(S=res$intermediate$polychor,
                                                            N = res$input$n.observations,
                                                    type="alpha",
                                                    conf.level = conf.level,
                                                    interval.type=intervalType);
        
        ### Compute CI and point estimates for ordinal omega
        res$intermediate$omega.ordinal.ci <- ci.reliability(S=res$intermediate$polychor,
                                                            N = res$input$n.observations,
                                                    type="omega",
                                                    conf.level = conf.level,
                                                    interval.type=intervalType);
        res$output$dat$alpha.ordinal <- res$intermediate$alpha.ordinal.ci$est;
        res$output$dat$omega.ordinal <- res$intermediate$omega.ordinal.ci$est;
        
      }
      
      ### Extract and store ci bounds
      res$output$alpha.ci <- c(res$intermediate$alpha.ci$ci.lower,
                               res$intermediate$alpha.ci$ci.upper);
      res$output$omega.ci <- c(res$intermediate$omega.ci$ci.lower,
                               res$intermediate$omega.ci$ci.upper);
      res$output$alpha.ordinal.ci <- c(res$intermediate$alpha.ordinal.ci$ci.lower,
                                       res$intermediate$alpha.ordinal.ci$ci.upper);
      res$output$omega.ordinal.ci <- c(res$intermediate$omega.ordinal.ci$ci.lower,
                                       res$intermediate$omega.ordinal.ci$ci.upper);
      ### Add to dat
      res$output$dat$alpha.ci.lo <- res$output$alpha.ci[1];
      res$output$dat$alpha.ci.hi <- res$output$alpha.ci[2];
      res$output$dat$omega.ci.lo <- res$output$omega.ci[1];
      res$output$dat$omega.ci.hi <- res$output$omega.ci[2];
      res$output$dat$alpha.ordinal.ci.lo <- res$output$alpha.ordinal.ci[1];
      res$output$dat$alpha.ordinal.ci.hi <- res$output$alpha.ordinal.ci[2];
      res$output$dat$omega.ordinal.ci.lo <- res$output$omega.ordinal.ci[1];
      res$output$dat$omega.ordinal.ci.hi <- res$output$omega.ordinal.ci[2];
      if ((interval.type %in% c("bi", "perc", "bca")) && (!silent)) {
        cat("-- FINISHED BOOTSTRAPPING TO COMPUTE CONFIDENCE INTERVALS! --\n");
      }
    }
  }
  
  class(res) <- c("scaleStructure", "scaleReliability");
  ### Return result
  return(res);
}

print.scaleStructure <- function (x, digits=x$input$digits, ...) {
  
  if (packageVersion('psych') < '1.5.4') {
    cat("Note: your version of package 'psych' is lower than 1.5.4 (",
        as.character(packageVersion('psych')), " to be precise). This means that you ",
        "might see errors from the 'fa' function above this notice. ",
        "You can safely ignore these.\n\n", sep="");
  }
  
  cat(paste0("Information about this analysis:\n",
             "\n                 Dataframe: ", x$input$dat.name,
             "\n                     Items: ", paste(x$input$items, collapse=", "),
             "\n              Observations: ", x$input$n.observations,
             "\n     Positive correlations: ", x$intermediate$cor.pos,
             " out of ", x$intermediate$cor.total, " (",
             round(100*x$intermediate$cor.proPos), "%)\n\n",
             "Estimates assuming interval level:\n"));
  if (x$input$n.items > 2) {
    cat(paste0("\n             Omega (total): ", round(x$output$omega, digits=digits),
               "\n      Omega (hierarchical): ",
               round(x$intermediate$omega.psych$omega_h, digits=digits)));
    if (x$input$omega.psych) {
      cat(paste0("\nOmega (from psych package): ", round(x$output$omega.psych, digits=digits)));
    }
    cat(paste0("\nGreatest Lower Bound (GLB): ", round(x$output$glb, digits=digits),
               "\n          Cronbach's alpha: ", round(x$output$cronbach.alpha, digits=digits), "\n"));
    if (x$input$ci & !is.null(x$output$alpha.ci)) {
      ### If confidence intervals were computed AND obtained, print them
      cat(paste0("Confidence intervals:\n             Omega (total): [",
                 round(x$output$omega.ci[1], digits=digits), ", ",
                 round(x$output$omega.ci[2], digits=digits), "]\n",
                 "          Cronbach's alpha: [", round(x$output$alpha.ci[1], digits=digits),
                 ", ", round(x$output$alpha.ci[2], digits=digits), "]\n"));
    }
    if (x$input$poly && x$intermediate$maxLevels < 9 && x$intermediate$maxRange < 9) {
      if (!is.null(x$intermediate$omega.ordinal)) {
        cat(paste0("\nEstimates assuming ordinal level:\n",
                   "\n     Ordinal Omega (total): ",
                   round(x$intermediate$omega.ordinal$omega.tot, digits=digits),
                   "\n Ordinal Omega (hierarch.): ",
                   round(x$intermediate$omega.ordinal$omega_h, digits=digits)));
        if (x$input$omega.psych) {
          cat(paste0("\nOrd. Omega (psych package): ", round(x$intermediate$omega.ordinal$omega.tot, digits=digits)));
        }
        cat(paste0("\n  Ordinal Cronbach's alpha: ",
                   round(x$intermediate$alpha.ordinal$total$raw_alpha, digits=digits), "\n"));
        if (x$input$ci & !is.null(x$output$alpha.ordinal.ci)) {
          ### If confidence intervals were computed AND obtained, print them
          cat(paste0("Confidence intervals:\n     Ordinal Omega (total): [",
                     round(x$output$omega.ordinal.ci[1], digits=digits), ", ",
                     round(x$output$omega.ordinal.ci[2], digits=digits), "]\n",
                     "  Ordinal Cronbach's alpha: [", round(x$output$alpha.ordinal.ci[1], digits=digits),
                     ", ", round(x$output$alpha.ordinal.ci[2], digits=digits), "]\n"));
        }
      } else {
        cat0("\n(Estimates assuming ordinal level not computed, as the polychoric ",
             "correlation matrix has missing values.)\n");
      }
    } else if (x$input$poly == TRUE){
      cat("\n(Estimates assuming ordinal level not computed, as at least one item seems to have more than 8 levels; ",
          "the highest number of distinct levels is ", x$intermediate$maxLevels, " and the highest range is ",
          x$intermediate$maxRange, ". This last number needs to be lower than 9 for the polychoric function to work. ",
          "If this is unexpected, you may want to check for outliers.)\n", sep="");
    }
    if (x$input$omega.psych) {
      cat(paste0("\nNote: the normal point estimate and confidence interval for omega are based on the procedure suggested by ",
                 "Dunn, Baguley & Brunsden (2013) using the MBESS function ci.reliability, whereas the psych package point estimate was ",
                 "suggested in Revelle & Zinbarg (2008). See the help ('?scaleStructure') for more information.\n"));
    } else {
      cat(paste0("\nNote: the normal point estimate and confidence interval for omega are based on the procedure suggested by ",
                 "Dunn, Baguley & Brunsden (2013). To obtain the (usually higher) omega point estimate using the procedure ",
                 "suggested by Revelle & Zinbarg (2008), use argument 'omega.psych=TRUE'. See the help ('?scaleStructure') ",
                 "for more information. Of course, you can also call the 'omega' function from the psych package directly.\n"));
    }
  } else if (x$input$n.items == 2) {
    cat(paste0("\nSpearman Brown coefficient: ", round(x$output$spearman.brown, digits=digits),
               "\n          Cronbach's alpha: ", round(x$output$cronbach.alpha, digits=digits),
               "\n       Pearson Correlation: ", round(x$intermediate$cor[1, 2], digits=digits), "\n\n"));
  }
  invisible();
}
