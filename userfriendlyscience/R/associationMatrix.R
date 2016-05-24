###########################################################
### Define main functions
###########################################################

### First, we define the functions that will compute the
### results. We have functions for the different combinations
### of measurement levels of the two variables.

### Function for the t-test
computeStatistic_t <- function(var1, var2, conf.level=.95) {
  if (nlevels(as.factor(var1)) == 2) {
    dichotomous <- factor(var1);
    interval <- var2;
  }
  else if (nlevels(as.factor(var2)) == 2) {
    dichotomous <- factor(var2);
    interval <- var1;
  }
  else {
    stop("Error: none of the two variables has only two levels!");
  }
  res <- list();
  res$object <- meanDiff(interval ~ dichotomous, conf.level = conf.level,
                         envir=environment());
  res$statistic <- res$object$t;
  res$statistic.type <- "t";
  res$parameter <- res$object$df;
  res$p.raw <- res$object$p;
  return(res);
}

### Function for the Pearson correlation (r)
computeStatistic_r <- function(var1, var2, conf.level=.95) {
  res <- list();
  res$object <- cor.test(var1, var2, use="complete.obs");
  res$statistic <- res$object$statistic;
  res$statistic.type <- "r";
  res$parameter <- res$object$parameter;
  res$p.raw <- res$object$p.value;
  return(res);
}

### Function for Anova (f)
computeStatistic_f <- function(var1, var2, conf.level=.95) {
  
  if (is.factor(var1) & is.numeric(var2)) {
    factor <- var1;
    dependent <- var2;
  } else if (is.factor(var2) & is.numeric(var1)) {
    factor <- var2;
    dependent <- var1;
  } else if (nlevels(as.factor(var1)) < nlevels(as.factor(var2))) {
    ## We will treat var1 as factor
    factor <- factor(var1);
    dependent <- as.numeric(var2);
  }
  else {
    factor <- factor(var2);
    dependent <- as.numeric(var1);
  }
  
  ### In the future perhaps include tests of
  ### homogeneity of variances:
  #bartlett.test(dependent ~ factor);
  #fligner.test(dependent ~ factor);
  #safeRequire('car');
  #levene.test(dependent ~ factor);
  
  res <- list();
  res$object <- aov(dependent ~ factor);
  
  res$statistic <- summary(res$object)[[1]][['F value']][1];
  res$statistic.type <- "f";
  res$parameter <- c(summary(res$object)[[1]][['Df']]);
#                    summary(res$object)[[2]][['Df']]);
  res$p.raw <- summary(res$object)[[1]][['Pr(>F)']][1];
  return(res);
}

### Function for chi-square (chisq)
computeStatistic_chisq <- function(var1, var2, conf.level=.95) {
  res <- list();
  res$object <- chisq.test(var1, var2, correct=FALSE);
  res$statistic <- res$object$statistic;
  res$statistic.type <- "chisq";
  res$parameter <- res$object$parameter;
  res$p.raw <- res$object$p.value;
  return(res);
}

### Effect size Cohens d
computeEffectSize_d <- function(var1, var2, conf.level=.95) {
  if (nlevels(as.factor(var1)) == 2) {
    dichotomous <- factor(var1);
    interval <- var2;
  }
  else if (nlevels(as.factor(var2)) == 2) {
    dichotomous <- factor(var2);
    interval <- var1;
  }
  else {
    stop("Error: none of the two variables has only two levels!");
  }
  res <- list();
  res$object <- meanDiff(interval ~ dichotomous, conf.level = conf.level,
                         envir=environment());
  res$es <- res$object$meanDiff.g;
  res$es.type <- "g";
  res$ci <- c(res$object$meanDiff.g.ci.lower,
              res$object$meanDiff.g.ci.upper);
  return(res);
}

### Effect size Pearson's r
computeEffectSize_r <- function(var1, var2, conf.level=.95) {
  res <- list();
  res$object <- cor.test(var1, var2, use="complete.obs");
  res$es <- res$object$estimate;
  res$es.type <- "r";
  res$ci <- res$object$conf.int;
  return(res);
}

### Function for eta squared (etasq)
computeEffectSize_etasq <- function(var1, var2, conf.level=.95) {
  
  if (is.factor(var1) & is.numeric(var2)) {
    factor <- var1;
    dependent <- var2;
  }
  else if (is.factor(var2) & is.numeric(var1)) {
    factor <- var2;
    dependent <- var1;
  } else if (nlevels(as.factor(var1)) < nlevels(as.factor(var2))) {
    ## We will treat var1 as factor
    factor <- factor(var1);
    dependent <- as.numeric(var2);
  }
  else {
    factor <- factor(var2);
    dependent <- as.numeric(var1);
  }
  
  res <- list();

  ### Confidence level should be doubled (i.e. unconfidence level
  ### should be doubled to be more precise), so .95 becomes .90 ->
  ### see http://daniellakens.blogspot.nl/2014/06/calculating-confidence-intervals-for.html
  ### for a brief explanation and links to more extensive explanations.
  
  res$realConfidence <- 1 - ((1-conf.level) * 2);
  
  res$object.aov <- aov(dependent ~ factor);
  
  df_num <- summary(res$object)[[1]][1,1];
  df_den <- summary(res$object)[[1]][2,1];
  f_val <- summary(res$object)[[1]][1,4];
  
  ### This is suggested by the page at
  ### http://yatani.jp/HCIstats/ANOVA#RCodeOneWay
  ### (capture.output used because this function for
  ###  some reason very tenaciously outputs results)
  ### (also note that we double the 'unconfidence' level,
  ###  e.g. conf.level=.95 becomes conf.level=.90, to
  ###  retain consistency with the NHST p-value; see
  ###  the Word doc by Karl Wuensch references above,
  ###  or the paper he cites:
  ###    Steiger, J. H. (2004). Beyond the F test:
  ###      Effect size confidence intervals and tests
  ###      of close fit in the analysis of variance and
  ###      contrast analysis. Psychological methods, 9(2),
  ###      164-82. doi:10.1037/1082-989X.9.2.164
  
  res$es <- df_num*f_val/(df_den + df_num*f_val);
  res$es.type <- "etasq";
  capture.output(res$object <- ci.pvaf(F.value=f_val, df.1=df_num, df.2=df_den,
                        N=(df_den+df_num+1), conf.level=res$realConfidence));
  
  res$ci <- c(res$object$Lower.Limit.Proportion.of.Variance.Accounted.for,
              res$object$Upper.Limit.Proportion.of.Variance.Accounted.for);
  
  return(res);
}  

### Function for Cramers V effect size (v)
computeEffectSize_v <- function(var1, var2, conf.level=.95,
                                     bootstrap=FALSE, samples=5000) {
  res <- list();
  if (bootstrap) {
    res$object <- confIntV(var1, var2,
                           method="bootstrap",
                           samples=samples);
    res$ci <- res$object$output$confIntV.bootstrap;
  } else {
    res$object <- confIntV(var1, var2, method="fisher");
    res$ci <- res$object$output$confIntV.fisher;
  }
  res$es <- res$object$intermediate$cramersV$output$cramersV
  res$es.type <- "V";
  return(res);
}

### This is the function that calls the functions
### to compute statistics and effect sizes, and
### organises the resulting objects in sets of
### lists. The elements of the first list are
### the 'rows' of the matrix. Each element (each
### 'row') is itself again a list, where each
### element corresponds to a 'cell' in the
### final 'matrix'. Each of these element (each
### of these cells) contains two objects; the one
### containing the statistic and the one
### containing the effect size.

associationMatrixStatDefaults <- list(dichotomous =
                                        list(dichotomous = "computeStatistic_chisq",
                                             nominal = "computeStatistic_chisq",
                                             ordinal = "computeStatistic_chisq",
                                             numeric = "computeStatistic_t"),
                                      nominal =
                                        list(dichotomous = "computeStatistic_chisq",
                                             nominal = "computeStatistic_chisq",
                                             ordinal = "computeStatistic_chisq",
                                             numeric = "computeStatistic_f"),
                                      ordinal =
                                        list(dichotomous = "computeStatistic_chisq",
                                             nominal = "computeStatistic_chisq",
                                             ordinal = "computeStatistic_chisq",
                                             numeric = "computeStatistic_f"),
                                      numeric =
                                        list(dichotomous = "computeStatistic_t",
                                             nominal = "computeStatistic_f",
                                             ordinal = "computeStatistic_f",
                                             numeric = "computeStatistic_r"));
associationMatrixESDefaults <- list(dichotomous =
                                      list(dichotomous = "computeEffectSize_v",
                                           nominal = "computeEffectSize_v",
                                           ordinal = "computeEffectSize_v",
                                           numeric = "computeEffectSize_d"),
                                    nominal =
                                      list(dichotomous = "computeEffectSize_v",
                                           nominal = "computeEffectSize_v",
                                           ordinal = "computeEffectSize_v",
                                           numeric = "computeEffectSize_etasq"),
                                    ordinal =
                                      list(dichotomous = "computeEffectSize_v",
                                           nominal = "computeEffectSize_v",
                                           ordinal = "computeEffectSize_v",
                                           numeric = "computeEffectSize_etasq"),
                                    numeric =
                                      list(dichotomous = "computeEffectSize_d",
                                           nominal = "computeEffectSize_etasq",
                                           ordinal = "computeEffectSize_etasq",
                                           numeric = "computeEffectSize_r"));

associationMatrix <- function(dat=NULL, x=NULL, y=NULL, conf.level = .95,
                              correction = "fdr", bootstrapV=FALSE,
                              info=c("full", "ci", "es"),
                              includeSampleSize = "depends",
                              bootstrapV.samples = 5000, digits = 2,
                              pValueDigits=digits + 1, colNames = FALSE,
                              type=c("R", "html", "latex"), file="",
                              statistic = associationMatrixStatDefaults,
                              effectSize = associationMatrixESDefaults) {
  
  ### Make object to store results
  res <- list(input = as.list(environment()),
              intermediate = list(),
              output = list());
  res$intermediate$statistics <- list();
  res$intermediate$effectSizes <- list();
  res$intermediate$sampleSizes <- list();
  
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
  
  ### If no variables are specified, take them all.
  if (is.null(x) && is.null(y)) {
    x <- names(dat);
  }
  
  ### If y was accidently specified, but x wasn't, copy y to x.
  if (is.null(x) && !is.null(y)) {
    x <- y;
  }
  
  ### Check whether the first vector of variable names has sufficient elements.
  if (length(x) < 1) {
    stop(paste0("Error: x vector has 0 elements or less; ",
                "make sure to specify at least one variable name!."));
  }
  
  ### Check and store the measurement level of the variables:
  ### dichotomous, nominal, ordinal, or interval
  measurementLevelsX <- vector();
  xCounter <- 1;
  for(curXvar in x) {
    if (var(as.numeric(dat[,curXvar]), na.rm=TRUE) == 0) {
      stop("Variable '", curXvar, "' has no variance (everybody scores the same)! ",
           "This prohibits the calculation of effect size measures, so I'm aborting.");
    }
    if (is.numeric(dat[,curXvar])) {
      measurementLevelsX[xCounter] <- "numeric";
    }
    else if (is.factor(dat[,curXvar])) {
      if (length(levels(dat[,curXvar])) == 2) {
        measurementLevelsX[xCounter] <- "dichotomous";
      }        
      else if (is.ordered(dat[,curXvar])) {
        measurementLevelsX[xCounter] <- "ordinal";
      }
      else {
        measurementLevelsX[xCounter] <- "nominal";      
      }
    }
    else {
      stop(paste0("Error: variable '", curXvar, "'' does not have ",
                  "nominal, ordinal, or interval measurement level!"));
    }
    xCounter <- xCounter + 1;
  }
  
  ### Check whether we have a second set of variables
  if (!is.null(y)) {
    ### Check whether the second vector of variable names has sufficient elements.
    if (length(y) < 1) {
      stop(paste0("Error: y vector has 0 elements or less; ",
                  "make sure to specify at least one variable name!."));
    }
    symmetric <- FALSE;
    ### Check and store the measurement level of the variables:
    ### dichotomous, nominal, ordinal, or interval
    measurementLevelsY <- vector();
    yCounter <- 1;
    for(curYvar in y) {
      if (var(as.numeric(dat[,curYvar]), na.rm=TRUE) == 0) {
        stop("Variable '", curYvar, "' has no variance (everybody scores the same)! ",
             "This prohibits the calculation of effect size measures, so I'm aborting.");
      }
      if (is.numeric(dat[,curYvar])) {
        measurementLevelsY[yCounter] <- "numeric";
      }
      else if (is.factor(dat[,curYvar])) {
        if (length(levels(dat[,curYvar])) == 2) {
          measurementLevelsY[yCounter] <- "dichotomous";
        }        
        else if (is.ordered(dat[,curYvar])) {
          measurementLevelsY[yCounter] <- "ordinal";
        }
        else {
          measurementLevelsY[yCounter] <- "nominal";      
        }
      }
      else {
        stop(paste0("Error: variable '", curYvar, "'' does not have ",
                    "nominal, ordinal, or interval measurement level!"));
      }
      yCounter <- yCounter + 1;
    }
  }
  else {
    symmetric <- TRUE;
    y <- x;
    measurementLevelsY <- measurementLevelsX;
  }
  
  ### Generate vectors with row and column names
  if (colNames) {
    rowNames <- x;
    columnNames <- y;
  }
  else {
    rowNames <- paste(1:length(x), x, sep=". ");
    columnNames <- paste0(1:length(y), ".");
  }
  
  ### Generate matrices for results and set row and column names
  res$output$matrix <- list();
  res$output$matrix$es <- matrix(nrow = length(x), ncol = length(y));
  rownames(res$output$matrix$es) <- rowNames;
  colnames(res$output$matrix$es) <- columnNames;
  res$output$matrix$sampleSizes <- matrix(nrow = length(x), ncol = length(y));
  rownames(res$output$matrix$sampleSizes) <- rowNames;
  colnames(res$output$matrix$sampleSizes) <- columnNames;
  res$output$matrix$ci <- matrix(nrow = length(x), ncol = length(y));
  rownames(res$output$matrix$ci) <- rowNames;
  colnames(res$output$matrix$ci) <- columnNames;
  res$output$matrix$full <- matrix(nrow = 2 * length(x), ncol = length(y));
  rownames(res$output$matrix$full) <- rep("", 2*length(rowNames));
  rownames(res$output$matrix$full)[seq(1, (2*length(rowNames)) - 1, by=2)] <- rowNames;
  colnames(res$output$matrix$full) <- columnNames;
  
  xCounter <- 1;
  for(curXvar in x) {
    ### For each row, create the object (list) that will
    ### contain the cells
    res$intermediate$statistics[[curXvar]] <- list();
    res$intermediate$effectSizes[[curXvar]] <- list();
    res$intermediate$sampleSizes[[curXvar]] <- list();
    yCounter <- 1;
    for(curYvar in y) {
      ### If a symmetric table was requested, don't do
      ### anything unless we're in the lower left half.
      if (!symmetric | (yCounter < xCounter)) {
        
        ### Call the function to compute the statistic.
        ### Which function this is, depends on the preferences
        ### of the user (or the defaults).
        tmpFun <- match.fun(statistic[[measurementLevelsX[xCounter]]]
                                     [[measurementLevelsY[yCounter]]]);
        res$intermediate$statistics[[curXvar]][[curYvar]] <- 
          tmpFun(dat[,curXvar], dat[,curYvar], conf.level = conf.level);
        ### We repeat the same trick for the effect sizes.
        tmpFun <- match.fun(effectSize[[measurementLevelsX[xCounter]]]
                            [[measurementLevelsY[yCounter]]]);
        res$intermediate$effectSizes[[curXvar]][[curYvar]] <- 
          tmpFun(dat[,curXvar], dat[,curYvar], conf.level = conf.level);
        res$intermediate$sampleSizes[[curXvar]][[curYvar]] <- nrow(na.omit(dat[,c(curXvar, curYvar)]));
      }
      yCounter <- yCounter + 1;
    }
    xCounter <- xCounter + 1;
  }
  
  ### Correct p-values for multiple testing
  ### First build a matrix with the raw p-values
  res$intermediate$pvalMatrix <- matrix(nrow=length(x), ncol=length(y), dimnames=list(x, y));
  for(curXvar in x) {
    for(curYvar in y) {
      if (!is.null(res$intermediate$statistics[[curXvar]][[curYvar]]$p.raw)) {
        res$intermediate$pvalMatrix[curXvar, curYvar] <- res$intermediate$statistics[[curXvar]][[curYvar]]$p.raw;
      }
    }
  }
  ### Adjust p-values
  res$intermediate$pvalMatrix.adj <- matrix(p.adjust(res$intermediate$pvalMatrix, method=correction),
                               nrow(res$intermediate$pvalMatrix), ncol(res$intermediate$pvalMatrix),
                               dimnames=dimnames(res$intermediate$pvalMatrix));
  ### Store adjusted p-values in objects
  for(curXvar in x) {
    for(curYvar in y) {
      if (!is.null(res$intermediate$statistics[[curXvar]][[curYvar]]$p.raw)) {
        res$intermediate$statistics[[curXvar]][[curYvar]]$p.adj <-
          res$intermediate$pvalMatrix.adj[curXvar, curYvar];
      }
    }
  }

  ### Run the loop again to create the output matrices (one with point
  ### estimates and p-values corrected for multiple testing; one with
  ### confidence intervals; and one with two rows for each variable,
  ### combining the information).

  for(rowVar in 1:length(x)) {
    for(colVar in 1:length(y)) {
      ### If a symmetric table was requested, only fill the cells if we're
      ### in the lower left half.
      if (!symmetric | (colVar < rowVar)) {
        ### Extract and set confidence interval and then es estimate & p value
        res$output$matrix$ci[rowVar, colVar] <- paste0(
          substr(res$intermediate$effectSizes[[rowVar]][[colVar]]$es.type, 1, 1), "=[",
               round(res$intermediate$effectSizes[[rowVar]][[colVar]]$ci[1], digits), "; ",
               round(res$intermediate$effectSizes[[rowVar]][[colVar]]$ci[2], digits), "]");
        res$output$matrix$es[rowVar, colVar] <-
          paste0(substr(res$intermediate$effectSizes[[rowVar]][[colVar]]$es.type, 1, 1), "=",
                 round(res$intermediate$effectSizes[[rowVar]][[colVar]]$es, digits), ", ",
                 formatPvalue(res$intermediate$statistics[[rowVar]][[colVar]]$p.adj, digits=pValueDigits, spaces=FALSE));
        res$output$matrix$sampleSizes[rowVar, colVar] <-
          res$intermediate$sampleSizes[[rowVar]][[colVar]];
        ### Convert x (row variable) to two row indices in combined matrix
        res$output$matrix$full[(rowVar*2)-1, colVar] <-
          res$output$matrix$ci[rowVar, colVar];
        res$output$matrix$full[(rowVar*2), colVar] <-
          res$output$matrix$es[rowVar, colVar];
        if (((includeSampleSize == "depends") &&
            (length(unique(unlist(res$intermediate$sampleSizes))) > 1)) ||
            (includeSampleSize == "always")) {
          res$output$matrix$full[(rowVar*2), colVar] <-
            paste0(res$output$matrix$full[(rowVar*2), colVar],
                   ", n=", res$output$matrix$sampleSizes[rowVar, colVar]);
        }
      }
      else {
        res$output$matrix$es[rowVar, colVar] <- "";
        res$output$matrix$ci[rowVar, colVar] <- "";
        ### Convert x (row variable) to two row indices in combined matrix
        res$output$matrix$full[(rowVar*2)-1, colVar] <- "";
        res$output$matrix$full[rowVar*2, colVar] <- "";
      }
    }
  }
  
  ### Set class & return result
  class(res) <- c("associationMatrix");
  return(res);
}

print.associationMatrix <- function (x, type = x$input$type,
                                     info = x$input$info, 
                                     file = x$input$file, ...) {
  
  ### Extract matrix to print (es, ci, or full)
  matrixToPrint <- x$output$matrix[[info[1]]];
  
  ### Either show in R, or convert to html or latex
  if (toupper(type[1])=="R") {
    if (file=="") {
      print(matrixToPrint, quote=FALSE);
    } else {
      write.table(matrixToPrint, file=file, sep="\t",
                  quote=FALSE, row.names=TRUE, col.names=TRUE);
    }
  }
  else {
    ### Replace even row names (currently empty) with unique comments
    ### for html and LaTeX
    if ((tolower(type[1])=="latex") && (info[1]=='full')) {
      rownames(matrixToPrint)[seq(2, nrow(matrixToPrint), by=2)] <-
        paste0("%%% Variable ", 1:(nrow(matrixToPrint)/2), "\n");
    } else if (info[1]=='full') {
      rownames(matrixToPrint)[seq(2, nrow(matrixToPrint), by=2)] <-
        paste0("<!-- Variable ", 1:(nrow(matrixToPrint)/2), " -->");
    }
    print(xtable(matrixToPrint, align=c('l', rep('c', ncol(matrixToPrint)))),
          type=type,
          html.table.attributes = "cellpadding='5' style='border=0 solid;'",
          sanitize.text.function=function(x) { return (x); },
          file=file);
  }
    
  invisible();
}
