regr <- function(formula, dat=NULL, conf.level=.95, digits=2,
                 pvalueDigits = 3, coefficients=c("raw", "scaled"),
                 plot=FALSE) {

  ### Generate object to store input, intermediate outcomes, and results
  res <- list(input = as.list(environment()),
              intermediate = list(),
              output = list());
  
  ### Extract variables from formula
  res$intermediate$variableNames <- all.vars(formula);
  
  ### Convert formula to a character string
  res$intermediate$formula.as.character <-
    paste0(as.character(formula)[c(2, 1, 3)], collapse=" ");
  
  if (is.null(dat)) {
    ### Extract variables from formula
    res$intermediate$variables <-
      as.character(as.list(attr(terms(formula), 'variables'))[-1]);
    
    ### Store variablesnames only for naming in raw dataframe
    res$intermediate$variables_namesOnly <- unlist(
      lapply(strsplit(res$intermediate$variables, "\\$"), tail, 1));
    
    ### Store variables in a dataframe
    res$intermediate$dat.raw <- list();
    for (varToGet in 1:length(res$intermediate$variables)) {
      res$intermediate$dat.raw[[res$intermediate$variables_namesOnly[varToGet]]] <-
        eval(parse(text=res$intermediate$variables[varToGet]), envir=parent.env(environment()));
    }
    ### Convert the list to a dataframe
    res$intermediate$dat.raw <- data.frame(res$intermediate$dat.raw);
    
    ### Convert variable names to bare variable names
    for (currentVariableIndex in 1:length(res$intermediate$variables)) {
      res$intermediate$formula.as.character <-
        gsub(res$intermediate$variables[currentVariableIndex],
             res$intermediate$variables_namesOnly[currentVariableIndex],
             res$intermediate$formula.as.character, fixed=TRUE);  
    }
    
  } else {
    ### Store variables in a dataframe
    res$intermediate$dat.raw <- dat[, res$intermediate$variableNames];
    res$intermediate$variables_namesOnly <- res$intermediate$variableNames;
  }

  res$intermediate$formula <- formula(res$intermediate$formula.as.character);
  
  ### Standardize variables
  res$intermediate$dat.scaled <- as.data.frame(scale(res$intermediate$dat.raw));
  
  ### Run and store lm objects
  res$intermediate$lm.raw <-
    lm(res$intermediate$formula, res$intermediate$dat.raw);
  res$intermediate$lm.scaled <-
    lm(res$intermediate$formula, res$intermediate$dat.scaled);
  
  ### R^2 confidence interval based on formula at
  ### http://www.danielsoper.com/statcalc3/calc.aspx?id=28
  res$intermediate$rsq <- rsq <- summary(res$intermediate$lm.raw)$r.squared;
  res$intermediate$k <- k <- length(res$intermediate$lm.raw$terms) - 1;
  res$intermediate$n <- n <- res$intermediate$lm.raw$df.residual + k + 1;
  res$intermediate$rsq.se <- sqrt((4*rsq*(1-rsq)^2*(n-k-1)^2)/
                                  ((n^2-1)*(3+n)));
  res$intermediate$rsq.t.crit <- qt(p=1-(1-conf.level)/2, df=n-k-1);
  res$output$rsq.ci <- c(res$intermediate$rsq -
                           res$intermediate$rsq.t.crit *
                           res$intermediate$rsq.se,
                         res$intermediate$rsq +
                           res$intermediate$rsq.t.crit *
                           res$intermediate$rsq.se);
  
  ### Run confint on lm object
  res$intermediate$confint.raw <-
    confint(res$intermediate$lm.raw, level=conf.level);
  res$intermediate$confint.scaled <-
    confint(res$intermediate$lm.scaled, level=conf.level);
  
  ### Run lm.influence on lm object
  res$intermediate$lm.influence.raw <-
    lm.influence(res$intermediate$lm.raw);
  res$intermediate$lm.influence.scaled <-
    lm.influence(res$intermediate$lm.scaled);

  ### get summary for lm objects
  res$intermediate$summary.raw <- summary(res$intermediate$lm.raw);
  res$intermediate$summary.scaled <- summary(res$intermediate$lm.scaled);
  
  ### Generate output
  res$output$coef.raw <- cbind(data.frame(res$intermediate$confint.raw[, 1],
                                          res$intermediate$confint.raw[, 2]),
                               res$intermediate$summary.raw$coefficients);
  res$output$coef.scaled <- cbind(data.frame(res$intermediate$confint.scaled[, 1],
                                             res$intermediate$confint.scaled[, 2]),
                                  res$intermediate$summary.scaled$coefficients);
  
  names(res$output$coef.raw) <- names(res$output$coef.scaled) <-
    c(paste0(conf.level*100,"% CI, lo"),
      paste0(conf.level*100,"% CI, hi"),
      'estimate', 'se', 't', 'p');
  
  if (plot) {
    if (length(res$intermediate$variables_namesOnly) == 2) {
      res$output$plot <- ggplot(res$intermediate$dat.raw,
                                aes_string(y=res$intermediate$variables_namesOnly[1],
                                        x=res$intermediate$variables_namesOnly[2])) +
        geom_point() + geom_smooth(method='lm') + theme_bw();
    } else {
      warning("You requested a plot, but for now plots are only available for regression analyses with one predictor.");
    }
  }
  
  ### Use Z = (b1 - b2) / sqrt(SE_b1^2 + SE_b2^2) to test whether
  ### the coefficients differ; add arguments such as
  ### compareCoefficients=FALSE, p.adjust="fdr" to enable user to
  ### request & correct this.
  
  class(res) <- 'regr';
  return(res);

}

print.regr <- function(x, digits=x$input$digits,
                       pvalueDigits=x$input$pvalueDigits, ...) {

  cat(paste0("Regression analysis for formula: ",
             x$intermediate$formula.as.character, "\n\n",
             "Significance test of the entire model (all predictors together):\n",
             "  Multiple R-squared: [",
             round(x$output$rsq.ci[1], digits), ", ",
             round(x$output$rsq.ci[2], digits),
             "] (point estimate = ",
             round(x$intermediate$summary.raw$r.squared, digits),
             ", adjusted = ",
             round(x$intermediate$summary.raw$adj.r.squared, digits), ")\n",
             "  Test for significance: F[",
             x$intermediate$summary.raw$fstatistic[2], ", ",
             x$intermediate$summary.raw$fstatistic[3], "] = ",
             round(x$intermediate$summary.raw$fstatistic[1], digits),
             ", ", formatPvalue(pf(x$intermediate$summary.raw$fstatistic[1],
                                   x$intermediate$summary.raw$fstatistic[2],
                                   x$intermediate$summary.raw$fstatistic[3],
                                   lower.tail=FALSE), digits=pvalueDigits), "\n"));
  if ("raw" %in% x$input$coefficients) {
    cat("\nRaw regression coefficients (unstandardized beta values, called 'B' in SPSS):\n\n");
    tmpDat <- round(x$output$coef.raw[, 1:5], digits);
    tmpDat[[1]] <- paste0("[", tmpDat[[1]], "; ", tmpDat[[2]], "]");
    tmpDat[[2]] <- NULL;
    names(tmpDat)[1] <- paste0(x$input$conf.level*100, "% conf. int.");
    tmpDat$p <- formatPvalue(x$output$coef.raw$p,
                             digits=pvalueDigits,
                             includeP=FALSE);
    print(tmpDat, ...);
  }
  if ("scaled" %in% x$input$coefficients) {
    cat("\nScaled regression coefficients (standardized beta values, called 'Beta' in SPSS):\n\n");
    tmpDat <- round(x$output$coef.scaled[, 1:5], digits);
    tmpDat[[1]] <- paste0("[", tmpDat[[1]], "; ", tmpDat[[2]], "]");
    tmpDat[[2]] <- NULL;
    names(tmpDat)[1] <- paste0(x$input$conf.level*100, "% conf. int.");
    tmpDat$p <- formatPvalue(x$output$coef.scaled$p,
                             digits=pvalueDigits,
                             includeP=FALSE);
    print(tmpDat, ...);
  }
  if (!is.null(x$output$plot)) {
    print(x$output$plot);
  }
  invisible();
  
}