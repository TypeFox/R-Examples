###########################################################
###
### R file with the function meanDiff.multi, which
### explores difference in means for a combination of 
### dichotomous and interval variables, generating a
### dataframe with the results and a number of plots.
###
### File created by Gjalt-Jorn Peters. Questions? You can
### contact me through http://behaviorchange.eu. Additional
### functions can be found at http://github.com/matherion
###
###########################################################

### Note: this is necessary to prevent Rcmd CHECK from throwing a note;
### otherwise it think these variable weren't defined yet.
utils::globalVariables(c('g', 'g.ci.lo', 'g.ci.hi', 'ci.lo', 'ci.hi',
                         'coord_flip', 'geom_hline'));

###########################################################
### Define functions
###########################################################

meanDiff.multi <- function(dat, y, x=NULL, var.equal = "yes",
                           conf.level = .95, digits = 2,
                           orientation = "vertical",
                           zeroLineColor = "grey",
                           zeroLineSize = 1.2,
                           envir = parent.frame()) {

  ### Check basic arguments
  if (!is.character(var.equal) || length(var.equal) != 1) {
    stop("Argument 'var.equal' must be either 'test', 'yes', or 'no'!");
  }
  var.equal <- tolower(var.equal);
  
  if (!is.numeric(conf.level) || length(conf.level) != 1 || conf.level < 0 || conf.level > 1) {
    stop("Argument 'conf.level' must be a number between 0 and 1!");
  }
  
  if (!is.numeric(digits) || length(digits) != 1 || digits < 0 || (round(digits) != digits)) {
    stop("Argument 'digits' must be a whole, positive number!");
  }
  
  ### Check whether y is a character vector
  if (!is.vector(y) | length(y) < 1 | !is.character(y)) {
    stop("y argument must be a character vector with at least ",
         "one variable name!");
  }
  
  ### Create object to store results
  res <- list();
  res$results.raw <- list();
  res$plots <- list();
  res$dlvPlots <- list();
  res$results.compiled <- data.frame();
  res$plots.compiled <- list();
  res$input <- list();
  res$input$dat <- dat;
  res$input$x <- x;
  res$input$y <- y;
  res$input$var.equal <- var.equal;
  res$input$conf.level <- conf.level;
  res$input$digits <- digits;
  res$input$envir <- envir;
  
  ### Check whether we have to do independent or dependent t-tests
  if (is.null(x)) {
    ### Dependent t-tests
    ### Check whether we have at least two variables to compare
    if (length(y) < 2) {
      stop("For a dependent test, the y argument must be a character ",
           "vector with at least two variable names!");
    }
    for (xVar in 1:length(y)) {
      for (yVar in 2:length(y)) {
        res$results.raw[[length(res$results.raw) + 1]] <-
          meanDiff(x=dat[[x[xVar]]], y=dat[[y[yVar]]], paired=TRUE,
                   conf.level=conf.level, digits=digits, envir=envir);
        tempData <- data.frame(x=factor(c(rep(0, length(dat[[x[xVar]]])),
                                          rep(1, length(dat[[y[yVar]]])))),
                               y=c(dat[[x[xVar]]], dat[[y[yVar]]]));
        ### Build dataframe for plot
        tempData <-
          data.frame(y = c(res$results.raw[[length(res$results.raw)]]$x,
                           res$results.raw[[length(res$results.raw)]]$y));
        tempData$x <- factor(c(rep(0, length(res$results.raw[[length(res$results.raw)]]$x)),
                               rep(1, length(res$results.raw[[length(res$results.raw)]]$y))),
                             levels=c(0,1), labels=c(res$results.raw[[length(res$results.raw)]]$groups[1],
                                                     res$results.raw[[length(res$results.raw)]]$groups[2]));
        ### Store plot
        res$plots[[length(res$plots) + 1]] <-
          ggplot(tempData, aes(y=y, x=x, group=x)) + geom_point() +
          stat_summary(fun.data = "mean_cl_boot", colour = "red");
        
        ### Extract information for compilation
        res$results.compiled[nrow(res$results.compiled)+1, 'x'] <- x[xVar];
        res$results.compiled[nrow(res$results.compiled), 'y'] <- y[yVar];
        
        res$results.compiled[nrow(res$results.compiled), 'group1'] <-
          res$results.raw[[length(res$results.raw)]]$groups[1];
        res$results.compiled[nrow(res$results.compiled), 'group2'] <-
          res$results.raw[[length(res$results.raw)]]$groups[2];
        res$results.compiled[nrow(res$results.compiled), 'mean1'] <-
          res$results.raw[[length(res$results.raw)]]$mean[1];
        res$results.compiled[nrow(res$results.compiled), 'mean2'] <-
          res$results.raw[[length(res$results.raw)]]$mean[2];
        res$results.compiled[nrow(res$results.compiled), 'sd1'] <-
          res$results.raw[[length(res$results.raw)]]$sd[1];
        res$results.compiled[nrow(res$results.compiled), 'sd2'] <-
          res$results.raw[[length(res$results.raw)]]$sd[2];
        res$results.compiled[nrow(res$results.compiled), 'n1'] <-
          res$results.raw[[length(res$results.raw)]]$n[1];
        res$results.compiled[nrow(res$results.compiled), 'n2'] <-
          res$results.raw[[length(res$results.raw)]]$n[2];
        
        res$results.compiled[nrow(res$results.compiled), 'g'] <-
          res$results.raw[[length(res$results.raw)]]$meanDiff.g;
        res$results.compiled[nrow(res$results.compiled), 'g.ci.lo'] <-
          res$results.raw[[length(res$results.raw)]]$meanDiff.g.ci.lower;
        res$results.compiled[nrow(res$results.compiled), 'g.ci.hi'] <-
          res$results.raw[[length(res$results.raw)]]$meanDiff.g.ci.upper;
        
        res$results.compiled[nrow(res$results.compiled), 'pwr.g'] <-
          res$results.raw[[length(res$results.raw)]]$power$power;
        res$results.compiled[nrow(res$results.compiled), 'pwr.small'] <-
          res$results.raw[[length(res$results.raw)]]$power.small$power;
        res$results.compiled[nrow(res$results.compiled), 'pwr.medium'] <-
          res$results.raw[[length(res$results.raw)]]$power.medium$power;
        res$results.compiled[nrow(res$results.compiled), 'pwr.large'] <-
          res$results.raw[[length(res$results.raw)]]$power.large$power;

        res$results.compiled[nrow(res$results.compiled), 't'] <-
          res$results.raw[[length(res$results.raw)]]$t;
        res$results.compiled[nrow(res$results.compiled), 'df'] <-
          res$results.raw[[length(res$results.raw)]]$df;
        res$results.compiled[nrow(res$results.compiled), 'p'] <-
          res$results.raw[[length(res$results.raw)]]$p;
        
      }
    }
  }
  else {
    ### Independent t-tests
    for (xVar in 1:length(x)) {
      for (yVar in 1:length(y)) {
        
        meanDiffFormula <- formula(paste0('dat$', y[yVar], " ~ ",
                                          'dat$', x[xVar]));
        
        res$results.raw[[length(res$results.raw) + 1]] <-
           meanDiff(x=meanDiffFormula, paired=FALSE,
                    var.equal=var.equal,
                    conf.level=conf.level, digits=digits, envir=envir);
        
        ### Build dataframe for plot
        tempData <-
          data.frame(y = c(res$results.raw[[length(res$results.raw)]]$x,
                           res$results.raw[[length(res$results.raw)]]$y));
        tempData$x <- factor(c(rep(0, length(res$results.raw[[length(res$results.raw)]]$x)),
                               rep(1, length(res$results.raw[[length(res$results.raw)]]$y))),
                             levels=c(0,1), labels=c(res$results.raw[[length(res$results.raw)]]$groups[1],
                                                     res$results.raw[[length(res$results.raw)]]$groups[2]));
        
        ### Store dlvPlot
        res$dlvPlots[[length(res$dlvPlots) + 1]] <- dlvPlot(tempData, y="y", x="x", conf.level=conf.level);

        ### Store regular plot with lines, use dlvPlot descriptives
        res$plots[[length(res$plots) + 1]] <- ggplot(tempData, aes(y=y, x=x, group=x)) +
          geom_point() +
          geom_pointrange(data=res$dlvPlots[[length(res$dlvPlots)]]$descr,
                          aes(x=x, y=mean, ymin=ci.lo, ymax=ci.hi), color="red");
                          
        ### Extract information for compilation
        res$results.compiled[nrow(res$results.compiled)+1, 'x'] <- x[xVar];
        res$results.compiled[nrow(res$results.compiled), 'y'] <- y[yVar];
        
        res$results.compiled[nrow(res$results.compiled), 'group1'] <-
          res$results.raw[[length(res$results.raw)]]$groups[1];
        res$results.compiled[nrow(res$results.compiled), 'group2'] <-
          res$results.raw[[length(res$results.raw)]]$groups[2];
        res$results.compiled[nrow(res$results.compiled), 'mean1'] <-
          res$results.raw[[length(res$results.raw)]]$mean[1];
        res$results.compiled[nrow(res$results.compiled), 'mean2'] <-
          res$results.raw[[length(res$results.raw)]]$mean[2];
        res$results.compiled[nrow(res$results.compiled), 'sd1'] <-
          res$results.raw[[length(res$results.raw)]]$sd[1];
        res$results.compiled[nrow(res$results.compiled), 'sd2'] <-
          res$results.raw[[length(res$results.raw)]]$sd[2];
        res$results.compiled[nrow(res$results.compiled), 'n1'] <-
          res$results.raw[[length(res$results.raw)]]$n[1];
        res$results.compiled[nrow(res$results.compiled), 'n2'] <-
          res$results.raw[[length(res$results.raw)]]$n[2];
        
        res$results.compiled[nrow(res$results.compiled), 'g'] <-
          res$results.raw[[length(res$results.raw)]]$meanDiff.g;
        res$results.compiled[nrow(res$results.compiled), 'g.ci.lo'] <-
          res$results.raw[[length(res$results.raw)]]$meanDiff.g.ci.lower;
        res$results.compiled[nrow(res$results.compiled), 'g.ci.hi'] <-
          res$results.raw[[length(res$results.raw)]]$meanDiff.g.ci.upper;
        
        res$results.compiled[nrow(res$results.compiled), 'pwr.g'] <-
          res$results.raw[[length(res$results.raw)]]$power$power;
        res$results.compiled[nrow(res$results.compiled), 'pwr.small'] <-
          res$results.raw[[length(res$results.raw)]]$power.small$power;
        res$results.compiled[nrow(res$results.compiled), 'pwr.medium'] <-
          res$results.raw[[length(res$results.raw)]]$power.medium$power;
        res$results.compiled[nrow(res$results.compiled), 'pwr.large'] <-
          res$results.raw[[length(res$results.raw)]]$power.large$power;
        
        res$results.compiled[nrow(res$results.compiled), 't'] <-
          res$results.raw[[length(res$results.raw)]]$t;
        res$results.compiled[nrow(res$results.compiled), 'df'] <-
          res$results.raw[[length(res$results.raw)]]$df;
        res$results.compiled[nrow(res$results.compiled), 'p'] <-
          res$results.raw[[length(res$results.raw)]]$p;
      }
      tempDat <- res$results.compiled[res$results.compiled$x==x[xVar], ];
      res$plots.compiled[[x[xVar]]] <- ggplot(tempDat,
                                              aes(x=y, y=g,
                                                  ymin=g.ci.lo,
                                                  ymax=g.ci.hi));
      if (tolower(orientation) == "vertical") {
        res$plots.compiled[[x[xVar]]] <- res$plots.compiled[[x[xVar]]] +
          coord_flip();
      }
      res$plots.compiled[[x[xVar]]] <- res$plots.compiled[[x[xVar]]] +
        geom_hline(yintercept=0, color=zeroLineColor, size=zeroLineSize) +
        geom_pointrange(size=1) + labs(x="interval variable",
                                         y="effect size g (unbiased estimate of Cohen's d)") +
        ggtitle(x[xVar]);;
    }
  }

  ### Set class & return result
  class(res) <- c("meanDiff.multi");
  return(res);
}

print.meanDiff.multi <- function (x, digits=x$digits, powerDigits=x$digits + 2, ...) {
  print(x$results.compiled, ...);
  print(x$plots.compiled, ...);
  invisible();
}
