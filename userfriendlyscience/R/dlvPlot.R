### XKCD styled plots
#require(xkcd)
#vignette("xkcd-intro")

### Histogram with dots
# http://stackoverflow.com/questions/16216312/how-to-plot-stacked-point-histograms-in-ggplot2-in-r


### Note: this is necessary to prevent Rcmd CHECK from throwing a note;
### otherwise it think these variables weren't defined yet.
utils::globalVariables(c("y_density", "yMaxFromY"));

### Theme used for the plots
dlvTheme <- function(base_size = 14, base_family = "", ...) {
  # Starts with theme_grey and then modify some parts
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text         = element_text(colour="#000000", size = rel(0.8)),
      axis.ticks        = element_line(colour = "black"),
      axis.title        = element_blank(),
      legend.text       = element_text(size = rel(0.6)),
      legend.key        = element_rect(colour = "grey80"),
      legend.position   = "top",
      legend.direction  = "horizontal",
      legend.key.size   = unit(6, "mm"),
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.border      = element_rect(fill = NA, colour = "grey50"),
      panel.grid.major  = element_line(colour = "grey90", size = 0.2),
      panel.grid.minor  = element_line(colour = "grey98", size = 0.5),
      strip.background  = element_rect(fill = "grey80", colour = "grey50"),
      strip.background  = element_rect(fill = "grey80", colour = "grey50"),
      panel.margin      = unit(c(.5), "cm"),
      ...
    )
}

dlvPlot <- function(dat, x = NULL, y, z = NULL, conf.level = .95,
                    jitter = "FALSE", binnedDots = TRUE, binwidth=NULL,
                    error="lines", dotsize="density", densityDotBaseSize=3,
                    normalDotBaseSize=1, violinAlpha = .2, dotAlpha = .4,
                    lineAlpha = 1, connectingLineAlpha = 1,
                    meanDotSize=5, posDodge=0.2, errorType = "both") {
  ### This function constructs a dot-line-violin plot.

  ### Create object to return results
  res <- list();
  
  ### Store data
  res$dat.raw <- dat;
  ### Remove irrelevant variables
  res$dat <- dat <- data.frame(dat[, c(x, y, z)]);
  
  ### Remove incomplete cases
  res$dat <- data.frame(dat[complete.cases(dat), ]);

  ### Replace names again
  names(dat) <- names(res$dat) <- c(x, y, z);
  
  if(!is.null(x) & !(is.factor(dat[, x]))) {
    warning("Error: variable x (', x,') is not of type factor. X must be a categorical ",
            "variable with a limited number of categories. If this is the case, but it's ",
            "simply stored as a numerical vector, use the 'factor' function to convert ",
            "it (see '?factor'). Trying to convert x myself now.");
    res$dat[[x]] <- factor(res$dat[[x]]);
  }

  if(!is.null(z) & !(is.factor(dat[, z]))) {
    warning("Error: variable z (', z,') is not of type factor. Z must be a categorical ",
            "variable with a limited number of categories. If this is the case, but it's ",
            "simply stored as a numerical vector, use the 'factor' function to convert ",
            "it (see '?factor'). Trying to convert z myself now.");
    res$dat[[z]] <- factor(res$dat[[z]]);
  }
  
  if(is.null(x)) {
    ### We have no predictor variable - this means we construct univariate plots.

    ### Now check whether we have to construct one or several.
    if(length(y)==1) {
      
      ###############################################################
      ### Constructing one univariate plot                        ###
      ###############################################################
      
      ### Store variable name in dataframe
      if (is.null(res$dat$variable)) {
        res$dat$variable <- y;
        xVarName <- 'variable';
      }
      else {
        res$dat$variable_dlvPlot <- y;
        xVarName <- 'variable_dlvPlot';
      }
      
      ### Store density at y value
      dens <- density(res$dat[[y]], na.rm=TRUE);
      res$dat$y_density <- approx(dens$x, dens$y, xout=res$dat[[y]])$y;
      ### Multiply so that points at average density have size 1
      res$dat$y_density <- res$dat$y_density *
        (densityDotBaseSize/mean(res$dat$y_density, na.rm=TRUE));
      
      ### Construct dataframe with confidence interval info
      n <- nrow(res$dat);
      mean <- mean(res$dat[, y]);
      sd <- sd(res$dat[, y]);
      se <- sd / sqrt(nrow(res$dat));
      criticalValue <- qt(1-((1-conf.level)/2), df=n-1); 
      ci.lo <- mean - criticalValue * se;
      ci.hi <- mean + criticalValue * se;
      res$descr <- data.frame(y = y,
                              n = n,
                              mean = mean, sd = sd,
                              se = se,
                              ci.lo = ci.lo,
                              ci.hi = ci.hi);
      res$yRange=c(min(res$dat[[y]][!is.na(res$dat[[y]])]),
                   max(res$dat[[y]][!is.na(res$dat[[y]])]));
      
      ### Generate plot
      res$plot <- ggplot(data=res$dat, aes_string(x=xVarName, y=y));
      res$plot <- res$plot + dlvTheme();
      res$plot <- res$plot + geom_violin(trim=FALSE, alpha=violinAlpha, fill="#BBBBBB", linetype="blank");
      if (jitter) {
        res$plot <- res$plot + geom_jitter(position=position_jitter(width=.1, height=.01), alpha=dotAlpha);
      }
      else {
        if (binnedDots) {
          tempBinwidth <- ifelse(is.null(binwidth), (res$yRange[2]-res$yRange[1])/30, binwidth);
          res$plot <- res$plot + geom_dotplot(alpha=dotAlpha, show.legend=FALSE,
                                              binaxis="y", binwidth=tempBinwidth, dotsize=normalDotBaseSize,
                                              stackdir="center", position=position_dodge(width=posDodge));
        }
        else if (dotsize=="density") {
          res$plot <- res$plot + geom_point(aes(size=y_density), color='grey60',
                                            alpha=dotAlpha, show.legend=FALSE);
        }
        else {
          res$plot <- res$plot + geom_point(alpha=dotAlpha, dotsize=normalDotBaseSize);
        }
      }
      if (error == "lines") {
        if (errorType=="ci") {
          res$plot <- res$plot + geom_pointrange(data=res$descr,
                                                 aes(x=y, y=mean, ymin=ci.lo, ymax=ci.hi),
                                                 size = 1, alpha=lineAlpha);
        } else if (errorType=="se") {
          res$plot <- res$plot + geom_pointrange(data=res$descr,
                                                 aes(x=y, y=mean, ymin=mean-se, ymax=mean+se),
                                                 size = 1, alpha=lineAlpha);
        } else if (errorType=="both") {
          res$plot <- res$plot + geom_pointrange(data=res$descr,
                                                 aes(x=y, y=mean, ymin=ci.lo, ymax=ci.hi),
                                                 size = 1, alpha=lineAlpha);
          res$plot <- res$plot + geom_errorbar(data=res$descr,
                                               aes(x=y, y=mean, ymin=mean-se, ymax=mean+se),
                                               size = 2, alpha=lineAlpha, width=0);
        }        
      }
      else if (error == "whiskers") {
        res$plot <- res$plot + geom_errorbar(data=res$descr,
                                             aes(x=y, y=mean, ymin=ci.lo, ymax=ci.hi),
                                             size = 1, width=.1, alpha=lineAlpha);
      }
      res$plot <- res$plot + geom_point(data=res$descr,
                                        aes(x=y, y=mean), size=meanDotSize, alpha=lineAlpha);
      
    }
    else {
      
      ###############################################################
      ### Constructing several univariate plots                   ###
      ###############################################################
      
      ### Apparently, we have to construct several plots.
      ### First generate a dataframe where the variables names
      ### are stored in another variable that we can use to
      ### make categories on the x axis
      
      ### Store original dataframe
      res$dat.original <- res$dat;
      res$dat <- data.frame();
      ### Create empty descriptives dataframe
      res$descr <- data.frame();
      
      ### Loop through original dataframe and construct new one
      for (currentVar in y) {
        tempDf <- data.frame(y = res$dat.original[, currentVar]);
        tempDf$x <-  currentVar;
        ### Store density for at y value
        dens <- density(tempDf$y, na.rm=TRUE);
        tempDf$y_density <- approx(dens$x, dens$y, xout=tempDf$y)$y;
        tempDf$y_density <- tempDf$y_density *
          (densityDotBaseSize/mean(tempDf$y_density, na.rm=TRUE));
        ### Store y values and name of y variable in res$dat dataframe
        res$dat <- rbind(res$dat, tempDf);
        ### Get mean and confidence interval for descriptives table
        n <- nrow(tempDf);
        mean <- mean(tempDf$y);
        sd <- sd(tempDf$y);
        se <- sd / sqrt(nrow(tempDf));
        criticalValue <- qt(1-((1-conf.level)/2), df=n-1); 
        ci.lo <- mean - criticalValue * se;
        ci.hi <- mean + criticalValue * se;
        ### Add descriptives
        res$descr <- rbind(res$descr, data.frame(y = currentVar,
                                                 n = n,
                                                 mean = mean, sd = sd,
                                                 se = se,
                                                 ci.lo = ci.lo,
                                                 ci.hi = ci.hi));
      }
      
      res$yRange=c(min(res$dat[['y']][!is.na(res$dat[['y']])]),
                   max(res$dat[['y']][!is.na(res$dat[['y']])]));
      
      res$plot <- ggplot(data=res$dat, aes(x=x, y=y));
      res$plot <- res$plot + dlvTheme();
      res$plot <- res$plot + geom_violin(trim=FALSE, alpha=violinAlpha, fill="#BBBBBB", linetype="blank");
      if (jitter) {
        res$plot <- res$plot + geom_jitter(position=position_jitter(width=.1, height=.01), alpha=dotAlpha);
      }
      else {
        if (binnedDots) {
          tempBinwidth <- ifelse(is.null(binwidth), (res$yRange[2]-res$yRange[1])/30, binwidth);
          res$plot <- res$plot + geom_dotplot(alpha=dotAlpha, show.legend=FALSE,
                                              binaxis="y", binwidth=tempBinwidth, dotsize=normalDotBaseSize,
                                              stackdir="center", position=position_dodge(width=posDodge)
                                              );
        }
        else if (dotsize=="density") {
          res$plot <- res$plot + geom_point(aes(size=y_density), color='grey60',
                                            alpha=dotAlpha, show.legend=FALSE);
        }
        else {
          res$plot <- res$plot + geom_point(alpha=dotAlpha, dotsize=normalDotBaseSize);
        }
      }
      if (error == "lines") {
        if (errorType=="ci") {
          res$plot <- res$plot + geom_pointrange(data=res$descr,
                                                 aes(x=x, y=mean, ymin=ci.lo, ymax=ci.hi),
                                                 size = 1, alpha=lineAlpha);
        } else if (errorType=="se") {
          res$plot <- res$plot + geom_pointrange(data=res$descr,
                                                 aes(x=x, y=mean, ymin=mean-se, ymax=mean+se),
                                                 size = 1, alpha=lineAlpha);
        } else if (errorType=="both") {
          res$plot <- res$plot + geom_pointrange(data=res$descr,
                                                 aes(x=x, y=mean, ymin=ci.lo, ymax=ci.hi),
                                                 size = 1, alpha=lineAlpha);
          res$plot <- res$plot + geom_errorbar(data=res$descr,
                                               aes(x=x, y=mean, ymin=mean-se, ymax=mean+se),
                                               size = 2, alpha=lineAlpha, width=0);
        }
      }
      else if (error == "whiskers") {
        res$plot <- res$plot + geom_errorbar(data=res$descr,
                                             aes(x=y, y=mean, ymin=ci.lo, ymax=ci.hi),
                                             size = 1, width=.1, alpha=lineAlpha);
      }
      res$plot <- res$plot + geom_point(data=res$descr,
                                        aes(x=y, y=mean), size=meanDotSize, alpha=lineAlpha);
      
    }
  }
  else {
    ### We have a predictor variable, so check whether we have a moderator
    if (is.null(z)) {
      
      ###############################################################
      ### Constructing multivariate plot without moderator        ###
      ###############################################################
      
      ### Construct dataframe with confidence interval info
      res$descr <- ddply(.data = res$dat, .variables = c(x),
                         .fun = function (dat, conf.level) {
                           dat <- dat[complete.cases(dat), ];
                           n <- nrow(dat);
                           mean <- mean(dat[, y]);
                           sd <- sd(dat[, y]);
                           se <- sd / sqrt(nrow(dat));
                           criticalValue <- qt(1-((1-conf.level)/2), df=n-1); 
                           ci.lo <- mean - criticalValue * se;
                           ci.hi <- mean + criticalValue * se;
                           rslt <- data.frame(x = dat[1, x],
                                              y = y,
                                              n = nrow(dat),
                                              mean = mean, sd = sd,
                                              se = se, ci.lo = ci.lo,
                                              ci.hi = ci.hi);
                           rslt <- rslt[complete.cases(rslt), ];
                           return(rslt);
                         }, conf.level=conf.level);
      ### Store densities; must be done for each group (value of x)
      ### separately
      res$dat <- ddply(.data = res$dat, .variables = c(x),
                       .fun = function (dat) {
                         ### Store density for at y value
                         dens <- density(dat[[y]], na.rm=TRUE);
                         dat$y_density <- approx(dens$x, dens$y, xout=dat[[y]])$y;
                         ### Multiply with densityDotBaseSize / mean (this allows
                         ### control over the size of the dots)
                         dat$y_density <- dat$y_density *
                           (densityDotBaseSize/mean(dat$y_density, na.rm=TRUE));
                         return(dat);
                       });
      
      res$yRange=c(min(res$dat[[y]][!is.na(res$dat[[y]])]),
                   max(res$dat[[y]][!is.na(res$dat[[y]])]));
      
      res$plot <- ggplot(data=res$dat, aes_string(x=x, y=y));
      res$plot <- res$plot + dlvTheme();
      res$plot <- res$plot + geom_violin(trim=FALSE, alpha=violinAlpha, fill="#BBBBBB", linetype="blank", position=position_dodge(width=posDodge));
      if (jitter) {
        res$plot <- res$plot + geom_jitter(position=position_jitter(width=.1, height=.01), alpha=dotAlpha);
      }
      else {
        if (binnedDots) {
          tempBinwidth <- ifelse(is.null(binwidth), (res$yRange[2]-res$yRange[1])/30, binwidth);
          res$plot <- res$plot + geom_dotplot(alpha=dotAlpha, show.legend=FALSE,
                                              binaxis="y", binwidth=tempBinwidth, dotsize=normalDotBaseSize,
                                              stackdir="center", position=position_dodge(width=posDodge));
        }
        else if (dotsize=="density") {
          res$plot <- res$plot + geom_point(aes(size=y_density),
                                            alpha=dotAlpha, show.legend=FALSE);
        }
        else {
          res$plot <- res$plot + geom_point(alpha=dotAlpha, dotsize=normalDotBaseSize);
        }
      }
      if (error == "lines") {
        if (errorType=="ci") {
        res$plot <- res$plot + geom_pointrange(data=res$descr,
                                               aes(x=x, y=mean, ymin=ci.lo, ymax=ci.hi),
                                               size = 1, alpha=lineAlpha);
        } else if (errorType=="se") {
          res$plot <- res$plot + geom_pointrange(data=res$descr,
                                                 aes(x=x, y=mean, ymin=mean-se, ymax=mean+se),
                                                 size = 1, alpha=lineAlpha);
        } else if (errorType=="both") {
          res$plot <- res$plot + geom_pointrange(data=res$descr,
                                                 aes(x=x, y=mean, ymin=ci.lo, ymax=ci.hi),
                                                 size = 1, alpha=lineAlpha);
          res$plot <- res$plot + geom_errorbar(data=res$descr,
                                               aes(x=x, y=mean, ymin=mean-se, ymax=mean+se),
                                               size = 2, alpha=lineAlpha, width=0);
        }
      }
      else if (error == "whiskers") {
        res$plot <- res$plot + geom_errorbar(data=res$descr,
                                             aes(x=x, y=mean, ymin=ci.lo, ymax=ci.hi),
                                             size = 1, width=.1, alpha=lineAlpha);
      }
      res$plot <- res$plot + stat_summary(fun.y=mean, geom="point", size=meanDotSize, alpha=lineAlpha);
      res$plot <- res$plot + geom_line(data=res$descr,
                                       aes(x=as.numeric(x), y=mean), size=1, alpha=connectingLineAlpha);
    }
    else {
      
      ###############################################################
      ### Constructing multivariate plot with moderator           ###
      ###############################################################
      
      ### Construct dataframe with confidence interval info
      res$descr <- ddply(.data = res$dat, .variables = c(x, z),
                     .fun = function (dat, conf.level) {
                       dat <- dat[complete.cases(dat), ];
                       n <- nrow(dat);
                       mean <- mean(dat[, y]);
                       sd <- sd(dat[, y]);
                       se <- sd / sqrt(nrow(dat));
                       criticalValue <- qt(1-((1-conf.level)/2), df=n-1); 
                       ci.lo <- mean - criticalValue * se;
                       ci.hi <- mean + criticalValue * se;
                       res <- data.frame(x = dat[1, x],
                                         y = y,
                                         z = dat[1, z],
                                         n = nrow(dat),
                                         mean = mean, sd = sd,
                                         se = se, ci.lo = ci.lo,
                                         ci.hi = ci.hi);
                       return(res[complete.cases(res), ]);
                     }, conf.level=conf.level);
      ### Store densities; must be done for each group (value of x)
      ### separately
      res$dat <- ddply(.data = res$dat, .variables = c(x, z),
                       .fun = function (dat) {
                         ### Store density for at y value
                         dens <- density(dat[[y]], na.rm=TRUE);
                         dat$y_density <- approx(dens$x, dens$y, xout=dat[[y]])$y;
                         ### Multiply with densityDotBaseSize / mean (this allows
                         ### control over the size of the dots)
                         dat$y_density <- dat$y_density *
                           (densityDotBaseSize/mean(dat$y_density, na.rm=TRUE));
                         return(dat);
                       });

      res$yRange=c(min(res$dat[[y]][!is.na(res$dat[[y]])]),
                   max(res$dat[[y]][!is.na(res$dat[[y]])]));
      
      res$plot <- ggplot(data=res$dat, aes_string(x=x, y=y, z=z, colour=z, group=paste0(x,":",z)));
      res$plot <- res$plot + dlvTheme();
      res$plot <- res$plot + geom_violin(data=res$dat, aes_string(fill=z), alpha=violinAlpha, trim=FALSE, linetype="blank", position=position_dodge(width=posDodge));
      if (jitter) {
        res$plot <- res$plot + geom_jitter(position=position_jitter(width=.1, height=.01), alpha=dotAlpha);
      }
      else {
        if (binnedDots) {
          tempBinwidth <- ifelse(is.null(binwidth), (res$yRange[2]-res$yRange[1])/30, binwidth);
          res$plot <- res$plot + geom_dotplot(alpha=dotAlpha, show.legend=FALSE,
                                              aes_string(fill=z), binaxis="y",
                                              binwidth=tempBinwidth, dotsize=normalDotBaseSize,
                                              stackdir="center", position=position_dodge(width=posDodge));
        }
        else if (dotsize=="density") {
          res$plot <- res$plot + geom_point(aes(size=y_density),
                                            alpha=dotAlpha, show.legend=FALSE, position=position_dodge(width=posDodge));
        }
        else {
          res$plot <- res$plot + geom_point(alpha=dotAlpha, dotsize=normalDotBaseSize, position=position_dodge(width=posDodge));
        }
      }
      if (error == "lines") {
        if (errorType=="ci") {
          res$plot <- res$plot + geom_pointrange(data=res$descr,
                                                 aes(x=x, y=mean, ymin=ci.lo, ymax=ci.hi, group=z),
                                                 size = 1, alpha=lineAlpha, position=position_dodge(width=posDodge));
        } else if (errorType=="se") {
          res$plot <- res$plot + geom_pointrange(data=res$descr,
                                                 aes(x=x, y=mean, ymin=mean-se, ymax=mean+se, group=z),
                                                 size = 1, alpha=lineAlpha, position=position_dodge(width=posDodge));
        } else if (errorType=="both") {
          res$plot <- res$plot + geom_pointrange(data=res$descr,
                                                 aes(x=x, y=mean, ymin=ci.lo, ymax=ci.hi, group=z),
                                                 size = 1, alpha=lineAlpha, position=position_dodge(width=posDodge));
          res$plot <- res$plot + geom_errorbar(data=res$descr,
                                               aes(x=x, y=mean, ymin=mean-se, ymax=mean+se, group=z),
                                               size = 2, alpha=lineAlpha, width=0, position=position_dodge(width=posDodge));
        }
      }
      else if (error == "whiskers") {
        res$plot <- res$plot + geom_errorbar(data=res$descr,
                                             aes(x=x, y=mean, ymin=ci.lo, ymax=ci.hi, group=z),
                                             size = 1, width=.1, alpha=lineAlpha, position=position_dodge(width=posDodge));
      }
      res$plot <- res$plot + stat_summary(fun.y=mean, geom="point", size=meanDotSize, position=position_dodge(width=posDodge));
      res$plot <- res$plot + geom_line(data=res$descr, aes(x=x, y=mean, group=z), size=1,
                                       alpha=connectingLineAlpha, position=position_dodge(width=posDodge));
    }
  }
  
  assign('yMaxFromY', max(res$plot$data[, res$plot$labels$y]), envir = res$plot$plot_env);
  res$plot <- res$plot + aes(ymax = yMaxFromY);
  
  ### Set class of result
  class(res) <- c('dlvPlot');
  
  ### Return result
  return(res);
}

print.dlvPlot <- function(x, ...) {
  print(x$plot, ...);
  invisible();
}