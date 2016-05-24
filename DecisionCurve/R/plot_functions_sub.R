

preparePlotData   <- function(x, curve.names, confidence.intervals){
#  if(missing(curve.names)) curve.nams <- NA
  #extract data from x,
  #there is only one DecisionCurve to plot
  if(class(x) == "decision_curve"){

    dc.data <- x$derived.data

    if(is.na(confidence.intervals[1])){
      confidence.intervals <- x$confidence.intervals
    }else{
      confidence.intervals <- ifelse(confidence.intervals, 1, "none")
    }

    if(is.na(curve.names[1])){
      predictors <- unique(dc.data$model)
      predictors <- predictors[!is.element(predictors, c("None", "All"))]
    }else{
      dc.data$model[!is.element(dc.data$model, c("None", "All"))] <- curve.names[1]
      predictors <- curve.names[1]
    }

  }else if(class(x)=="list"){
    xx <- NULL
    #check to make sure each element of the list is a decision_curve object,
    # or else we get funky results.
    if(!all(sapply(x, FUN = function(xx) class(xx) == "decision_curve")) ){
      stop("One or more elements of the list provided is not an object of class 'decision_curve' (output from the function 'decision_curve').")
    }

    if(is.na(confidence.intervals[1])){
      confidence.intervals <- x[[1]]$confidence.intervals
    }else{
      confidence.intervals <- ifelse(confidence.intervals, 1, "none")
    }

    message("Note: When multiple decision curves are plotted, decision curves for 'All' are calculated using the first DecisionCurve in the list provided.")

    model <- NULL #appease check
    #multiple dc's to plot
    #pull the "all' and 'none' curves from the first element in x
    dc.data <- subset(x[[1]]$derived.data, is.element(model, c("All", "None")))

    #fill in ci variables for if confidence intervals weren't calculated using DecisionCurve
    if(ncol(dc.data) == 9 ) dc.data <- add.ci.columns(dc.data)

    predictors <- NULL

    #loop through the remaining curves
    i = 0
    for(ll in x){
      i = i + 1
      #extract data to add
      newdata <-  subset(ll$derived.data, !is.element(model, c("All", "None")))
      #predictor name

      if(is.na(curve.names[1])){
        #check to make sure the name is different
        newpred <- unique(newdata$model)
        if(is.element(newpred, predictors)) stop("After extracting the curve names from the decision_curve object, the names of the decision curves provided are the same for two or more decision_curve objects. Please set curve.names to avoid errors in plotting.")
      }else{
        newdata$model <- curve.names[i]
        newpred <- unique(newdata$model)
      }

      predictors <- c(predictors, newpred)

      #if confidence intervals weren't calculated
      if(ncol(newdata) == 9 ){
        if(confidence.intervals) warning(paste("confidence interval plotting were requested for curve '", newpred, "' but not calculated using decision_curve", sep = ''))
        #fill in ci variables for if confidence intervals weren't calculated using DecisionCurve
        newdata <- add.ci.columns(newdata)
      }

      dc.data <- rbind(dc.data, newdata)
    }


  }

  return(list(dc.data = dc.data,
              predictors = predictors,
              confidence.intervals = confidence.intervals))

}

########################################################################################
########################################################################################
########################################################################################
########################################################################################



plot_generic<- function(xx, predictors, value, plotNew,
                        standardize, confidence.intervals,
                        cost.benefit.axis = TRUE, cost.benefits, n.cost.benefits,
                        cost.benefit.xlab, xlab, ylab,
                        col, lty, lwd,
                        xlim, ylim, legend.position,
                        lty.fpr = 2, lty.tpr = 1,
                        tpr.fpr.legend = FALSE,
                        impact.legend = FALSE,
                        population.size = 1000, ...){
## xx is output from get_DecisionCurve,
## others are directly from the function call

  #save old par parameters and reset them once the function exits.
  old.par<- par("mar"); on.exit(par(mar = old.par))


  xx.wide <- reshape::cast(xx, thresholds~model, value =  value, add.missing = TRUE, fill = NA)
  xx.wide$thresholds <- as.numeric(as.character(xx.wide$thresholds))

  if(is.numeric(confidence.intervals)){

    val_lower <- paste(value, "lower", sep = "_")
    val_upper <- paste(value, "upper", sep = "_")

    xx.lower <- cast(xx, thresholds~model, value = val_lower, add.missing = TRUE, fill = NA)
    xx.upper <- cast(xx, thresholds~model, value = val_upper, add.missing = TRUE, fill = NA)
    xx.lower$thresholds <- as.numeric(as.character(xx.lower$thresholds))
    xx.upper$thresholds <- as.numeric(as.character(xx.upper$thresholds))
  }


  # adjust margins to add extra x-axis
  if(cost.benefit.axis) par(mar = c(7.5, 4, 3, 2) + 0.1)

  #set default ylim if not provided


  #initial call to plot and add gridlines
  if(plotNew){

  plot(xx.wide$thresholds, xx.wide$None, type = "n", ylim = ylim,
       col = "black", xlim = xlim,  xlab = "", ylab = ylab, frame.plot = FALSE, ...)

  grid(lty = 1, col = "grey92")
  }

  if(is.element(value, c("NB", "sNB"))){
  #plot none and all
  lines(xx.wide$thresholds, xx.wide$None, type = "l",
        col = col[length(predictors)+ 2],
        lty = lty[length(predictors)+ 2],
        lwd = lwd[length(predictors)+ 2])

  lines(xx.wide$threshold, xx.wide$All, type = "l",
        col = col[length(predictors)+ 1],
        lty = lty[length(predictors)+ 1],
        lwd = lwd[length(predictors)+ 1])

  if(is.numeric(confidence.intervals)){
    lines(xx.lower[,c("thresholds", "All")],
          col = col[length(predictors)+ 1],
          lty = lty[length(predictors)+ 1],
          lwd = lwd[length(predictors)+ 1]/2)

    lines(xx.upper[,c("thresholds", "All")],
          col = col[length(predictors)+ 1],
          lty = lty[length(predictors)+ 1],
          lwd = lwd[length(predictors)+ 1]/2)

  }
  }

  #the clinical impact plots are on a different scale
  if(is.element(value, c("DP" , "prob.high.risk"))){
    #population.size
    ps <- population.size
  }else{
    ps <- 1
  }

  #plot each predictor
  for(i in 1:length(predictors)){
    #plot ci's if asked for

    j <- ifelse(is.element(value, c("TPR", "prob.high.risk")), 1, i)
    j <- ifelse(is.element(value, c("FPR", "DP")), 2, i)
    if(is.numeric(confidence.intervals)){
      #get rid of cases missing for that predictor; this sometimes
       #happens due to different thresholds for each predictor
      cc <- complete.cases(xx.lower[,c("thresholds", predictors[i])])

      lines(x = xx.lower[cc, c("thresholds")],
            y = xx.lower[cc, c(predictors[i])]*ps,
            type = "l",  col = col[j], lty = lty[i], lwd = lwd[i]/2)

      cc <- complete.cases(xx.upper[,c("thresholds", predictors[i])])

      lines(x = xx.upper[cc, c("thresholds")],
            y = xx.upper[cc, c(predictors[i])]*ps,
            type = "l",  col = col[j], lty = lty[i], lwd = lwd[i]/2)

    }
     cc <- complete.cases(xx.wide[,c("thresholds", predictors[i])])

     lines(x = xx.wide[cc, c("thresholds")],
           y = xx.wide[cc, c(predictors[i])]*ps,
           type = "l",  col = col[j], lty = lty[i], lwd = lwd[i])


  }

  #add legend
  if(is.element(legend.position, c("bottomright", "topright", "bottomleft", "topleft", "right", "left", "top", "bottom"))){

    if(value == "NB" | value == "sNB"){
     legend(legend.position, lty = lty, col = col, lwd = lwd, legend = c(predictors, "All", "None"))
    }else if(tpr.fpr.legend){
      n.preds <- length(predictors)
      legend(legend.position,
             lty = c( lty.tpr, lty.fpr),
             col = col,
             lwd = lwd, legend = c("True positive rate", "False positive rate"))


    } else if(impact.legend){
      legend(legend.position,
             lty = c( 1, 2),
             col = col,
             lwd = lwd, legend = c("Number high risk", "Number high risk with outcome"))

    }

  }

  #add cost benefit axis if wanted
  if(cost.benefit.axis){
    tmp <- Add_CostBenefit_Axis(xlim = xlim,
                                cost.benefits = cost.benefits,
                                n.cost.benefits = n.cost.benefits,
                                line = 4)
    mtext(xlab, 1, 2.2)
    mtext(cost.benefit.xlab, side = 1, 6.1)
  }else{
    mtext(xlab, side = 1, 3)
  }
}

####################################################################################
####################################################################################
####################################################################################
####################################################################################

