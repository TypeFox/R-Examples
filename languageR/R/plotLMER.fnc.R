`plotLMER.fnc` <-
function(model, 
   xlabel = NA,     # label for X axis, if other than column name
   xlabs  = NA,     # vector of xlabels for multipanel plot
   ylabel = NA,     # label for Y axis, if other than dependent variable name
   ylimit = NA,     # ylimit to be set for plot or all subplots
   ilabel = NA,     # if not NA, then this will be the label for the interaction on side=3
   fun=NA,          # transform predicted values using function fun; if specified, 
                    # fun should be a function object, not the name (string) for the function
                    # for binomial models, fun is taken to be plogis by default
   pred=NA,         # if not NA, single plot will be produced for pred
   control = NA,    # if not NA,  list(predictor, val) for further control of the plot of the partial effect for pred
   ranefs = NA,     # if not NA,  list(group, level, predictor, predname2) for plotting curve for specific subject, item, ...
   n=100,           # number of points on X-axis for numerical predictors
   intr=NA,         # a list specifying an interaction, with as
                    # obligatory arguments the name of the interaction variable,
                    # followed by a vector of values for that variable, followed by
                    # the position for interaction labels ("beg", "mid", or "end", or NA, 
                    # in which case labels will not be shown),
                    # optionally followed by a list with a vector of colors and a vector of line types
   lockYlim=TRUE,   # lock ylim for all subplots to fit range of all predictors and
                    # ci's
   addlines=FALSE,  # for interaction plots for two factors, add lines to plot
   withList=FALSE,  # if true, a list with data frames for the individual panels is produced
   cexsize=0.5,     # the character expansion for text in plots for interactions
   linecolor=1,     # the color for the lines in the display, by default black
   addToExistingPlot = FALSE,  # if TRUE, plot will be added to previous plot, but only if pred is specified
   verbose = TRUE,  # provide information on stdout
   ...) {           # further graphical parameters




   ###############################################################################
   # validate the input
   ###############################################################################
   if (!is.na(xlabel[1])) {
     if (!is.character(xlabel)) 
       stop("xlabel should be a string\n")
   }
   if (!is.na(ylabel)) {
     if (!is.character(ylabel)) 
       stop("ylabel should be a string\n")
   }
   if (!is.na(ylimit[1])) {
     if ((!is.numeric(ylimit)) | (length(ylimit)!=2)) 
       stop("ylimit should be a two-element numeric vector\n")
   }
   if (!is.na(intr[1])) {
     if (!is.list(intr))
       stop("intr should be a list\n")
   }
   if (!is.numeric(n)) {
       stop("n should be an integer\n")
   }
   if (!is.na(pred)) {
     if (!is.character(pred)) 
       stop("pred should be a string\n")
   }
   if (!is.function(fun)) {
     if (!is.na(fun)) {
       stop("fun should be a function (not the name of a function)\n")
     }
   } 
   if ((length(grep("^glmer", as.character(model@call))) == 1) &
         (length(grep("binomial", as.character(model@call))) == 1)) {
       if (!is.function(fun)) {
         fun = plogis
         if (verbose == TRUE) cat("log odds are back-transformed to probabilities\n")
       }
   }

   if (is.na(pred)) addToExistingPlot = FALSE

   ###############################################################################
   # define variables for displaying interactions
   ###############################################################################

   # set some defaults
   conditioningPred = ""         # the name of the second predictor in the interaction
   conditioningVals = NULL       # the values of this second predictor, can be factor
                                 # or numeric
   conditioningPos = NA

   conditioningColors = 1        # optionally, the colors to be displayed, defaults to 1
   conditioningLines = 1
   if (!is.na(intr[[1]])) {
     conditioningPred = intr[[1]]
     conditioningVals = intr[[2]]
     conditioningPos = intr[[3]]
     if (length(intr)==4) {
       conditioningColors = intr[[4]][[1]]
       if (length(conditioningColors) != length(conditioningVals)) {
         stop("number of colors and number of conditioning values mismatch")
       }
       conditioningLines = intr[[4]][[2]]
       if (length(conditioningLines) != length(conditioningLines)) {
         stop("number of line types and number of conditioning values mismatch")
       }
     }
   }
   
   if (length(ylimit)>1) {
     lockYlim=FALSE              
     # user-specified ylimits override limits calculated from data
   }

   if (!is.na(control[[1]])) {
     if (!((length(control) == 2) & is.list(control))) {
       stop("control should be a two-element list\n")
     }
   }

   if (!is.na(ranefs[[1]])) {
       if (!((length(ranefs) == 4) & is.list(ranefs))) {
         stop("ranefs should be a four-element list\n")
       }
   }

   ###############################################################################
   # set label for Y-axis  (taken by default from model formula)
   ###############################################################################

   if (is.na(ylabel)) ylabel = as.character(eval(model@call[2]$formula))[2]

   ###############################################################################
   # check for multiple predictors for multiple subplots
   ###############################################################################

   if (is.na(pred)) {   # no predictor specified, so subplots for each predictor
                        # in the model
     predictors = colnames(model@frame)
     ranefnames=unique(names(ranef(model)))
     depvar = as.character(eval(model@call[2]$formula))[2]
     # predictors for random intercepts are also listed in the model frame, so
     # we have to restrict the predictors to those up to the first random effect name
     predictors = predictors[1: (which(predictors==ranefnames[1])-1)]
     predictors = predictors[!predictors %in% c(ranefnames,depvar)]
   } else {             # only one plot is produced, with pred as predictor
     predictors = pred
   }

   if (!is.na(xlabs[1])) {
     if (length(xlabs) != length(predictors)) {
       stop("number of labels in xlabs is not the same as the number of predictors\n")
     }
   }


   ###############################################################################
   # we now loop over the predictors and make the (sub)plot(s)
   ###############################################################################

   plots = list()

   for (i in 1:length(predictors)) {

     # cat("preparing panel for", predictors[i], "\n")
     if (length(predictors) == 1) xlabelShow = xlabel
     else {
       if (!is.na(xlabs[1])) {
         xlabelShow = xlabs
       } else {
         xlabelShow = NA
       }
     }
     # prepare label for X-axis if not specified by the user; for
     # multiple subplots take name of predictor 
     if (is.na(xlabel[1]) | length(predictors)>1) {
       xlabel = predictors[i]
     }
     
     # we now destinguish between plot with interaction (if) and simple (sub)plot (else)
     if ((length(predictors)==1) & (!is.null(conditioningVals))) {
       #------------------------------------------------
       # we first prepare the colors for the interaction
       #------------------------------------------------
       if (is.null(conditioningColors)) {
         colors = rep(1, length(conditioningVals))  # just black
         lineTypes = rep(1, length(conditioningVals)) # just solid lines
       } else {
         colors = conditioningColors
         lineTypes = conditioningLines
         if (length(colors) < length(conditioningVals)) {
           nc = (length(conditioningVals) %% length(colors))+1
           colors=rep(colors, nc)
         }
         if (length(lineTypes) < length(conditioningVals)) {
           nc = (length(conditioningLines) %% length(lineTypes))+1
           lineTypes=rep(lineTypes, nc)
         }
       }
       ################################################################################
       # create a model matrix with only main effects instantiated, with
       # medians for numeric predictors, and zeros for all factor columns
       # and interactions, except for the predictor (or factor level) in the interaction
       # given by val
       ###############################################################################
       val = conditioningVals[1]
       m = makeDefaultMatrix.fnc(model, n, conditioningPred, val, control)
       #----------------------------------------------------------------------------
       # note: the matrix is upgraded within preparePredictor.fnc for the interactions
       #----------------------------------------------------------------------------
       subplots = list()
       dfr = preparePredictor.fnc(predictors[i], model, m, ylabel, fun, 
         val, xlabel=xlabel, ranefs, lty=1, col=0, ...)
       subplots[[1]] = dfr
       if (verbose==TRUE) {
         cat("effect sizes (ranges) for the interaction of ", 
           predictors[i], " and ", conditioningPred, ":\n")
         cat("   ", conditioningPred, " = ", val, ": ", max(dfr$Y)-min(dfr$Y), "\n")
       }
       #----------------------------------------------------------------------------
       # finally loop over the different values of the interacting predictor, either
       # the factor levels or the numerical values (often quantiles)
       #----------------------------------------------------------------------------
       for (j in 2:length(conditioningVals)) {
         val = conditioningVals[j]
         m = makeDefaultMatrix.fnc(model, n, conditioningPred, val, control)
         dfr = preparePredictor.fnc(predictors[i], model, m, ylabel, fun, 
           val,  ranefs, lty=j, xlabel=xlabel, ...)
         subplots[[j]] = dfr
         if (verbose==TRUE) {
           cat("   ", conditioningPred, " = ", val, ": ", max(dfr$Y)-min(dfr$Y), "\n")
         }
       }
       plots[[i]] = subplots
     } else {  # a straightforward (sub)plot
       ################################################################################
       # first create a model matrix with only main effects instantiated, with
       # medians for numeric predictors, and zeros for all factor columns
       # and interactions
       ###############################################################################
       lineTypes = 1
       m = makeDefaultMatrix.fnc(model, n, "", NULL, control)
       # and make the plot
       dfr = preparePredictor.fnc(predictors[i], model, m, ylabel, fun, 
         val=NA, xlabel=xlabel, ranefs, ...)
       plots[[i]] = dfr
       if (verbose==TRUE) {
         cat("effect size (range) for ", predictors[i], "is ", max(dfr$Y)-min(dfr$Y), "\n")
       }
     }
   }
   names(plots) = predictors

   if (!is.na(ilabel)) {
      intrName = ilabel
   } else {
      intrName = conditioningPred
   }
   plotAll.fnc(plots, sameYrange=lockYlim, ylabel, xlabel = xlabelShow, intrName=intrName, 
   pos=conditioningPos, ylimit=ylimit, addlines=addlines, cexsize = cexsize, 
   conditioningVals=conditioningVals, conditioningColors = colors, conditioningLines=lineTypes, 
   lineColor=linecolor, addToExistingPlot, ...)

   if (withList) return(plots)

}

