# ---------------------------------------------------------------------------
#
# Abbreviating items (from questionnaire or other) measures 
# using genetic algorithms (GAs)
#
# version 1.1 Feb 2016
#
# Luca Scrucca
# University of Perugia, Italy
# luca@stat.unipg.it
#
# ---------------------------------------------------------------------------
#
# Based on the code gaabbreviate.R v 0.11 (2/10/2011) of Tal Yarkoni 
# <tyarkoni@gmail.com>
# http://www.shortermeasures.com/
# http://www.talyarkoni.org/blog/2010/03/31/abbreviating-personality-measures-in-r-a-tutorial/
#
# ---------------------------------------------------------------------------

GAabbreviate <- function(items = NULL, scales = NULL, 
                         itemCost = 0.05, maxItems = 5, 
                         maxiter = 100, popSize = 50, ...,
                         plot = FALSE, verbose = TRUE, 
                         crossVal = TRUE, impute = FALSE, 
                         pairwise = FALSE, minR = 0, 
                         sWeights = NULL, nSample = NULL) 
{
  
  if(is.null(items) & is.null(scales)) 
    stop("Invalid input provided! Either the items and scales, or a 'GAabbreviate' object  arguments must be provided.")

  if(!is.numeric(maxItems) || length(maxItems) > 1 || maxItems < 1) 
    stop("maxItems must be a positive integer!")
      
  items <- data.matrix(items)
  scales <- data.matrix(scales)
  nSubjects <- nrow(items)
  nItems <- ncol(items)
  nScales <- ncol(scales)
  if(nSubjects != nrow(scales)) 
    stop("Number of subjects in item and scale matrices do not match!")
  
  # Equal weighting of scales by default
  if(is.null(sWeights)) sWeights <- rep(1, nScales)
  if(length(sWeights) != nScales) 
    stop("Length of scale weights vector must match the number of scales!")
  
  # Can sample only a subset of observations
  if(!is.null(nSample) && is.numeric(nSample) && nSample > 0) 
    { if(verbose)
         message(paste("Using", nSample, "of", nSubjects, 
                       "subjects to generate measure..."))
         samp <- sample(nSubjects, nSample, replace = FALSE)  
         nSubjects <- nSample
         items <- items[samp,,drop=FALSE]
         scales <- scales[samp,,drop=FALSE]
  }
    
  # If cross-validation is TRUE, randomly split sample in two halves
  trainSubjects <- 1:nSubjects
  if(crossVal) 
    { ntrain <- round(nSubjects/2)
      trainSubjects <- sample(1:nSubjects, ntrain)
  }
  crossvalSubjects <- setdiff(seq(nSubjects), trainSubjects)

  # Check for missing values 
  impute <- as.logical(impute)
  if(any(is.missing(items)) || any(is.missing(scales)))
    { if(!impute)
        { message("There seem to be missing values in either 'items' or 'scales' matrix. Please ensure you impute missing values in items and create scale scores using a method most appropriate for these data before running GAabbreviate(). Note that using 'impute = TRUE' argument in GAabbreviate() would use mean imputation, which may or may not be appropriate for these data.")
          return() }
      # Impute the mean for missing values
      if(any(is.missing(items)))  items <- impute(items)
      if(any(is.missing(scales))) scales <- impute(scales)
  }
        
  # Create new 'GAabbreviate' object if none given
  obj <- list(data = list(items = items,
                          scales = scales,
                          nItems = nItems,
                          nScales = nScales,
                          trainSubjects = trainSubjects,
                          crossvalSubjects = crossvalSubjects),
              settings = list(maxiter = maxiter,
                              itemCost = itemCost,
                              maxItems = maxItems,
                              minR = minR,
                              popSize = popSize,
                              crossVal = crossVal,
                              impute = impute,
                              pairwise = pairwise,
                              sWeights = sWeights,
                              verbose = verbose),
              results = list(solution = NULL, # the solution
                             nItems = NULL,   # number of items in measure
                             cost = NULL,     # total cost of measure
                             meanR2 = NULL    # mean R^2 of measure
                            ),        
              # properties of the single best measure so far.
              # mainly used for display purposes in monitoring function.
              best = list(cost = Inf, 
                          fit = NULL),    # R^2s for all scales
              GA = NULL,
              measure = NULL
  )
  class(obj) <- "GAabbreviate"

  # env <- asNamespace("GAabbreviate")
  env <- globalenv()
  # env <- parent.frame()
  assign("TMPGAabbreviateObject", obj, envir = env)
  
  # Run GA and store in gaa
  if(verbose & interactive())
    cat("Starting GA run...\n")

  obj$GA <- ga(type = "binary", 
               fitness = fitness.GAabbreviate, 
               nBits = obj$data$nItems, 
               popSize = obj$settings$popSize,
               maxiter = obj$settings$maxiter, 
               keepBest = TRUE,
               monitor = function(...) 
                         monitor.GAabbreviate(..., 
                                              verbose = verbose, 
                                              plot = plot),
               ...)

  # get intermediate results saved at each iteration
  obj$results <- get("TMPGAabbreviateObject", envir = env)$results
  obj$best <- get("TMPGAabbreviateObject", envir = env)$best
  # generate measures
  obj$measure <- measure.GAabbreviate(obj)

  # remove object from global environment
  assign("TMPGAabbreviateObject", NULL, envir = env)
  # and return the object from function call 
  return(obj)
}

# Fitness function
fitness.GAabbreviate <- function(x) 
{
  # get values from object stored in package's environment
  # env <- asNamespace("GAabbreviate")
  env <- globalenv()
  # env <- parent.frame()
  obj <- get("TMPGAabbreviateObject", envir = env)
  items <- obj$data$items
  scales <- obj$data$scales
  nScales <- obj$data$nScales
  itemCost <- obj$settings$itemCost
  minR <- obj$settings$minR
  maxItems <- obj$settings$maxItems
  sWeights <- obj$settings$sWeights

  # Select items to include in measure
  active <- which(x == 1)
  activeItems <- items[,active,drop=FALSE]
  # Calculate item cost
  totalItemCost <- length(active) * itemCost
  # Calculate variance cost
  totalVarianceCost <- 0
  use <- ifelse(obj$settings$pairwise, "pairwise.complete.obs", "everything")
  cors <- cor(activeItems, scales, use = use) # NAs will crash in 'e' mode
  signs <- sign(cors)
  abscors <- abs(cors)
  varCosts <- rep(1,nScales)
  for(i in 1:nScales) 
     { # combine criteria: use minR as threshold, then rank-order items
       # and select the top absolute correlations, up to maxItems.
       cols <- order(abscors[,i], decreasing=TRUE)
       npass <- length(which(abscors[cols,i] >= minR))
       cols <- cols[1:min(npass, maxItems)]
       # If no items are left, we're totally unable to predict scale
       if(length(cols) == 0) next
       y <- scales[,i]
       x <- activeItems[,cols,drop=FALSE] %*% signs[cols,i]
       varCosts[i] <- (1-cor(x,y, use=use)^2)
  }
  
  totalCost <- totalItemCost + sum(sWeights*varCosts)
  attr(totalCost, "itemCost") <- totalItemCost
  attr(totalCost, "varCosts") <- varCosts
  
  # we want to minimize fitness, but GA maximize it, 
  # so returns negative total cost
  return(-totalCost)
}

# Monitoring function

monitor.GAabbreviate <- function(object, verbose = TRUE, plot = FALSE, 
                                 digits = getOption("digits"), ...)
{
  fitness <- -1*na.exclude(object@fitness)
  iter <- object@iter 
  
  if(verbose & interactive())
    { cat(paste("Iter =", iter, 
                " | Mean =", format(mean(fitness), digits = digits), 
                " | Best =", format(min(fitness), digits = digits),  "\n")) 
  }
  
  # get results at current iteration 
  # env <- asNamespace("GAabbreviate")
  env <- globalenv()
  # env <- parent.frame()
  gaa <- get("TMPGAabbreviateObject", envir = env)
  gaa$results$solution <- rbind(gaa$results$solution, 
                                object@bestSol[[iter]][1,])
  nItems <- sum(object@bestSol[[iter]][1,])
  gaa$results$nItems <- c(gaa$results$nItems, nItems)  
  gaa$results$cost <- c(gaa$results$cost, min(fitness))
  fitAttr <- attributes(fitness.GAabbreviate(gaa$results$solution[iter,]))
  gaa$best$cost <- min(fitness)
  gaa$best$fit <- fitAttr$varCosts
  meanR2 <- 1 - (gaa$best$cost - nItems*gaa$settings$itemCost)/sum(gaa$settings$sWeights)
  gaa$results$meanR2 <- c(gaa$results$meanR2, meanR2)

  if(plot) plot.GAabbreviate(gaa)

  # record results at current iteration 
  assign("TMPGAabbreviateObject", gaa, envir = env)
  # 
  invisible(gaa)
}

# Creates scoring key, generates measure statistics, etc.
measure.GAabbreviate <- function(obj) 
{  
  items <- which(obj$GA@solution[1,]==1)
  nItems <- length(items)
  key <- makeKey.GAabbreviate(obj)
  nScaleItems <- apply(key, 2, function(x) sum(abs(x)))
  obj$measure <- list(items = items,
                      nItems = nItems,
                      key = key,
                      nScaleItems = nScaleItems)  
  obj$measure$alpha <- alphas.GAabbreviate(obj)
  obj$measure$ccTraining <- convCorr.GAabbreviate(obj, "train")
  obj$measure$ccValidation <- if(obj$settings$crossVal & 
                                 length(obj$data$crossvalSubjects) > 0)
                                 convCorr.GAabbreviate(obj, "validation") 
                              else NULL
  
  return(obj$measure)
}

# Generate a scoring key for the best measure in the gaa object.
# Doesn't need to be called directly, as it will be called internally
# and set as the measure$key attribute in every gaa object.
# Returns an item x scale matrix where 1 or -1 indicates 
# direction of scoring. Compatible with scoreItems() in psych.
# Always uses the training dataset, never the validation data.
makeKey.GAabbreviate <- function(obj) 
{
  active <- which(obj$GA@solution[1,]==1)
  subs <- obj$data$trainSubjects
  items <- obj$data$items[subs,active,drop=FALSE]
  scales <- obj$data$scales[subs,,drop=FALSE]
  use <- ifelse(obj$settings$pairwise, "pairwise.complete.obs", "everything")
  cors <- cor(items, scales, use = use)
  key <- matrix(0, nrow = ncol(items), ncol = ncol(scales))
  for(i in 1:ncol(scales)) 
     { icors <- abs(cors[,i])
       xvars <- order(icors, decreasing=TRUE)
       npass <- length(which(icors[xvars] >= obj$settings$minR))
       xvars <- xvars[1:min(npass, obj$settings$maxItems)]
       key[xvars,i] = 1*sign(cors[xvars,i])
  }
  return(key)
}

# Compute alpha coefficients for all scales using scoreItems
# in the psych package. Runs on validation sample if available,
# otherwise on training sample.
alphas.GAabbreviate <- function(obj) 
{
  subs <- if(obj$settings$crossVal & length(obj$data$crossvalSubjects) > 0) 
            obj$data$crossvalSubjects
          else
          { if(obj$settings$verbose)
              warning("No validation subjects found. Computing coefficient alphas on training subjects instead. This means the estimates may be biased!")
            obj$data$trainSubjects
          }
  oldwarn <- getOption("warn") # this is needed to avoid a warning
  options(warn = -1)           # not related to computing alphas
  scores <- scoreItems(keys = obj$measure$key, 
                       items = obj$data$items[subs,obj$measure$items,drop=FALSE])
  options(warn = oldwarn)

  return(scores$alpha)
}

# Returns convergent correlations between GA-based predicted scores and
# true scores. Can be applied to either the training data or (if present)
# the validation data.
convCorr.GAabbreviate <- function(obj, mode = c("train", "validation")) 
{
  # Basic validation
  mode <- match.arg(mode)
  if(is.null(obj$measure$key)) 
    stop("The 'GAabbreviate' object has no scoring key!")
  else 
    if(mode == "validation" & length(obj$data$crossvalSubjects)==0) 
      stop("There is no validation data for this measure; all subjects were used for training.")
  
  subs <- if(mode == "train") obj$data$trainSubjects 
          else obj$data$crossvalSubjects
  active <- which(obj$GA@solution[1,]==1)
  if(length(active) == 0) 
    stop("There are no items left in the abbreviated measure! Aborting.")
  scores <- obj$data$items[subs,active,drop=FALSE] %*% obj$measure$key
  cc <- diag(cor(scores, obj$data$scales[subs,], use = "pairwise"))
  if(!is.null(colnames(obj$data$scales))) 
    names(cc) <- colnames(obj$data$scales)
  return(cc)
}

# Return basic information about the best solution.
print.GAabbreviate <- function(x, ...)
{ print(str(x,2), ...) }

summary.GAabbreviate <- function(object, digits = getOption("digits"), verbose = FALSE, ...)
{
  out <- capture.output(print(summary(object$GA), digits = digits))
  out <- out[seq(grep("Fitness function value", out, value = FALSE)-1)]
  out <- c(out, 
           paste(format("Total cost", justify = "left", width = 22), "=", 
                 signif(object$best$cost, digits = digits)),
           paste(format("Number of items in initial set", justify = "left", width = 22), "=", 
                 signif(object$data$nItems, digits = digits)),
           paste(format("Number of items in final set", justify = "left", width = 22), "=", 
                 signif(object$measure$nItems, digits = digits)),
           paste(format("Mean coefficient alpha", justify = "left", width = 22), "=", 
                 signif(mean(object$measure$alpha), digits = digits-3)),
           paste(format("Mean convergent correlation (training)", justify = "left", width = 22), "=", 
                 signif(mean(object$measure$ccTraining), digits = digits-3))
           )
  
  if(object$settings$crossVal)
    { out <- c(out, 
               paste(format("Mean convergent correlation (validation)", justify = "left", width = 22), "=", 
                 signif(mean(object$measure$ccValidation), digits = digits-3))) }

  cat(out, sep = "\n")
  
  if(verbose)
    { cat("\nSelected items:\n")
      i <- which(object$results$solution[length(object$results$nItems),] == 1)  
      names(i) <- colnames(object$data$items)[i]
      print(i)  
  }
}

# Plot method 
plot.GAabbreviate <- function(x, ...) 
{
  object <- x  # Argh.  Really want to use 'object' anyway
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  layout(cbind(c(1,2,3), c(5,5,5), c(4,4,4)))
  # layout.show(5)
  par(mgp = c(2,0.7,0), plt = c(0.1,0.85,0.2,0.85), 
      mar = c(3,3,2,.5)+0.1)
  
  # Plot line graph of total cost
  iter <- length(object$results$cost)
  plot(seq(iter), object$results$cost, type="o",
       col = "red3", cex = 0.5, xaxt = "n",
       ylim = extendrange(object$results$cost, f = 0.05),
       xlab = "GA generation", ylab = "Minimum cost",
       main = "Total cost", panel.first = grid())
  axis(side = 1, at = { at <- pretty(seq(iter), n = 5, high.u.bias = 2)
                        at[at==0] <- 1; at } )
  
  # Plot line graph of no. of items
  plot(seq(iter), object$result$nItems, type="o",
       col = "darkorange2", cex = 0.5,
       ylim = extendrange(object$results$nItems, f = 0.05),
       xlab = "GA generation", ylab = "No. of items",
       main = "Measure length", yaxt = "n", xaxt = "n",
       panel.first = grid())
  axis(side = 2, at = sort(unique(object$result$nItems)))
  axis(side = 1, at = { at <- pretty(seq(iter), n = 5, high.u.bias = 2)
                        at[at==0] <- 1; at } )
  
  # Plot line graph of mean R^2 (across scales) 
  plot(seq(iter), object$result$meanR2, type = "o",
       col = "dodgerblue3", cex = 0.5, ylim = c(0,1), xaxt = "n",
       xlab = "GA generation", ylab = "Mean R^2 per scale",
       main = "Mean R^2", panel.first = grid())
  axis(side = 1, at = { at <- pretty(seq(iter), n = 5, high.u.bias = 2)
                        at[at==0] <- 1; at } )
  
  # Plot history of best items
  l <- length(object$results$solution[1,])
  par(mar = c(3,1,2,3)+0.1, mgp = c(2,0.5,0))
  image(y = (0:iter), x = 1:l, z = t(object$results$solution),
        ylim = c(iter,0), xlim = c(1, l),
        col = c("grey95", "grey25"), xlab = "", ylab = "",
        xaxs = "r", xaxt = "n", yaxt = "n", bty = "n")
  title(main = "Items included in measure")
  box()
  axis(side = 1, at = { at <- pretty(seq(l), n = 5, high.u.bias = 2)
                        at[at==0] <- 1; at } )
  atIter <- { at <- unique(round(pretty(seq(iter), n = 5, high.u.bias = 2)))
              at[at==0] <- 1; at[at<=iter] }
  axis(side = 4, at = atIter, labels = atIter)
  usr <- par("usr")
  mtext(side = 4, text = "GA generation", at = mean(usr[3:4]),
        line = 2, cex = 0.7)
  # axis(4, at = atIter, labels = FALSE)
  # text(usr[2]+max(strwidth(atIter))*1.5, atIter,
  #      labels = atIter, srt = 270, xpd = TRUE)
  # text(usr[2]+max(strwidth(atIter))*2.5, mean(usr[3:4]),
  #      labels = "GA generation", srt = 270, xpd = TRUE)

  # Plot details of current best solution
  nScales <- object$data$nScales
  par(mgp = c(2,0.7,0), mar = c(3,4,2,1)+0.1)
  plot(0, 0, type = "n", xlim = c(0, 1),
       ylim = c(nScales, 1), yaxs = "r", yaxt = "n", # xaxs="i",
       xlab = "Variance explained (R^2)", ylab = "",
       main = "Best solution", xaxt = "n")
  axis(side = 1, at = pretty(0:1, n = 5, high.u.bias = 2))
  points(x = (1-object$best$fit), y = 1:nScales, type = "p", pch = 15)
  abline(v = (1:4)/5, lty = 2, col = "darkgray", lwd = 0.5)
  labels <- colnames(object$data$scales)
  if(is.null(labels)) labels <- 1:nScales
  mtext(labels, side = 2, line = 0.5, at = 1:nScales, cex = 0.7, las = 2)

  invisible(object)
}

is.missing <- function(x) 
{
  x <- as.matrix(x)
  (is.na(x) | !is.finite(x))
}

# Helper function; imputes all NAs with the column-wise mean.
impute <- function(x) 
{
  cnames <- names(x)
  x <- apply(x, 2, function(y) { y[is.na(y)] <- mean(y, na.rm=TRUE); y })
  colnames(x) <- cnames
  x
}
