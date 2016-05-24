
#' A3 Results for Arbitrary Model
#'
#' This function calculates the A3 results for an arbitrary model construction algorithm (e.g. Linear Regressions, Support Vector Machines or Random Forests). For linear regression models, you may use the \code{\link{a3.lm}} convenience function.
#'
#' @param formula the regression formula.
#' @param data a data frame containing the data to be used in the model fit.
#' @param model.fn the function to be used to build the model.
#' @param model.args a list of arguments passed to \code{model.fn}.
#' @param ... additional arguments passed to \code{\link{a3.base}}.
#' @return S3 \code{A3} object; see \code{\link{a3.base}} for details
#' @references Scott Fortmann-Roe (2015). Consistent and Clear Reporting of Results from Diverse Modeling Techniques: The A3 Method. Journal of Statistical Software, 66(7), 1-23. <http://www.jstatsoft.org/v66/i07/>
#' @examples
#' \donttest{
#'  ## Standard linear regression results:
#'  
#'  summary(lm(rating ~ ., attitude))
#'  
#'  ## A3 Results for a Linear Regression model:
#'  
#'  # In practice, p.acc should be <= 0.01 in order
#'  # to obtain finer grained p values.
#'  
#'  a3(rating ~ ., attitude, lm, p.acc = 0.1)
#'  
#'  
#'  ## A3 Results for a Random Forest model:
#'  
#'  # It is important to include the "+0" in the formula
#'  # to eliminate the constant term.
#'  
#'  require(randomForest)
#'  a3(rating ~ .+0, attitude, randomForest, p.acc = 0.1)
#'  
#'  # Set the ntrees argument of the randomForest function to 100
#'  
#'  a3(rating ~ .+0, attitude, randomForest, p.acc = 0.1, model.args = list(ntree = 100))
#'  
#'  # Speed up the calculation by doing 5-fold cross-validation.
#'  # This is faster and more conservative (i.e. it should over-estimate error)
#'  
#'  a3(rating ~ .+0, attitude, randomForest, n.folds = 5, p.acc = 0.1)
#'  
#'  # Use Leave One Out Cross Validation. The least biased approach,
#'  # but, for large data sets, potentially very slow.
#'  
#'  a3(rating ~ .+0, attitude, randomForest, n.folds = 0, p.acc = 0.1)
#'  
#'  ## Use a Support Vector Machine algorithm.
#'  
#'  # Just calculate the slopes and R^2 values, do not calculate p values.
#'  
#'  require(e1071)
#'  a3(rating ~ .+0, attitude, svm, p.acc = NULL)
#'  }

a3 <- function(formula, data, model.fn, model.args = list(), ...){
  
  model.fn.w.args <- function(y, x){
    dat <- data.frame(cbind(y, x))
    
    names(dat) <-  c("y", paste("x", 1:ncol(x), sep=""))
    
    new.model.args = list(formula = y ~ . + 0, data = dat)
    for(n in names(model.args)){
      new.model.args[[n]] = model.args[[n]]
    }
    
    return(do.call(model.fn, new.model.args))
  }
  
  simulate.fn <- function(y, x, new.x, ...){
    reg <- model.fn.w.args(y, x)
    new.data <- data.frame(new.x)
    if(ncol(new.data) != ncol(x)){
      new.data <- data.frame(t(new.data))
    }
    names(new.data) <- paste("x", 1:ncol(x), sep="")
    return(predict(reg, new.data))
  }
  
  a3.base(formula, data, model.fn.w.args, simulate.fn, ...)
}

#' A3 for Linear Regressions
#'
#' This convenience function calculates the A3 results specifically for linear regressions. It uses R's \code{\link{glm}} function and so supports logistic regressions and other link functions using the \code{family} argument. For other forms of models you may use the more general \code{\link{a3}} function.
#'
#' @param formula the regression formula.
#' @param data a data frame containing the data to be used in the model fit.
#' @param family the regression family. Typically 'gaussian' for linear regressions.
#' @param ... additional arguments passed to \code{\link{a3.base}}.
#' @return S3 \code{A3} object; see \code{\link{a3.base}} for details
#' @examples
#' \donttest{
#'  ## Standard linear regression results:
#'  
#'  summary(lm(rating ~ ., attitude))
#'  
#'  ## A3 linear regression results:
#'  
#'  # In practice, p.acc should be <= 0.01 in order
#'  # to obtain fine grained p values.
#'  
#'  a3.lm(rating ~ ., attitude, p.acc = 0.1)
#'  
#'  # This is equivalent both to:
#'  
#'  a3(rating ~ ., attitude, glm, model.args = list(family = gaussian), p.acc = 0.1)
#'  
#'  # and also to:
#'  
#'  a3(rating ~ ., attitude, lm, p.acc = 0.1)
#'  }

a3.lm <- function(formula, data, family = gaussian, ...){
  a3(formula, data, model.fn = glm, model.args = list(family = family), ...)
}


#' Base A3 Results Calculation
#'
#' This function calculates the A3 results. Generally this function is not called directly. It is simpler to use \code{\link{a3}} (for arbitrary models) or \code{\link{a3.lm}} (specifically for linear regressions).
#'
#' @param formula the regression formula.
#' @param data a data frame containing the data to be used in the model fit.
#' @param model.fn function used to generate a model.
#' @param simulate.fn function used to create the model and generate predictions.
#' @param n.folds the number of folds used for cross-validation. Set to 0 to use Leave One Out Cross Validation.
#' @param data.generating.fn the function used to generate stochastic noise for calculation of exact p values.
#' @param p.acc the desired accuracy for the calculation of exact p values. The entire calculation process will be repeated \eqn{1/p.acc} times so this can have a dramatic affect on time required. Set to \code{NULL} to disable the calculation of p values.
#' @param features whether to calculate the average slopes, added \eqn{R^2} and p values for each of the features in addition to the overall model.
#' @param slope.sample if not NULL the sample size for use to calculate the average slopes (useful for very large data sets).
#' @param slope.displacement the amount of displacement to take in calculating the slopes. May be a single number in which case the same slope is applied to all features. May also be a named vector where there is a name for each feature.
#' @return S3 \code{A3} object containing:
#' \item{model.R2}{The cross validated \eqn{R^2} for the entire model.}
#' \item{feature.R2}{The cross validated \eqn{R^2}'s for the features (if calculated).}
#' \item{model.p}{The p value for the entire model (if calculated).}
#' \item{feature.p}{The p value for the features (if calculated).}
#' \item{all.R2}{The \eqn{R^2}'s for the model features, and any stochastic simulations for calculating exact p values.}
#' \item{observed}{The observed response for each observation.}
#' \item{predicted}{The predicted response for each observation.}
#' \item{slopes}{Average slopes for each of the features (if calculated).}
#' \item{all.slopes}{Slopes for each of the observations for each of the features (if calculated).}
#' \item{table}{The A3 results table.}
#' 
a3.base <- function(formula, data, model.fn, simulate.fn,  n.folds = 10, data.generating.fn = replicate(ncol(x), a3.gen.default), p.acc = 0.01, features = TRUE, slope.sample = NULL, slope.displacement = 1){
  if(! is.null(p.acc)){
    if(p.acc <= 0 || p.acc >=1){
      stop("p.acc must be between 0 and 1. Set p.acc to NULL to disable the calculation of p values.")
    }
  }
  if(n.folds < 2 && n.folds != 0){
    stop("n.folds must be >= 2. Set n.folds to 0 to use Leave One Out Cross Validation.")
  }
  
  n.reps <- 0
  if(! is.null(p.acc)){
    n.reps <- ceiling(1/p.acc)
  }
  
  res <- list()
  mf <- model.frame(formula, data, drop.unused.levels = TRUE)
  x <- model.matrix(formula, mf)
  y <- model.response(mf)
  
  if(length(data.generating.fn) != ncol(x)){
    stop("data.generating.fn must be a list of functions one for each column in the model matrix")
  }
  
  if(n.folds == 0){
    n.folds <- length(y)
  }
  
  my.apply <- lapply
  if( ! is.null(p.acc) ){
    # if( library(pbapply, quietly = TRUE, logical.return = TRUE) == TRUE ){ # not needed due to depends
      my.apply <- pblapply # Show a progress bar if available
    # }
  }
  
  # Calculate the groups for cross validation
  cv.folds <- split(sample(1:length(y)), rep(1:n.folds, length = length(y)))
  
  # Generate random data series for p values
  new.data <- lapply(1:ncol(x), function(c){
    data.generating.fn[[c]](x[,c], n.reps)
  })
  
  
  r2.formatter <- function(x){
    signs <- sign(x)
    x <- abs(x)*100
    res <- paste(format(round(x, 1), digits = 3), "%")
    signs <- sapply(signs, function(x){
      if(x == -1){
        return("- ")
      }else{
        return("+ ")
      }
    })
    if(signs[1] == "+ "){
      signs[1] <- "  " # removed the plus sign for the overall model accuracy
    }
    res <- paste(signs, res, sep="")
    return(res)
  }
  
  p.formatter <- function(x){
    if(length(x) == 0){
      return(c()) 
    }
    
    res <- format(x, digits = 4)
    for(i in 1:length(x)){
      if(x[i] == 0){
        res[i] <- paste("<", p.acc)
      }
    }
    return(res)
  }
  
  slope.formatter <- function(x){
    if(length(x) == 0){
      return(c()) 
    }
    
    format(x, digits = 3)
  }
  
  # Setup iterations
  # "default" is initial simulation without any randomized data
  # Each rep after that has some form of randomized data
  iterations <- "default"
  if(! is.null(p.acc)){
    iterations <- c(iterations, 1:n.reps)
  }
  
  top <- 0
  if(features){
    top <- ncol(x)
  }
  
  # Iterate through each rep and the default
  # outputs[[1]] will have the set of default cases in it
  # outputs[[>1]] will have the randomized data cases
  outputs <- my.apply(iterations, function(rep){
    # Calculate R2's for the rep
    # We calculate for the model (0) and then for each column of data by numerical index
    out <- lapply(0:top, function(c){
      
      new.x <- x
      if(rep != "default"){
        # If we aren't on the default case, we add some form of randomization
        if(c==0){
          # We are doing the overall model
          # So randomize all the data
          
          for(j in 1:ncol(x)){
            new.x[,j] <- new.data[[j]][[as.numeric(rep)]]
          }
        }else{
          # We're looking at a specific column, so just randomize that data
          new.x[,c] <- new.data[[c]][[as.numeric(rep)]]	
        }
      }
      
      # Remove a column of data if we are at the un-randomized case
      if((c != 0) && (rep == "default")){
        if(top == 1){
          return(list(R2=0)); # if there is only one feature column added R^2 should be full value
        }
        new.x <- as.data.frame(new.x[,-c])
      }
      
      res <- a3.r2(y, new.x, simulate.fn, cv.folds)
      return(res)
      
    }
    )
    r2 <- sapply(out, function(x){x$R2})
    
    return(list( R2 = r2, predicted = out[[1]]$predicted, observed = out[[1]]$observed ))
  })
  
  predicted <- outputs[[1]]$predicted
  observed <- outputs[[1]]$observed
  
  outputs <- lapply(outputs, function(x){x$R2})
  
  if(features){
    get.names <- function(formula, data){
      #t <- terms(formula, data=data)
      #l <- attr(t, "term.labels")
      #if(attr(t, "intercept")==1){
      #  l <- c("(Intercept)",l)
      #}
      return(attr(x, "dimnames")[[2]])
    }
    entry.names <- c("-Full Model-", get.names(formula, data = data))
    
    getSlopes <- function (reg, data){
      slopes <- list()
      for(col in 2:ncol(data)){
        slopes[[as.character(col)]] <- c()
        span <- range(data[,col])
        span <- span[2] - span[1]
        
        for(row in 1:nrow(data)){
          
          point <- data[row,]
          at.point <- predict(reg, point)
          
          if(length(slope.displacement) == 1){
            dist <- slope.displacement
          }else{
            dist <- slope.displacement[entry.names[col]]
          }
          
          slope <- 0
          
          #while(TRUE){
            
            above.point <- point
            above.point[col] <- point[col] + dist
            at.above <- predict(reg, above.point)[[1]]
            below.point <- point
            below.point[col] <- point[col] - dist
            at.below <- predict(reg, below.point)[[1]]
            new.slope <- (at.above - at.below)/(dist*2)
#             
#             if(new.slope == 0){
#               dist <- dist * 2
#               if(slope != 0){
#                 break
#               }
#               if(dist > span){
#                 break
#               }
#             }else{
#               if(abs( (slope-new.slope) / new.slope) < epsilon){
#                 break
#               }
#               dist <- dist / 2
#             }
            
            slope <- new.slope
          #}
          slopes[[as.character(col)]] <- c(slopes[[as.character(col)]], slope)
        }
      }
      slopes
    }
    # print(model.fn(y, x))
    slope.data <- data.frame(cbind(y, x))
    names(slope.data) <-  c("y", paste("x", 1:ncol(x), sep=""))
	if(! is.null(slope.sample)){
		slope.data <- slope.data[sample(1:nrow(slope.data), slope.sample),]
	}
    res[["all.slopes"]] <- getSlopes(model.fn(y, x), slope.data)
    res[["all.slopes"]] <- lapply(res[["all.slopes"]], function(x){ round(x, digits = 8)})
    res[["slopes"]] = sapply(res[["all.slopes"]], function(x){
      return(median(x))
#       r <- unique(round(range(x), digits=8))
#       if(length(r) == 1){
#         return(r)
#       }else{
#         return(paste(r, collapse = " - "))
#       }
    } )
    names(res[["all.slopes"]]) <- entry.names[-1]
    names(res[["slopes"]]) <- entry.names[-1]
    
  }else{
    entry.names <- c("-Full Model-")
  }
  
  # Now take the data and calculate R2 and p values
  res[["predicted"]] <- predicted
  res[["observed"]] <-  observed
  if(! is.null(p.acc)){
    # we did a set of repitions so we should calculate  p value
    r2 <- c()
    p.values <- c()
    res$all.R2 <- list()
    
    # item 1 is the overall model
    # item > 1 is a specific column
    for(i in 1:(top+1)){
      # Get the R2 for the specific item
      items <- sapply(outputs, function(x){x[i]}) 
      
      # check if we are on a column item
      if(i > 1){
        # if so replace the first element of the items list with the value of the overall model R2
        items[1] <- r2[1]
      }
      
      names(items) <- c("Base",paste("Rep", 1:n.reps))
      res$all.R2[[paste(entry.names[i])]] <- items
      
      # rank the items by R2
      dist <- rank(items)
      
      # find the position of the first item (overall model R2) within the list of randomly derived p values
      # this is the emprical R2
      p.values <- c(p.values, 1 - (dist[1]-1)/(length(dist)-1))
      
      # the R2 of the model as generated with random data
      new.null <- mean(items[-1])
      
      if(i==1){
        # the R2 of the model as generated with random data
        new.null <- mean(items[-1])
        
        # if the R2 with stochasticity is better than 0, we will use as a baseline to scale our R2
        #if(new.null > 0){ # Adjust R^2 based on what was observed in stochastic series XXX reenable?
        #  r2 <- (items[1]-new.null) / (1-new.null)
        #}else{
          r2 <- items[1]
        #}
      }else{
        # get data series again (we overwrote the item[1] position earlier)
        items <- sapply(outputs, function(x){x[i]})
        
        # see how the results improve compared to the baseline
        r2 <- c(r2, r2[1] - items[1])#max(new.null, items[1])) # Adjust R^2 based on what was observed in stochastic series XXX reenable?
      } 
    }
    
    res$table <- data.frame(`Average Slope` = c("", slope.formatter(res$slopes)), `CV R^2` = r2.formatter(r2), `p value` = p.formatter(p.values), check.names=F)
    res$model.R2 <- r2[1]
    res$feature.R2 <- r2[-1]
    res$model.p <- p.values[1]
    names(res$model.p) <- entry.names[1]
    res$feature.p <- p.values[-1] 
    names(res$feature.p) <- entry.names[-1]
  }else{
    # we didn't do repetitions so no p values
    
    r2 <- outputs[[1]]
    r2[-1] <- r2[1] - r2[-1] # get delta to overall model R2
    
    res$table <- data.frame(`Average Slope` = c("", slope.formatter(res$slopes)), `CV R^2` = r2.formatter(r2), check.names=F)
    res$all.R2 <- r2
    res$model.R2 <- r2[1]
    res$feature.R2 <- r2[-1]
  }
  
  names(res$model.R2) <- entry.names[1]
  names(res$feature.R2) <- entry.names[-1]
  
  row.names(res$table) <- entry.names
  
  class(res) <- "A3"
  
  return(res)
}

#' Plot A3 Results
#'
#' Plots an 'A3' object results. Displays predicted versus observed values for each observation along with the distribution of slopes measured for each feature.
#'
#' @param x an A3 object.
#' @param ... additional options provided to \code{\link{plotPredictions}}, \code{\link{plotSlopes}} and \code{\link{plot}} functions.
#' @method plot A3
#' @examples
#' \donttest{
#'  data(housing)
#'  res <- a3.lm(MED.VALUE ~ NOX + ROOMS + AGE + HIGHWAY + PUPIL.TEACHER, housing, p.acc = NULL)
#'  plot(res)
#'  }

plot.A3 <- function(x, ...){
  if(class(x) != "A3"){
    stop("'x' must be of class 'A3'.")
  }
  
  plotPredictions(x, ...)
  old.par <- par(ask=T)
  plotSlopes(x, ...)
  par(old.par)
}

#' Plot Predicted versus Observed
#'
#' Plots an 'A3' object's values showing the predicted versus observed values for each observation.
#'
#' @param x an A3 object,
#' @param show.equality if true plot a line at 45-degrees.
#' @param xlab the x-axis label.
#' @param ylab the y-axis label.
#' @param main the plot title.
#' @param ... additional options provided to the \code{\link{plot}} function.
#' 
#' @examples
#'  data(multifunctionality)
#'  x <- a3.lm(MUL ~ ., multifunctionality, p.acc = NULL, features = FALSE)
#'  plotPredictions(x)

plotPredictions <- function(x, show.equality = TRUE, xlab = "Observed Value", ylab = "Predicted Value", main = "Predicted vs Observed", ...){
  if(class(x) != "A3"){
    stop("'x' must be of class 'A3'.")
  }
  
  plot(x$observed, x$predict, xlab=xlab, ylab=ylab, main=main, ...)
  abline(h=0, col="Gray"); abline(v=0, col="Gray");
  if(show.equality){
    abline(coef = c(0, 1), col = "Blue", lty = 2)
  }
}

#' Plot Distribution of Slopes
#'
#' Plots an 'A3' object's distribution of slopes for each feature and observation. Uses Kernel Density Estimation to create an estimate of the distribution of slopes for a feature.
#'
#' @param x an A3 object.
#' @param ... additional options provided to the \code{\link{plot}} and \code{\link{density}} functions.
#' 
#' @examples
#' \donttest{
#'  require(randomForest)
#'  data(housing)
#'  
#'  x <- a3(MED.VALUE ~ NOX + PUPIL.TEACHER + ROOMS + AGE + HIGHWAY + 0,
#'    housing, randomForest, p.acc = NULL, n.folds = 2)
#'  
#'  plotSlopes(x)
#'  }

plotSlopes <- function(x,  ...){
  if(class(x) != "A3"){
    stop("'x' must be of class 'A3'.")
  }
  size <- length(x$slopes)
  if(size == 0 ){
    stop("no slopes to plot")
  }
  width <- ceiling(sqrt(size))
  height <- floor(sqrt(size))
  if(width*height < size){
    width <- width+1
  }
  
  old.par <- par(mfrow = c(height, width), mar = .55*c(5, 4, 4, 2) + 0.2)
  
  for(s in names(x$slopes)){
    if(length(unique(x$all.slopes[[s]]))==1){
      plot(x$slopes[[s]], 0, main = s)
    }else{
      plot(density(x$all.slopes[[s]],...), xlab = "", ylab="", main = s, ...)
      rug(x$all.slopes[[s]], col="Blue")
    }
    abline(h=0, col="Gray")
    abline(v=0, col="Gray")
  }
  
  par(old.par)
}

#' Print Fit Results
#'
#' Prints an 'A3' object results table.
#'
#' @param x an A3 object.
#' @param ... additional arguments passed to the \code{\link{print}} function.
#' @method print A3
#' @examples
#'  x <- a3.lm(rating ~ ., attitude, p.acc = NULL)
#'  print(x)

print.A3 <- function(x, ...){
  if(class(x) != "A3"){
    stop("'x' must be of class 'A3'.")
  }
  
  print(x$table, ...)
}

#' Nicely Formatted Fit Results
#'
#' Creates a LaTeX table of results. Depends on the \pkg{xtable} package.
#'
#' @param x an A3 object.
#' @param ... additional arguments passed to the \code{\link{print.xtable}} function.
#' @method xtable A3
#' @examples
#'  x <- a3.lm(rating ~ ., attitude, p.acc = NULL)
#'  xtable(x)

xtable.A3 <- function(x, ...){
  # require(xtable) # not needed due to depends
  if(class(x) != "A3"){
    stop("'x' must be of class 'A3'.")
  }
  
  data <- x$table
  names(data) <- gsub("p value","Pr(>R^2)", names(data), fixed=TRUE) 
  names(data) <- gsub("R^2","$R^2$", names(data), fixed=TRUE)
  print(xtable(data, align=c("l","|", rep("r", ncol(data))), ...), sanitize.colnames.function = function(x){x}, sanitize.text.function = function(x){
    return(sapply(x, function(x){
      trimmed = gsub("(^ +)|( +$)", "", x)
      if((! is.na(suppressWarnings(as.numeric(x)))) && suppressWarnings(as.numeric(x))==trimmed){
        return(paste0("$", x, "$"))
      }else if(substr(trimmed, nchar(trimmed), nchar(trimmed))=="%"){
        return(paste0("$", gsub("%", "\\%", trimmed, fixed=T), "$"));
      }else if(substr(trimmed, 1, 1)=="<"){
        return(paste0("$", trimmed, "$"));
      }else{
        return (gsub("_", "\\_", x, fixed=T));
      }
    }))
  })
}


#' Cross-Validated \eqn{R^2}
#'
#' Applies cross validation to obtain the cross-validated \eqn{R^2} for a model: the fraction of the squared error explained by the model compared to the null model (which is defined as the average response). A pseudo \eqn{R^2} is implemented for classification.
#'
#' @param y a vector or responses.
#' @param x a matrix of features.
#' @param simulate.fn a function object that creates a model and predicts y.
#' @param cv.folds the cross-validation folds.
#' 
#' @return A list comprising of the following elements:
#' \item{R2}{the cross-validated \eqn{R^2}}
#' \item{predicted}{the predicted responses}
#' \item{observed}{the observed responses}

a3.r2 <- function(y, x, simulate.fn, cv.folds){
  
  errors <- lapply(cv.folds, function(fold){
    test.y <- y[fold]
    test.x <- x[fold,]
    train.y <- y[-fold]
    train.x <- as.data.frame(x[-fold,])
    new.y <- simulate.fn(train.y, train.x, test.x)
    
    if(is.factor(new.y)){
      # classification
      return(list( type="classification", correct = (new.y == test.y), predicted = new.y, observed = test.y ))
    }else{
      # regression
      y.null <- mean(train.y)
      return(list( type="regression", ss.null = sum((test.y-y.null)^2), ss.model = sum((test.y-new.y)^2), predicted = new.y, observed = test.y ))
    }
  }
  )
  
  if(errors[[1]]$type == "regression"){
    # regression
    ss.model <- sum(unlist(sapply(errors, function(x){x$ss.model})))
    ss.null <- sum(unlist(sapply(errors, function(x){x$ss.null})))
    return( list(R2 = 1 - ss.model/ss.null, predicted = unlist(sapply(errors, function(x){x$predicted})), observed = unlist(sapply(errors, function(x){x$observed}))) )
  }else{
    # classification
    corrects <- unlist(sapply(errors, function(x){x$correct}))
    
    null.count <- max(table(as.factor(y))) # return number of observations in the largest class
    
    return( list( R2 = (sum(corrects)-null.count) / (length(corrects)-null.count), predicted = unlist(sapply(errors, function(x){x$predicted})), observed = unlist(sapply(errors, function(x){x$observed}))) )
  }
}

#' Stochastic Data Generators
#'
#' The stochastic data generators generate stochastic noise with (if specified correctly) the same properties as the observed data. By replicating the stochastic properties of the original data, we are able to obtain the exact calculation of p values.
#' 
#' Generally these will not be called directly but will instead be passed to the \code{data.generating.fn} argument of \code{\link{a3.base}}.
#'
#' @name a3.gen.default
#' @aliases a3.gen.default a3.gen.bootstrap a3.gen.resample a3.gen.normal a3.gen.autocor
#' 
#' @param x the original (observed) data series.
#' @param n.reps the number of stochastic repetitions to generate.
#' 
#' @return A list of of length \code{n.reps} of vectors of stochastic noise. There are a number of different methods of generating noise:
#' \item{a3.gen.default}{The default data generator. Uses \code{a3.gen.bootstrap}.}
#' \item{a3.gen.resample}{Reorders the original data series.}
#' \item{a3.gen.bootstrap}{Resamples the original data series with replacement.}
#' \item{a3.gen.normal}{Calculates the mean and standard deviation of the original series and generates a new series with that distribution.}
#' \item{a3.gen.autocor}{Assumesa first order autocorrelation of the original series and generates a new series with the same properties.}
#' 
#' @examples
#' \donttest{
#'  # Calculate the A3 results assuming an auto-correlated set of observations.
#'  # In usage p.acc should be <=0.01 in order to obtain more accurate p values.
#'  
#'  a3.lm(rating ~ ., attitude, p.acc = 0.1,
#'    data.generating.fn = replicate(ncol(attitude), a3.gen.autocor))
#'  }
#'  
#'  ## A general illustration:
#'  
#'  # Take x as a sample set of observations for a feature
#'  x <- c(0.349, 1.845, 2.287, 1.921, 0.803, 0.855, 2.368, 3.023, 2.102, 4.648)
#'  
#'  # Generate three stochastic data series with the same autocorrelation properties as x
#'  rand.x <- a3.gen.autocor(x, 3)
#'  
#'  plot(x, type="l")
#'  for(i in 1:3) lines(rand.x[[i]], lwd = 0.2)

# Default generator, use bootstrap 
a3.gen.default <- function(x, n.reps){
  if(length(unique(x))==1){
    #it's a constant, such as an intercept
    return(a3.gen.normal(x, n.reps))
  }
  a3.gen.bootstrap(x, n.reps)
}

# Generates a bootstrap random data series
a3.gen.bootstrap <- function(x, n.reps){  
  res <- lapply(1:n.reps, function(r) {sample(x, length(x), replace=TRUE)})
  res$default <- x
  res
}

# Generates a resampled random data series
a3.gen.resample <- function(x, n.reps){	
  res <- lapply(1:n.reps, function(r) {sample(x, length(x), replace=FALSE)})
  res$default <- x
  res
}

# Generates a normally distributed random data series
a3.gen.normal <- function(x, n.reps){
  mu <- mean(x)
  sd <- sd(x)
  if(sd == 0){
    sd <- 1
  }
  res <- lapply(1:n.reps, function(r) {rnorm(length(x), mu, sd)})
  res$default <- x
  res
}

# Generates a first order autocorrelated random data series
a3.gen.autocor <- function(x, n.reps){
  mu <- mean(x)
  sd <- sd(x)
  if(sd == 0){
    r <- 1
  }else{
    r <- cor(x[-1], x[-length(x)])
  }
  res <- lapply(1:n.reps,
                function(rep) {
                  dat <- rnorm(length(x), mu, sd)
                  for(i in 2:length(x)){
                    dat[i] <- dat[i-1]*r + dat[i]*(1-r)
                  }
                  return(dat)
                }
  )
  res$default <- x
  res
}

