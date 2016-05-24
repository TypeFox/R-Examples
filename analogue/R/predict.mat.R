###########################################################################
##                                                                       ##
## predict.mat() - 'predict' method for MAT models                       ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 27-May-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## object        - object on which method dispatch applied (Only 'mat')  ##
## newdata       - data frame, required only if bootstrap predictions &  ##
##                 sample specific errors are required. Must have same   ##
##                 number of columns, in same order, as 'x' in mat().    ##
##                 See example in ?join for how to do this. If not       ##
##                 provided, the fitted values are returned.             ##
## k             - number of analogues to use. If missing 'k' is chosen  ##
##                 automatically as the 'k' that achieves lowest RMSE.   ##
## weighted      - Logical. Should the analysis use weighted mean of env ##
##                 data of analogues as fitted/estimated values?         ##
## bootstrap     - Should bootstrap-derived estimates and samples        ##
##                 specific erros be calculated - ignored if 'newdata'   ##
##                 is missing.                                           ##
## n.boot        - number of bootstrap samples to take.                  ##
##                                                                       ##
###########################################################################
predict.mat <- function(object, newdata, k, weighted = FALSE,
                        bootstrap = FALSE, n.boot = 1000,
                        probs = c(0.01, 0.025, 0.05, 0.1), ...)
  {
    ## if no newdata then return the relevant fitted values
    if(missing(newdata))
      return(fitted(object, weighted = weighted))
    ## if not doing bootstrapping, just return the predictions
    if(!bootstrap) {
      ## automatically choose k?
      auto <- FALSE
      if(missing(k))
        {
          auto <- TRUE
          if(weighted)
            k <- which.min(object$weighted$rmsep)
          else
            k <- which.min(object$standard$rmsep)
        }
      Dij <- distance(x = object$orig.x, y = newdata,
                      method = object$method)
      minDC <- apply(Dij, 2, function(x) {sort(x)[1]})
      quantiles <- quantile(as.dist(object$Dij), probs = probs)
      if(weighted)
        predicted <- apply(Dij, 2, cumWmean, object$orig.y, drop = FALSE)
      else
        predicted <- apply(Dij, 2, cummean, object$orig.y, drop = FALSE)
      est <- fitted(object, weighted = weighted)$estimated
      obs <- object$orig.y
      resi <- resid(object, k = k, weighted = weighted)
      if(weighted) {
        model <- object$weighted$rmsep
        r2 <- object$weighted$r.squared
        avg.bias <- object$weighted$avg.bias
        max.bias <- object$weighted$max.bias
      } else {
        model <- object$standard$rmsep
        r2 <- object$standard$r.squared
        avg.bias <- object$standard$avg.bias
        max.bias <- object$standard$max.bias
      }
      res <- list(observed = obs,
                  model = list(estimated = est,
                    residuals = resi, r.squared = r2, avg.bias = avg.bias,
                    max.bias = max.bias, rmsep = model, k = k),
                  weighted = weighted, auto = auto,
                  method = object$method,
                  quantiles = quantiles,
                  predictions = list(model =
                    list(predicted = predicted,
                         k = k)),
                  minDC = minDC,
                  Dij = Dij)
    } else {
      if(missing(k))
        k <- NULL
      res <- bootstrap(object, newdata, k = k, weighted = weighted,
                       n.boot = n.boot)
    }
    class(res) <- "predict.mat"
    return(res)
  }

print.predict.mat <- function(x, digits = max(3, getOption("digits") - 3),
                              ...)
  {
    cat("\n")
    writeLines(strwrap("Modern Analogue Technique predictions", prefix = "\t"))
    cat("\n")
    cat(paste("Dissimilarity:", x$method, "\n"))
    cat(paste("k-closest analogues: ", x$predictions$model$k,
              ",\tChosen automatically? ", x$auto, "\n",
              sep = ""))
    cat(paste("Weighted mean:", x$weighted, "\n"))
    cat(paste("Bootstrap estimates:", ifelse(is.null(x$bootstrap),
                                             FALSE, TRUE),
              "\n"))
    if(!is.null(x$bootstrap)){
      cat(paste("Number of bootstrap cycles:", x$n.boot, "\n"))
      cat("\nModel and bootstrap-derived error estimates:\n\n")
      boot.errors <- with(x$bootstrap,
                          c(rmsep[k], s1[k], s2[k], r.squared[k],
                            avg.bias[k], max.bias[k]))
      model.errors <- with(x$model, c(rmsep[k], NA, NA, r.squared[k],
                                            avg.bias[k], max.bias[k]))
      errors <- rbind(model.errors, boot.errors)
      colnames(errors) <- c("RMSEP", "S1", "S2", "r.squared",
                            "avg.bias", "max.bias")
      rownames(errors) <- c("Model", "Bootstrap")
      print(errors, digits = digits, na.print = "-")
    } else {
      cat("\nModel error estimates:\n")
      model.errors <- with(x$model, c(rmsep[k], r.squared[k],
                                            avg.bias[k], max.bias[k]))
      names(model.errors) <- c("RMSEP", "r.squared",
                                  "avg.bias", "max.bias")
      print(model.errors, digits = digits, na.print = "-")
    }
    if(!is.null(x$predictions)) {
      cat("\nPredicted values:\n")
      with(x$predictions$model, print(predicted[k,], digits = digits))
      if(!is.null(x$bootstrap)) {
        txt <- paste("\nBootstrap-derived estimated values based on a",
                     ifelse(x$weighted, " weighted", ""),
                     " model with ", x$predictions$bootstrap$k,
                     "-closest analogues:", sep = "")
        cat("\n")
        writeLines(strwrap(txt))
        with(x$predictions$bootstrap, print(predicted[,k],
                                            digits = digits))
      }
    }
    invisible(x)
  }
