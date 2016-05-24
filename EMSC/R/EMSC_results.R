#' Predict Method for EMSC
#' 
#' Prediction for \code{EMSC} ojects. Corrections are calculated for the new
#' \code{matrix} based on the EMSC model used in the input object.
#' 
#' @param object An object fitted by the \code{EMSC} function.
#' @param newdata A \code{matrix} or object convertable to a matrix containing observations as rows.
#' @param ... unused.
#' 
#' @seealso \code{\link{EMSC}} \code{\link{EMSC_model}}
#' 
#' @examples
#' data(milk)
#' Raman.cal <- milk$Raman[  1:90,  850:3300]
#' Raman.val <- milk$Raman[-(1:90), 850:3300]
#' EMSC.cal  <- EMSC(Raman.cal)
#' EMSC.val  <- predict(EMSC.cal, Raman.val)
#' identical(EMSC.cal$model, EMSC.val$model) # Same model, reference spectrum, etc.
#' 
#' matplot(t(EMSC.cal$corrected), type = 'l', col = 'black', lty = 1, ylab = 'Intensity')
#' matplot(t(EMSC.val$corrected), type = 'l', col = 'red', lty = 2, add = TRUE)
#' legend('topleft', legend = c('Calibration','Validation'), lty = 1:2, col = 1:2)
#' 
#' @importFrom stats predict
#' @importFrom graphics axTicks matplot
#' @export
predict.EMSC <- function(object, newdata = NULL, ...){
  if(is.null(newdata))
    return(object)
  
  newdata <- unclass(as.matrix(newdata))
  EMSC(newdata, object$model)
}


#' @aliases print.EMSC summary.EMSC
#' 
#' @title Plot, print and summary methods for EMSC
#' 
#' @description Plotting routine for \code{EMSC} ojects. The default behaviour is to
#' plot raw spectra, reference spectrum, polynomials, interferents,
#' constituents, replicate model and corrected spectra. This can be
#' tweaked by changing the parameters.
#' 
#' @param object An object fitted by the \code{EMSC} function.
#' @param x An object fitted by the \code{EMSC} function.
#' @param y Unused parameter to conform to generic \code{plot}.
#' @param spec Parameter specifying if all spectra should be plotted (default)
#' or a subset (\code{numeric} vector).
#' @param what \code{character} vector defining what to plot.
#' @param where \code{integer} vector defining which elements should be plotted
#' in which subplot.
#' @param revX Reverse x axis (default = FALSE).
#' @param labels "names" or "numbers" uses column names for x axis labelling.
#' @param type plotting type (line, points, ...).
#' @param lty line type.
#' @param lwd line width.
#' @param pch plot character.
#' @param cex symbol/line scaling.
#' @param col symbol/line colour.
#' @param xlab x label.
#' @param ylab y label.
#' @param pretty.xlabels Use pretty x labels (default = TRUE).
#' @param xlim x limits.
#' @param ... Additional arguments to \code{matplot}.
#' 
#' @details The parameters \code{what} and \code{where} must match
#' so that the parts of the EMSC model end up in the correct subplot.
#' There are limits to the freedom of this function.
#' 
#' \code{print} and \code{summary} return minimal information on the
#' EMSC object.
#' 
#' @return No return.
#' 
#' @author Kristian Hovde Liland
#' 
#' @seealso \code{\link{EMSC}} \code{\link{EMSC_model}} \code{\link{plot.EMSC}}
#' 
#' @examples
#' data(milk, package = "EMSC")
#' Raman      <- milk$Raman[, 850:3300]
#' EMSC.rep   <- EMSC(Raman, degree = 6, reference = Raman[30, ],
#'                    replicates = milk$replicates)
#' plot(EMSC.rep)
#'                    
#' @importFrom graphics axis box grid lines par plot points polygon text
#' @export
plot.EMSC <- function(x, y, spec = "all", 
                      what = c("raw", "reference", "polynomials", "interferents",
                               "constituents", "replicates", "corrected"),
                      where = c(1,1, 2, 3,3, 4, 5),
                      revX = FALSE, labels, type, lty, lwd = NULL,
                      pch, cex = NULL, col, xlab, ylab, 
                      pretty.xlabels = TRUE, xlim, ...){
  if (!missing(y))
    warning("Argument 'y' is ignored")
  
  # Which spectra
  if(!is.numeric(spec) && spec == "all")
    spec <- 1:x$model$sizes[1]
  n <- dim(x$X)[1]
  
  # Abscissa
  xnum <- 1:x$model$sizes[1]
  xnam <- colnames(x$X)
  if (missing(type)) 
    type <- "l"
  if (missing(lty)) 
    lty <- 1:n
  if (missing(pch)) 
    pch <- 1:n
  if (missing(col)) 
    col <- 1:n
  if (missing(xlab)) 
    xlab <- "variable"
  if (missing(ylab)) 
    ylab <- "Relativ intensity"
  if(missing(labels)){
    xaxt <- par("xaxt")
    # xl <- 1:x$model$sizes[1]
  } else {
    xaxt <- "n"
    switch(match.arg(labels, c("names", "numbers")), 
           names = {
             labels <- xnam
           }, numbers = {
             if (length(grep("^[-0-9.]+[^0-9]*$", xnam)) == 
                 length(xnam)) {
               labels <- sub("[^0-9]*$", "", xnam)
               if (isTRUE(pretty.xlabels)) {
                 xnum <- as.numeric(labels)
                 xaxt <- par("xaxt")
               }
             } else {
               stop("Could not convert variable names to numbers.")
             }
           })
    if (!missing(labels) && xaxt == "n") {
      if (isTRUE(pretty.xlabels)) {
        ticks <- axTicks(1)
        ticks <- ticks[ticks >= 1 & ticks <= length(labels)]
      }
      else {
        ticks <- 1:length(labels)
      }
      axis(1, ticks, labels[ticks], ...)
    }
  }
  
  # Check what can and should be plotted
  terms <- rownames(x$model$model)
  k  <- 1
  po <- 1
  pos <- list()
  if("raw" %in% what){
    pos[[po]] <- TRUE
    po <- po+1
  }
  if("reference" %in% what){
    pos[[po]] <- 1
    k  <- k+1
    po <- po+1
  }
  if("polynomials" %in% what){
    if((p <- x$model$sizes[2]) > 0){
      pos[[po]] <- k:(k+p)
      k  <- k+p+1
      po <- po+1
    } else { # Remove if not present
      where <- where[-match("polynomials", what)]
      what  <- what [-match("polynomials", what)]
    }
  }
  if("interferents" %in% what){
    if((p <- x$model$sizes[3]) > 0){
      pos[[po]] <- k:(k+p-1)
      k  <- k+p+1
      po <- po+1
    } else { # Remove if not present
      where <- where[-match("interferents", what)]
      what  <- what [-match("interferents", what)]
    }
  }
  if("replicates" %in% what){
    if((p <- x$model$sizes[5]) > 0){
      pos[[po]] <- k:(k+p-1)
      k  <- k+p+1
      po <- po+1
    } else { # Remove if not present
      where <- where[-match("replicates", what)]
      what  <- what [-match("replicates", what)]
    }
  }
  if("constituents" %in% what){
    if((p <- x$model$sizes[4]) > 0){
      pos[[po]] <- k:(k+p-1)
      k  <- k+p+1
      po <- po+1
    } else { # Remove if not present
      where <- where[-match("constituents", what)]
      what  <- what [-match("constituents", what)]
    }
  }
  if("corrected" %in% what){
    pos[[po]] <- TRUE
  }
  
  # Set up plotting
  lu <- length(u <- unique(where))
  if(lu <= 3){
    old.par <- par(mfrow = c(lu,1))
  } else {
    nr <- ceiling(sqrt(lu))
    nc <- ceiling(lu/nr)
    old.par <- par(mfrow = c(nr, nc))
  }
  
  # Plot
  if (missing(xlim)) 
    xlim <- xnum[c(1, length(xnum))]
  if(revX)
    xlim <- rev(xlim)
  for(i in 1:lu){
    ind <- which(where == u[i])
    if(any(what[ind] == "raw")){  # Plot raw spectra
      matplot(xnum, t(x$X), xlab = xlab, ylab = ylab, type = type, 
              lty = lty, lwd = lwd, pch = pch, cex = cex, col = col, 
              xaxt = xaxt, xlim = xlim, main = "Raw spectra", ...)
      if (!missing(labels) && xaxt == "n") {
        if (isTRUE(pretty.xlabels)) {
          ticks <- axTicks(1)
          ticks <- ticks[ticks >= 1 & ticks <= length(labels)]
        }
        else {
          ticks <- 1:length(labels)
        }
        axis(1, ticks, labels[ticks], ...)
      }
      # Add reference?
      if(any(what[ind] == "reference") && where[which(what[ind] == "raw")] == where[re <- which(what[ind] == "reference")]){
        matplot(xnum, x$model$model[pos[[re]],], add = TRUE, lwd = 3)
      }
    } else {
      # Only reference
      if(any(what[ind] == "reference")){
        matplot(xnum, x$model$model[pos[[which(what == "reference")]],], xlab = xlab, ylab = ylab, type = type, 
                lty = lty, lwd = lwd, pch = pch, cex = cex, col = col, 
                xaxt = xaxt, xlim = xlim, main = "Reference spectrum", ...)
        if (!missing(labels) && xaxt == "n") {
          if (isTRUE(pretty.xlabels)) {
            ticks <- axTicks(1)
            ticks <- ticks[ticks >= 1 & ticks <= length(labels)]
          }
          else {
            ticks <- 1:length(labels)
          }
          axis(1, ticks, labels[ticks], ...)
        }
      }
    }
    if(any(what[ind] == "polynomials")){
      matplot(xnum, t(x$model$model[pos[[which(what == "polynomials")]],]), xlab = xlab, ylab = ylab, type = type, 
              lty = lty, lwd = lwd, pch = pch, cex = cex, col = col, 
              xaxt = xaxt, xlim = xlim, main = "Polynomials", ...)
      if (!missing(labels) && xaxt == "n") {
        if (isTRUE(pretty.xlabels)) {
          ticks <- axTicks(1)
          ticks <- ticks[ticks >= 1 & ticks <= length(labels)]
        }
        else {
          ticks <- 1:length(labels)
        }
        axis(1, ticks, labels[ticks], ...)
      }
    }
    if(any(what[ind] == "interferents")){
      matplot(xnum, t(x$model$model[pos[[which(what == "interferents")]],]), xlab = xlab, ylab = ylab, type = type, 
              lty = lty, lwd = lwd, pch = pch, cex = cex, col = col, 
              xaxt = xaxt, xlim = xlim, main = "Interferents", ...)
      if (!missing(labels) && xaxt == "n") {
        if (isTRUE(pretty.xlabels)) {
          ticks <- axTicks(1)
          ticks <- ticks[ticks >= 1 & ticks <= length(labels)]
        }
        else {
          ticks <- 1:length(labels)
        }
        axis(1, ticks, labels[ticks], ...)
      }
    }
    if(any(what[ind] == "constituents")){
      matplot(xnum, t(x$model$model[pos[[which(what == "constituents")]],]), xlab = xlab, ylab = ylab, type = type, 
              lty = lty, lwd = lwd, pch = pch, cex = cex, col = col, 
              xaxt = xaxt, xlim = xlim, main = "Constituents", ...)
      if (!missing(labels) && xaxt == "n") {
        if (isTRUE(pretty.xlabels)) {
          ticks <- axTicks(1)
          ticks <- ticks[ticks >= 1 & ticks <= length(labels)]
        }
        else {
          ticks <- 1:length(labels)
        }
        axis(1, ticks, labels[ticks], ...)
      }
    }
    if(any(what[ind] == "replicates")){
      matplot(xnum, t(x$model$model[pos[[which(what == "replicates")]],]), xlab = xlab, ylab = ylab, type = type, 
              lty = lty, lwd = lwd, pch = pch, cex = cex, col = col, 
              xaxt = xaxt, xlim = xlim, main = "Replicate model", ...)
      if (!missing(labels) && xaxt == "n") {
        if (isTRUE(pretty.xlabels)) {
          ticks <- axTicks(1)
          ticks <- ticks[ticks >= 1 & ticks <= length(labels)]
        }
        else {
          ticks <- 1:length(labels)
        }
        axis(1, ticks, labels[ticks], ...)
      }
    }
    if(any(what[ind] == "corrected")){
      matplot(xnum, t(x$corrected), xlab = xlab, ylab = ylab, type = type, 
              lty = lty, lwd = lwd, pch = pch, cex = cex, col = col, 
              xaxt = xaxt, xlim = xlim, main = "Corrected spectra", ...)
      if (!missing(labels) && xaxt == "n") {
        if (isTRUE(pretty.xlabels)) {
          ticks <- axTicks(1)
          ticks <- ticks[ticks >= 1 & ticks <= length(labels)]
        }
        else {
          ticks <- 1:length(labels)
        }
        axis(1, ticks, labels[ticks], ...)
      }
    }
  }
  
  par(old.par)
  return(invisible(NULL))
}

#' @rdname plot.EMSC
#' @export
print.EMSC <- function(x, ...){
  N <- dim(x$corrected)
  cat("EMSC object")
  cat("\n", N[1], " spectra of length ", N[2], sep = "")
  sizes <- x$model$sizes[-1]
  cat("\nPolynomial degree: ", sizes[1], sep = "")
  if(sizes[2]>0)
    cat("\nInterferents: ", sizes[2], sep = "")
  if(sizes[3]>0)
    cat("\nConstituents: ", sizes[3], sep = "")
  if(sizes[4]>0)
    cat("\nReplicate dimensions: ", sizes[4], sep = "")
}

#' @rdname plot.EMSC
#' @export
summary.EMSC <- function(object, ...){
  print(object)
  cat("\nModel terms:\n")
  cat(rownames(object$parameters))
}  
