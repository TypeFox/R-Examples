cover.plot <- function(object, row, column, x = NULL, CI = 0.95,
                      medians = TRUE, col = NULL, ylim = c(0,1), 
                      ylab, lty = par("lty"), lwd = par("lwd"), ...) { 

  if (any(class(object) == "eiMD")) {
    Beta <- object$draws$Beta
    if (is.null(Beta)) 
      stop("precinct-level draws for Beta not stored in output from ei.MD.bayes")
  }
  else if (any(class(object) == "eiMD.beta")) {
    Beta <- object
  }
  else
     stop("coverplot only works with output from ei.MD.bayes or read.beta")
  if (is.mcmc(Beta)) { 
    tnames <- strsplit(colnames(Beta), "beta.")
    get2 <- function(x) x[2]
    idx <- strsplit(sapply(tnames, get2), ".", fixed = TRUE)
    idx <- as.list(as.data.frame(matrix(unlist(idx), byrow = TRUE,
                                        nrow = length(idx), ncol = length(idx[[1]]))))
    idx <- lapply(idx, as.character)
    idx <- lapply(idx, unique)
    idx[[4]] <- 1:length(Beta[,1])
    names(idx) <- c("rows", "columns", "precincts", "sims")
    Betas <- array(t(object$draws$Beta), dim = sapply(idx, length),
                 dimnames = idx)
  }
  else {
    idx <- dimnames(Beta)
  }

  if (missing(row) | missing(column)) {
    stop("you must select a row and column marginal.")
  }
  frow <- idx$rows[grep(row, idx$rows)]
  fcolumn <- idx$columns[grep(column, idx$columns)]
  if (length(frow) < 1)
    stop(paste(row, "not among available row marginals:",
               paste(idx$rows, collapse = " ")))
  if (length(frow) > 2)
    stop(paste(row, "matches more than one row marginal"))
  if (length(fcolumn) < 1)
    stop(paste(column, "not among available column marginals:",
               paste(idx$columns, collapse = " ")))
  if (length(fcolumn) > 2)
    stop(paste(column, "matches more than one column marginal"))
  
  usebetas <- t(Betas[frow, fcolumn, ,])
  quant.fcn <- function(x){quantile(x, c(0.5 - CI/2, 0.5 + CI/2))}
  seglims <- apply(usebetas, 2, quant.fcn)
  meds <- apply(usebetas, 2, median)

  if (is.null(x)) 
    x <- as.integer(names(meds))

  if (missing(ylab))
    ylab <- paste("Proportion of", row, "in", column)

  if (length(col) >= 2) {
    col1 <- col[1]
    col2 <- col[2]
  }
  else { 
    if (is.null(col))  
      col1 <- col2 <- par("fg")
    else 
      col1 <- col2 <- col
  }
  if (medians) type <- "p"
  else type <- "n"

  plot(x, meds, type = type, ylim = ylim, ylab = ylab, col = col1, ...)
  segments(x, seglims[1,], x, seglims[2,], col = col2, lty = lty, lwd = lwd)
}
  
