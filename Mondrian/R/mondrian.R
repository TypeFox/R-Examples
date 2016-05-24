mondrian <- function(data, labels = colnames(data), xlab = "", ylab = "" , main = "", col = NULL , pop = NULL, indiv = FALSE, ...) {
  
  ## Initial checking
  data <- as.data.frame(data)
  data <- data[complete.cases(data), ]  ## delete rows (individuals) with missing data
  
  ## Default values
  if(is.null(labels))
    labels <- colnames(data)
  
  if(is.null(col)) {
    if(length(labels) == 2) 
      col <- c("blue", "red")
    else
      col <- brewer.pal(length(labels), "Set1")
  } else { 
    col <- rep(col, length.out = length(labels))
  }
  
  main <- paste(main, " (n = ", nrow(data), ")", sep = "")
  
  if(!is.null(pop)) { ## Individuals are grouped in sub-populations
    
    labelpop <- levels(unique(data[, pop]))
    nbpop <- length(labelpop)
    
    ## Graphic window management 
    nrow <- floor(sqrt(nbpop + 1))
    par(mfrow = c(nrow, ceiling(nbpop / nrow)))
    
    outpop <- list()
    ## Results for each sub-population
    subpop <- by(data, data[, pop], function(x) mondrian(x[, - pop], pop = NULL, xlab = xlab, ylab = ylab, main = unique(x[, pop]), col = col, indiv = indiv, ...))
    outpop <- lapply(subpop, function(x) x)
    names(outpop) <- labelpop
    
    ## Result if all individuals belong to the same population
    outpop$pop <- mondrian(data[, - pop], xlab = xlab, ylab = ylab , main = "Total population", col = col, pop = NULL, indiv = indiv, ...)
    par(mfrow = c(1, 1))
    invisible(outpop)
    
    
  } else {  ## Individuals are defined as belonging to the same population
    
    ## Percents matrix building
    counts_profiles <- rev(table(apply(data, 1, paste, collapse = ""))) ## counts table for profiles
    percents_profiles <- counts_profiles / sum(counts_profiles)  ## percents table for profiles
    mat_profiles <- matrix(rep(percents_profiles, length(col)), byrow = FALSE, ncol = length(col), nrow = length(percents_profiles))
    dimnames(mat_profiles) <- list(names(percents_profiles), labels)
    
    ## Presence-absence profiles
    profiles <- t(data.frame(strsplit(rownames(mat_profiles), "")))
    
    ## Display
    par(las = 3, mar = c(4, 3, 2.5, 1), mgp = c(2, 0, 0), font = 2, cex.axis = 1.2)
    plot(ncol(mat_profiles), 1, xlim = c(0, ncol(mat_profiles)), ylim = c(0, 1), type = "n", axes = FALSE, xlab = xlab, ylab = ylab, main = main, ...)
    axis(2, at = c(0, 1), c("0", "100%"), tick = FALSE, lwd = 0.5)
    par(las = 1, font = 2, cex.axis = 0.9)
    axis(1, at = 0.5:ncol(mat_profiles), labels, tick = FALSE)
    
    dimrect <- rbind(rep(0, ncol(mat_profiles)), apply(mat_profiles, 2, cumsum))
    
    out <- sapply(1:ncol(mat_profiles), function(Y) {
      sapply(1:nrow(mat_profiles), function(X) {
        rect(Y - 1, dimrect[X], Y, dimrect[X + 1], col = ifelse(profiles[X, Y] == "1", col[Y], "white"))
      })
    })
    
    if(indiv)
      sapply(1:nrow(data), function(i) lines(0:ncol(data), rep(i/nrow(data), 1 + ncol(data))))
    
    invisible(percents_profiles)
  }
}
