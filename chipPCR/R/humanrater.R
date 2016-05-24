humanrater <-
  function(x, cyc = 1, repeats = 1, 
           designations = list(y = "yes", a = "ambiguous", n = "not"), 
           shuffle = TRUE, ...) {
    if (!is.numeric(repeats) || length(repeats) > 1)
      stop("'repeats' must be a numeric vector of length 1.")
    if (repeats < 1)
      stop("'repeats' must be have value 1 or bigger")
    if (length(designations) < 2)
      stop("'designations' must have length at least 2.")
    if (!is.list(designations))
      stop("'designations' must be a list.")
    #     if (!is.character(designations) || length(designations) != 3)
    #       stop("'designations' must be a character vector of length 3.")
    allowed.symbols <- names(designations)
    
    if (any(table(allowed.symbols) > 1))
      stop("Do not use repeated names.")
    if ("" %in% allowed.symbols)
      stop("All elements of 'designations' list must be named.")
    prompt.line <- paste(sapply(1L:length(designations), function(i)
      paste0("[", allowed.symbols[i], "] if ", designations[[i]])), collapse = ", ")
    
    all.ratings <- sapply(1L:repeats, function(j) {
      if(repeats > 1)
        cat(paste0("Repeat:", j, "\n"))
      rating.order <- 2L:ncol(x)
      if(shuffle)
        rating.order <- sample(rating.order)
      all.ratings <- sapply(rating.order, function(i) {
        
        plot(x[, 1], x[, i], main = paste0("Experiment ", i), type = "b", pch = 19, lwd = 2, 
             xlab = "Cycle", ylab = "Fluorescence", ...)
        #declare dummy variable - without it while loop does not work
        #the dummy cannot belong to the set of designations
        res <- "dummy which surely would not be a designation"
        while((!(res %in% allowed.symbols))) {
          res <- readline(prompt = prompt.line)
          #check correctness
          if (!res %in% allowed.symbols)
            cat("Wrong input. Try again.\n")
        }
        res
      })
      all.ratings[order(rating.order)]
    })
    if (repeats > 1) {
      #check conformity of a human assessment
      check <- apply(all.ratings, 1, function(k) {
        ifelse(length(unique(k)) == 1, TRUE, FALSE)
      }) 
    } else {
      check <- rep("not tested", nrow(all.ratings))
    }
    data.frame(test.result = all.ratings, conformity = check)
  }
