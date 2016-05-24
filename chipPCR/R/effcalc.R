effcalc <- function(x, y, logx = TRUE, RSD = FALSE, rob = FALSE, level = 0.95) {
  testxy(x, y, txt.x = "Enter Concentration", txt.y = "Enter Cq data", 
         length = FALSE)
  
  if(!is.matrix(y)) 
    y <- as.matrix(y)
  
  # Removing all NA Cq rows
  i <- apply(y, 1, function(yrow) all(is.na(yrow)))
  if (TRUE %in% i) {
    x <- x[!i]
    y <- y[!i, ]
    warning(
      sprintf("Row %i was removed because it does not contain Cq data.\n",
              which(i == TRUE)))
  }
  
  if (logx) {
    x.tmp <- log10(x)
  } else {
    x.tmp <- x
  }
  
  if (rob) {
    loc.fct <- median
    dev.fct <- mad
  } else {
    loc.fct <- mean
    dev.fct <- sd
  }
  
  if (ncol(data.frame(y)) > 1) {
    y.m <- apply(y, 1, function(x) loc.fct(x, na.rm = TRUE))
    y.sd <- apply(y, 1, function(x) dev.fct(x, na.rm = TRUE))
  } else {
    y.m <- y
    y.sd <- rep(0, length(y))
  }
  
  if (RSD) {
    y.cv <- (y.sd / y.m) * 100
  } else {
    y.cv <- y.sd / y.m
  }
  
  # Aggregate calculated data and check for consistency
  # of the concentration (check if Na of Inf values
  # were produce and exclude these from further
  # calculation, check if at least two values for the
  # the linear regression are present.).
  res <- data.frame(x.tmp, y.m, y.sd, y.cv)
  res <- res[which(is.finite(res[, 1])), ]
  
  if (nrow(res) < 2) {
    stop("Cannot perform calculation. At least two
	  dilutions required.")
  }
  if (nrow(res) == 2) {
    warning("At least three dilutions should be used to 
	    determine an amplificantion efficiency.")
  }
  
  # Decide which type of measure (e.g., mean vs. median, 
  # standard deviation vs. standard error) are shown in
  # the plot.
  
  names(res) <- c("Concentration", "Location (Mean)", 
                  "Deviation (SD)", 
                  "Coefficient of Variance (RSD [%])")
  if (rob)
    names(res)[2L:3] <- c("Location (Median)", 
                          "Deviation (MAD)")
  
  if (RSD)
    names(res)[4] <- "Coefficient of Variance (RSD)"
  
  # Perform a linear regression based on the values of the 
  # calculated mean/median
  lm.res <- lm(res[, 2] ~ res[, 1])
  # Calculate goodness of fit
  
  # Calculate the amplification efficiency (in percent)
  AE <- round(10^(-1/coef(lm.res)[2])/ 2 * 100, 1)
  
  # Calculate correlation between the concentration and 
  # the Cq values along with the significance level
  cortest <- cor.test(res[, 1], res[, 2], conf.level = level)
  
  cortest[["data.name"]] <- paste0("~ `", names(res)[2], "` and ", names(res)[1])
  
  # TO DO: change call to nicer format  
  #   lm.res[["call"]] <- paste0("lm(formula = ", 
  #                              paste0("`", names(res)[2], "` ~ ", names(res)[1]),
  #                              ", data = res)")
  
  new("eff", .Data = data.matrix(res),
      amplification.efficiency = AE,
      regression = lm.res, correlation.test = cortest)
}

setGeneric("effcalc")

setMethod("effcalc", signature(x = "data.frame", y="missing"), 
          function(x, y, logx = TRUE, RSD = FALSE, rob = FALSE, level = 0.95) { 
            if (ncol(x) != 2) 
              stop("'x' must have two columns.")
            effcalc(x[, 1], x[, 2], logx = TRUE, RSD = FALSE, rob = FALSE,
                    level = 0.95)
          })

setMethod("effcalc", signature(x = "matrix", y="missing"), 
          function(x, y, logx = TRUE, RSD = FALSE, rob = FALSE, level = 0.95) { 
            if (ncol(x) != 2) 
              stop("'x' must have two columns.")
            effcalc(x[, 1], x[, 2], logx = TRUE, RSD = FALSE, rob = FALSE,
                    level = 0.95)
          })
