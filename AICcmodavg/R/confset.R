confset <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, method = "raw", level = 0.95, delta = 6, c.hat = 1) {

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    aic.table <- aictab(cand.set = cand.set, modnames = modnames, sort = TRUE, c.hat = c.hat,
                        second.ord = second.ord, nobs = nobs)

    ##add check to see whether response variable is the same for all models - these two lines not compatible with unmarked objects
    #check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    #if(length(unique(check.resp)) > 1) stop("You must use the same response variable for all models\n")
    
    nmods <- nrow(aic.table)
    
#method based on simply summing the Akaike weights until a given value is reached
    if(method=="raw") {
      iter <- 1
      sum.wt <- 0
      while(level >= sum.wt) {
        sum.wt <- aic.table[iter,6]+sum.wt
        if(sum.wt < level)    iter <- iter + 1
      }
      confset.tab <- list("method" = method, "level" = level, "table" = aic.table[1:iter,])
    }

#method based on ranking models ordinally according to their delta AIC values
    if (method=="ordinal") {
      substantial <- aic.table[which(aic.table[,4] <= 2),]
      some <- aic.table[which(aic.table[,4] > 2 & aic.table[,4] <= 7),]
      little <- aic.table[which(aic.table[,4] > 7 & aic.table[,4] <= 10), ]
      none <- aic.table[which(aic.table[,4] > 10), ]
      confset.tab <- list("method" = method, "substantial" = substantial, "some" = some, "little" = little, "none" = none)
    }

#method based on relative likelihoods for some cutoff value
    if (method=="ratio") {
      cutoff <- exp(-1*delta/2)
  #if a given cutoff value is requested, one can compute the corresponding delta AIC value as:
  #-2*log(cutoff); -2*log(0.125) = 4.16
      ratios <- matrix(NA, nrow=nmods, ncol=1)
      for (i in 1:nmods) {
        ratios[i] <- aic.table[i, 5]/aic.table[1, 5]
      }
      confset.tab <- list("method" = method, "cutoff" = cutoff, "delta" = delta,
                          "table" = aic.table[which(ratios> cutoff),])
    }

    class(confset.tab) <- c("confset", "list")
    return(confset.tab)
    
  }


print.confset <-
  function(x, digits = 2, ...) {
    if(x$method == "raw") {
      cat("\nConfidence set for the best model\n\n")
      cat("Method:\t", "raw sum of model probabilities\n\n")
      perc <- x$level*100
      perc.char <- paste(perc, "%", sep = "")
      cat(perc.char, "confidence set:\n")
      nice.tab <- cbind(x$table[, 2], x$table[, 3], x$table[, 4], x$table[, 6])
      colnames(nice.tab) <- colnames(x$table)[c(2, 3, 4, 6)]
      rownames(nice.tab) <- x$table[, 1]
      print(round(nice.tab, digits = digits)) #select rounding off with digits argument
      cat("\n")
      sum.wt <- round(sum(nice.tab[, 4]), digits = digits)
      cat("Model probabilities sum to", sum.wt, "\n\n")
    }

    if(x$method == "ordinal") {
      cat("\nConfidence set for the best model\n\n")
      cat("Method:\t", "ordinal ranking based on delta AIC\n\n")
      cat("Models with substantial weight:\n")
      nice.tab <- cbind(x$substantial[, 2], x$substantial[, 3], x$substantial[, 4], x$substantial[, 6])
      colnames(nice.tab) <- colnames(x$substantial)[c(2, 3, 4, 6)]
      rownames(nice.tab) <- x$substantial[,1]
      print(round(nice.tab, digits = digits)) #select rounding off with digits argument
      cat("\n\n")
      cat("Models with some weight:\n")
      nice.tab <- cbind(x$some[, 2],  x$some[, 3], x$some[, 4], x$some[, 6])
      colnames(nice.tab) <- colnames(x$some)[c(2, 3, 4, 6)]
      rownames(nice.tab) <- x$some[, 1]
      print(round(nice.tab, digits = digits)) #select rounding off with digits argument
      cat("\n\n")
      cat("Models with little weight:\n")
      nice.tab <- cbind(x$little[, 2], x$little[, 3], x$little[, 4], x$little[, 6])
      colnames(nice.tab) <- colnames(x$little)[c(2, 3, 4, 6)]
      rownames(nice.tab) <- x$little[, 1]
      print(round(nice.tab, digits = digits)) #select rounding off with digits argument
      cat("\n\n")
      cat("Models with no weight:\n")
      nice.tab <- cbind(x$none[, 2], x$none[, 3], x$none[, 4], x$none[, 6])
      colnames(nice.tab) <- colnames(x$none)[c(2, 3, 4, 6)]
      rownames(nice.tab) <- x$none[, 1]
      print(round(nice.tab, digits = digits)) #select rounding off with digits argument
      cat("\n\n")
    }
    
    if(x$method == "ratio") {
      cat("\nConfidence set for the best model\n\n")
      cat("Method:\t", "ranking based on relative model likelihood\n\n")
      round.cut <- round(x$cutoff, digits = digits)
      cat("Cutoff value:\t", round.cut, "(corresponds to delta AIC of", paste(x$delta,")", sep = ""),"\n\n")
      cat("Confidence set for best model:\n")
      nice.tab <- cbind(x$table[, 2], x$table[, 3], x$table[, 4], x$table[, 6])
      colnames(nice.tab) <- colnames(x$table)[c(2, 3, 4, 6)]
      rownames(nice.tab) <- x$table[, 1]
      print(round(nice.tab, digits = digits)) #select rounding off with digits argument
      cat("\n\n")
    }
  }


