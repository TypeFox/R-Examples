excessRisk <-
function (form = weight ~ seatbelt + airbag, response = "dead", 
            margin = "airbag", data = DAAG::nassCDS, decpl = 4, printResults = TRUE) 
{
  vars <- all.vars(form)
  isweight <- length(form) > 2
  deadweight <- data[, response]
  resplev <- levels(deadweight)
  if (is.null(resplev)) 
    resplev <- paste(0:1)
  if (is.factor(deadweight)) 
    deadweight <- (unclass(deadweight) - 1)
  if (!all(deadweight %in% 0:1)) 
    stop(paste("The parameter 'response' must either be", 
               "a 2-level factor or a 0/1 variable"))
  if (isweight) {
    deadweight <- deadweight * data[, vars[1]]
    rhs <- deparse(form[[3]])
  }
  else rhs <- deparse(form[[2]])
  formdead <- formula(paste("deadweight", " ~ ", rhs, sep = ""))
  total <- with(data, as.data.frame(xtabs(form, data = data)))
  dead <- with(data, as.data.frame(xtabs(formdead, data = data)))
  nc <- match("Freq", names(total))
  nway <- nc - 1
  nair <- match(margin, names(total))
  marglev <- levels(total[, margin])
  if (is.null(marglev)) 
    marglev <- sort(unique(total[, margin]))
  marglab <- as.character(marglev)
  lev1 <- with(total, total[, margin] == marglev[1])
  lev2 <- with(total, total[, margin] == marglev[2])
  nobag_d <- paste(marglab[1], "_", resplev[2], sep = "")
  nobag_tot <- paste(marglab[1], "_tot", sep = "")
  bag_d <- paste(marglab[2], "_", resplev[2], sep = "")
  bag_tot <- paste(marglab[2], "_tot", sep = "")
  nobagProp <- paste(marglab[1], "Prop", sep = "")
  bagProp <- paste(marglab[2], "Prop", sep = "")
  df <- cbind(dead[lev1, (1:nway)[-nair], drop = FALSE], nobag_d = dead[lev1, 
                                                           nc], nobag_tot = total[lev1, nc], bag_d = dead[lev2, 
                                                                                               nc], bag_tot = total[lev2, nc], nobagProp = dead[lev1, 
                                                                                                                                 nc]/total[lev1, nc], bagProp = dead[lev2, nc]/total[lev2, 
                                                                                                                                                        nc])
  names(df)[nway:(nway + 5)] <- c(nobag_d, nobag_tot, bag_d, 
                                  bag_tot, nobagProp, bagProp)
  if (printResults) {
    printdf <- df
    numcols <- c(nobag_d, nobag_tot, bag_d, bag_tot)
    fraccols <- c(nobagProp, bagProp)
    printdf[, numcols] <- round(printdf[, numcols])
    printdf[, fraccols] <- round(printdf[, fraccols], decpl)
    diffdf <- data.frame(df[, bag_d] - df[, bag_tot] * df[, nobagProp])  
    dimnames(diffdf) <- list(row.names(printdf),
                             paste("Difference:", bag_d, "-", nobag_d))
    print(printdf)
    cat("\n")
    print(diffdf)
    cat("Differences in expected number of deaths are calculated",
        paste("\nrelative to the level ",
              "'", marglev[1], "'", " of the factor ", "'",
              margin, "'", ".", sep=""), "\n")
  }
  invisible(df)
}

