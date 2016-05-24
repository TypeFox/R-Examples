# MHLC scale scoring script
# Creation date: 2013-06-20
# --------------

testChar <- list(acronym = "MHLC",
                 name = "Multidimensional Health Locus of Control",
                 ref = "Wallston KA, Wallston BS & DeVelis R, 1978",
                 n.items = 18,
                 valid = c(1, 2, 3, 4, 5, 6),
                 miss = c(0),
                 comm = "Public domain: www.nursing.vanderbilt.edu/faculty/kwallston/mhlcscales.htm")

scoring.fun <- function(answers, sex, age, id, date.test, comm) {
  # "answer" is a *character* vector as introduced through the keyboard.
  # "sex", "age", "id", "date.test" & "comm" used if appropiate.
  answers <- as.numeric(answers) # transform to numeric for easier scoring
  answers[answers %in% c(0)] <- NA # missing characters to NA
  blanks <- sum(is.na(answers)) # compute number of missing itemss
  pcnt.blanks <- round((blanks / 18) * 100) # compute % of missing items
  results <- data.frame(NULL) # Initialize null data frame for results

  # I scale scoring commands
  # --------------
  results[1, "Acronym"] <- "I" # acronym
  results[1, "Scale"] <- "Internal" # name of the scale
  # Items making up the scale
  items <- c(1, 6, 8, 12, 13, 17)
  results[1, "Miss"] <- sum(is.na(answers[items])) # number of missing items
  if (results[1, "Miss"] == length(items)) {  # all answers are missing items
    results[1, "Raw"] <- NA
  } else {
    results[1, "Raw"] <- sum(answers[items], na.rm=TRUE) # sum answered items
  }

  # C scale scoring commands
  # --------------
  results[2, "Acronym"] <- "C" # acronym
  results[2, "Scale"] <- "Chance" # name of the scale
  # Items making up the scale
  items <- c(2, 4, 9, 11, 15, 16)
  results[2, "Miss"] <- sum(is.na(answers[items])) # number of missing items
  if (results[2, "Miss"] == length(items)) {  # all answers are missing items
    results[2, "Raw"] <- NA
  } else {
    results[2, "Raw"] <- sum(answers[items], na.rm=TRUE) # sum answered items
  }

  # P scale scoring commands
  # --------------
  results[3, "Acronym"] <- "P" # acronym
  results[3, "Scale"] <- "Powerful Others" # name of the scale
  # Items making up the scale
  items <- c(3, 5, 7, 10, 14, 18)
  results[3, "Miss"] <- sum(is.na(answers[items])) # number of missing items
  if (results[3, "Miss"] == length(items)) {  # all answers are missing items
    results[3, "Raw"] <- NA
  } else {
    results[3, "Raw"] <- sum(answers[items], na.rm=TRUE) # sum answered items
  }

  # Vector for writing scores to a file
  # --------------------
  results.scores <- unlist(results[-c(1, 2)])  # not Acronym & Scale columns
  names <- paste(results$Acronym, names(results.scores), sep=".")
  names(results.scores) <- sub("[0-9]+$", "", names)  # delete ending numbers

  # Output in form of list
  # ------------------
  results.lst <- list(paste("Total number of missing items: ",
                            blanks,
                            " (",
                            pcnt.blanks,
                            "%)",
                            sep=""))

  # Return results
  # ------------------
  list(results.lst = results.lst,
              results.df = results,
              results.scores = results.scores)

} # end of scoring.fun