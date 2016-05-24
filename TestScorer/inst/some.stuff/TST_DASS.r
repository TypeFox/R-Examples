# DASS scale scoring script
# Creation date: 2013-06-20
# --------------

testChar <- list(acronym = "DASS",
                 name = "Depression Anxiety Stress Scales",
                 ref = "Lovibond SH & Lovibond PF, 1995",
                 n.items = 42,
                 valid = c(0, 1, 2, 3),
                 miss = c(4),
                 comm = "Public domain: www2.psy.unsw.edu.au/dass")

scoring.fun <- function(answers, sex, age, id, date.test, comm) {
  # "answer" is a *character* vector as introduced through the keyboard.
  # "sex", "age", "id", "date.test" & "comm" used if appropiate.
  answers <- as.numeric(answers) # transform to numeric for easier scoring
  answers[answers %in% c(4)] <- NA # missing characters to NA
  blanks <- sum(is.na(answers)) # compute number of missing items
  pcnt.blanks <- round((blanks / 42) * 100) # compute % of missings
  results <- data.frame(NULL) # Initialize null data frame for results

  # D scale scoring commands
  # --------------
  results[1, "Acronym"] <- "D" # acronym
  results[1, "Scale"] <- "Depression" # name of the scale
  # Items making up the scale
  items <- c(3, 5, 10, 13, 16, 17, 21, 24, 26, 31, 34, 37, 38, 42)
  results[1, "Miss"] <- sum(is.na(answers[items])) # number of missing items
  if (results[1, "Miss"] == length(items)) {  # all answers are missings
    results[1, "Raw"] <- NA
    results[1, "Centil"] <- NA
  } else {
    results[1, "Raw"] <- sum(answers[items], na.rm=TRUE) # sum answered items
    # Vector for score transformation
    trans.table <- c("5-20", "25-35", "40-45", "50-55", "60", "65", "70", "75", "76-78", "79-81", "82-83", "84-85", "86", "87", "88-89", "90", "91", "92", "93", "93", "94", "94", "95", "95", "96", "96", "96", "97", "97", "97", "97", "98", "98", "98", "98", "98", "99", "99", "99", "99", "99", "99", "99")
    index <- results[1, "Raw"]
    index <- 1 + index  # 1, because raw scores begin at 0
    results[1, "Centil"] <- trans.table[index]
  }

  # A scale scoring commands
  # --------------
  results[2, "Acronym"] <- "A" # acronym
  results[2, "Scale"] <- "Anxiety" # name of the scale
  # Items making up the scale
  items <- c(2, 4, 7, 9, 15, 19, 20, 23, 25, 28, 30, 36, 40, 41)
  results[2, "Miss"] <- sum(is.na(answers[items])) # number of missing items
  if (results[2, "Miss"] == length(items)) {  # all answers are missings
    results[2, "Raw"] <- NA
    results[2, "Centil"] <- NA
  } else {
    results[2, "Raw"] <- sum(answers[items], na.rm=TRUE) # sum answered items
    # Vector for score transformation
    trans.table <- c("5-30", "35-45", "50-55", "60-65", "70-75", "76-79", "80-83", "84-86", "87-89", "90", "91", "92", "93", "94", "94", "95", "95", "96", "96", "96", "97", "97", "98", "98", "98", "98", "99", "99", "99", "99", "99", "99", "99", "99", "99", "99", "99", "99", "99", "99")
    index <- results[2, "Raw"]
    index <- 1 + index  # 1, because raw scores begin at 0
    results[2, "Centil"] <- trans.table[index]
  }

  # S scale scoring commands
  # --------------
  results[3, "Acronym"] <- "S" # acronym
  results[3, "Scale"] <- "Stress" # name of the scale
  # Items making up the scale
  items <- c(1, 6, 8, 11, 12, 14, 18, 22, 27, 29, 32, 33, 35, 39)
  results[3, "Miss"] <- sum(is.na(answers[items])) # number of missing items
  if (results[3, "Miss"] == length(items)) {  # all answers are missings
    results[3, "Raw"] <- NA
    results[3, "Centil"] <- NA
  } else {
    results[3, "Raw"] <- sum(answers[items], na.rm=TRUE) # sum answered items
    # Vector for score transformation
    trans.table <- c("5", "10", "15", "20-25", "30", "35", "40", "45", "50-55", "60", "65", "65", "70", "75-77", "78-80", "81-82", "83-84", "85-86", "87-88", "89", "90", "91", "92", "93", "93", "94", "95", "95", "96", "96", "97", "97", "97", "97", "98", "98", "98", "99", "99", "99", "99", "99", "99")
    index <- results[3, "Raw"]
    index <- 1 + index  # 1, because raw scores begin at 0
    results[3, "Centil"] <- trans.table[index]
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
                            sep=""),
                      paste("\nNormative data from a non-clinical British sample. Data from", ### CODE INSERTED MANUALLY
                            "JR Crawford & JD Henry, Br J Clin Psychol 2003, 42:111-131.",    ###
                            sep='\n')                                                         ###
                     )

  # Return results
  # ------------------
  list(results.lst = results.lst,
       results.df = results,
       results.scores = results.scores)

} # end of scoring.fun