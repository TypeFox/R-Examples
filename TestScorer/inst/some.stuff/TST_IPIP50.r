# IPIP50 scale scoring script
# Creation date: 2016-02-18
# --------------

testChar <- list(acronym = "IPIP50",
                 name = "International Personality Item Pool - 50 items",
                 ref = "Goldberg et al. J Res Pers, 2006:40, 84-96.",
                 n.items = 50,
                 valid = c(1, 2, 3, 4, 5),
                 miss = c(0),
                 comm = "Public domain: www.ipip.ori.org/New_IPIP-50-item-scale.htm")

scoring.fun <- function(answers, sex, age, id, date.test, comm) {
  # "answer" is a *character* vector as introduced through the keyboard.
  # "sex", "age", "id", "date.test" & "comm" can be used if appropiate.
  if(!(age %in% 18:86)) {                                    ### Manually inserted code
    library(tcltk)                                           ###
    msg <- 'Age should be a number between 18 and 86.'       ###
	tkmessageBox(message=msg, icon = 'warning')              ###
	stop(msg)                                                ###
	}                                                        ###
  answers <- as.numeric(answers) # transform to numeric for easier scoring   
  answers[answers %in% c(0)] <- NA # missing characters to NA
  blanks <- sum(is.na(answers)) # compute number of missing items
  pcnt.blanks <- round((blanks / 50) * 100) # compute % of missings
  # Items which should be reversed
  reversed.items=c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 29, 30, 32, 34, 36, 38, 39, 44, 46, 49)
  answers[reversed.items] <- (5 + 1) - answers[reversed.items]
  results <- data.frame(NULL) # Initialize null data frame for results

  toT <- function(raw.score, mean, sd) {  # compute T score
    T.score <- round(((raw.score - mean) / sd) * 10 + 50)
    T.score
  } # end toT

  makeGraph <- function(T.score) {  # make a graph
    template  <- "|    :    :    |    :    |    :    |    :    :    |"
    options(warn=-1)
    T.score <- as.integer(T.score)  # NA in case of non-numerical string
    options(warn=0)
    if (!is.na(T.score)) {
      if (T.score < 0) T.score <- 0
        else if (T.score > 100) T.score <- 100
      position <- round((T.score / 2) + 1)
      graph <- paste(substr(template, 1, position-1),
                     substr(template, position + 1, nchar(template)),
                     sep="o")  # "o" marks the position
    } else {
      graph <- "Not graphicable"
    }
    graph
  } # end makeGraph

  # E scale scoring commands
  # --------------
  results[1, "Acronym"] <- "E" # acronym
  results[1, "Scale"] <- "Extraversion" # name of scale
  # Items making up the scale
  items <- c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46)
  results[1, "Miss"] <- sum(is.na(answers[items])) # number of missings
  if (results[1, "Miss"] == length(items)) {  # all answers are missings
    results[1, "Raw"] <- NA
    results[1, "T"] <- NA
    results[1, "Graph"] <- "All items missing"
  } else {
    results[1, "Raw"] <- sum(answers[items], na.rm=TRUE) # sum answered items
	# compute T score                                                                    ###  Manually inserted code
    if (sex=="Male") {                                                                   ###
	  if (age >= 18 & age < 30) results[1, "T"] <- toT(results[1, "Raw"], 25.4, 7.6)     ###
	  if (age >= 30 & age < 65) results[1, "T"] <- toT(results[1, "Raw"], 21.7, 8)       ###
	  if (age >= 65 & age <= 86) results[1, "T"] <- toT(results[1, "Raw"], 20.7, 7.7)    ###
	} else {                                                                             ###
	  if (age >= 18 & age < 30) results[1, "T"] <- toT(results[1, "Raw"], 24.3, 7.1)     ###
	  if (age >= 30 & age < 65) results[1, "T"] <- toT(results[1, "Raw"], 21.9, 8.5)     ###
	  if (age >= 65 & age <= 86) results[1, "T"] <- toT(results[1, "Raw"], 21.6, 7.3)    ###
	}                                                                                    ###
    results[1, "Graph"] <- makeGraph(results[1, "T"]) # makes the graph
  }

  # A scale scoring commands
  # --------------
  results[2, "Acronym"] <- "A" # acronym
  results[2, "Scale"] <- "Agreeableness" # name of scale
  # Items making up the scale
  items <- c(2, 7, 12, 17, 22, 27, 32, 37, 42, 47)
  results[2, "Miss"] <- sum(is.na(answers[items])) # number of missings
  if (results[2, "Miss"] == length(items)) {  # all answers are missings
    results[2, "Raw"] <- NA
    results[2, "T"] <- NA
    results[2, "Graph"] <- "All items missing"
  } else {
    results[2, "Raw"] <- sum(answers[items], na.rm=TRUE) # sum answered items
	# compute T score                                                                         ### Manually inserted code
    if (sex=="Male") {                                                                        ###
	  if (age >= 18 & age < 30) results[2, "T"] <- toT(results[2, "Raw"], 28.3, 4.7)          ###
	  else if (age >= 30 & age < 65) results[2, "T"] <- toT(results[2, "Raw"], 29.5, 6.2)     ###
	  else if (age >= 65 & age <= 86) results[2, "T"] <- toT(results[2, "Raw"], 29.8, 5.2)    ###
	} else {                                                                                  ###
	  if (age >= 18 & age < 30) results[2, "T"] <- toT(results[2, "Raw"], 30.9, 4.6)          ###
	  else if (age >= 30 & age < 65) results[2, "T"] <- toT(results[2, "Raw"], 32.5, 6.1)     ###
	  else if (age >= 65 & age <= 86) results[2, "T"] <- toT(results[2, "Raw"], 33.2, 4.8)    ###
	}
    results[2, "Graph"] <- makeGraph(results[2, "T"]) # makes the graph
  }

  # C scale scoring commands
  # --------------
  results[3, "Acronym"] <- "C" # acronym
  results[3, "Scale"] <- "Conscientiousness" # name of scale
  # Items making up the scale
  items <- c(3, 8, 13, 18, 23, 28, 33, 38, 43, 48)
  results[3, "Miss"] <- sum(is.na(answers[items])) # number of missings
  if (results[3, "Miss"] == length(items)) {  # all answers are missings
    results[3, "Raw"] <- NA
    results[3, "T"] <- NA
    results[3, "Graph"] <- "All items missing"
  } else {
    results[3, "Raw"] <- sum(answers[items], na.rm=TRUE) # sum answered items
 	# compute T score                                                                        ###  Manually inserted code
   if (sex=="Male") {                                                                        ###
	  if (age >= 18 & age < 30) results[3, "T"] <- toT(results[3, "Raw"], 21.6, 6.7)         ###
	  else if (age >= 30 & age < 65) results[3, "T"] <- toT(results[3, "Raw"], 26.7, 6.1)    ###
	  else if (age >= 65 & age <= 86) results[3, "T"] <- toT(results[3, "Raw"], 28.5, 6.2)   ###
	} else {                                                                                 ###
	  if (age >= 18 & age < 30) results[3, "T"] <- toT(results[3, "Raw"], 22.6, 6.6)         ###
	  else if (age >= 30 & age < 65) results[3, "T"] <- toT(results[3, "Raw"], 26.6, 6)      ###
	  else if (age >= 65 & age <= 86) results[3, "T"] <- toT(results[3, "Raw"], 28.4, 6.2)   ###
	}                                                                                        ###
    results[3, "Graph"] <- makeGraph(results[3, "T"]) # makes the graph
  }

  # ES scale scoring commands
  # --------------
  results[4, "Acronym"] <- "ES" # acronym
  results[4, "Scale"] <- "Emotional stability" # name of scale
  # Items making up the scale
  items <- c(4, 9, 14, 19, 24, 29, 34, 39, 44, 49)
  results[4, "Miss"] <- sum(is.na(answers[items])) # number of missings
  if (results[4, "Miss"] == length(items)) {  # all answers are missings
    results[4, "Raw"] <- NA
    results[4, "T"] <- NA
    results[4, "Graph"] <- "All items missing"
  } else { results[4, "Raw"] <- sum(answers[items], na.rm=TRUE) # sum answered items
	# compute T score 
    if (sex=="Male") {                                                                      ### Manually inserted code
	  if (age >= 18 & age < 30) results[4, "T"] <- toT(results[4, "Raw"], 22, 6.7)          ###
	  else if (age >= 30 & age < 65) results[4, "T"] <- toT(results[4, "Raw"], 21.0, 8.5)   ###
	  else if (age >= 65 & age <= 86) results[4, "T"] <- toT(results[4, "Raw"], 24.7, 8.5)  ###
	} else {                                                                                ###
	  if (age >= 18 & age < 30) results[4, "T"] <- toT(results[4, "Raw"], 18.6, 7.3)        ###
	  else if (age >= 30 & age < 65) results[4, "T"] <- toT(results[4, "Raw"], 20.2, 8.3)   ###
	  else if (age >= 65 & age <= 86) results[4, "T"] <- toT(results[4, "Raw"], 24, 7.9)    ###
	}                                                                                       ###
    results[4, "Graph"] <- makeGraph(results[4, "T"]) # makes the graph
  }

  # I scale scoring commands
  # --------------
  results[5, "Acronym"] <- "I" # acronym
  results[5, "Scale"] <- "Intellect/Imagination" # name of scale
  # Items making up the scale
  items <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
  results[5, "Miss"] <- sum(is.na(answers[items])) # number of missings
  if (results[5, "Miss"] == length(items)) {  # all answers are missings
    results[5, "Raw"] <- NA
    results[5, "T"] <- NA
    results[5, "Graph"] <- "All items missing"
  } else {
    results[5, "Raw"] <- sum(answers[items], na.rm=TRUE) # sum answered items
 	# compute T score                                                                      ### Manually inserted code
   if (sex=="Male") {                                                                      ###
	  if (age >= 18 & age < 30) results[5, "T"] <- toT(results[5, "Raw"], 27.9, 5.7)       ###
	  else if (age >= 30 & age < 65) results[5, "T"] <- toT(results[5, "Raw"], 28.5, 6.1)  ###
	  else if (age >= 65 & age <= 86) results[5, "T"] <- toT(results[5, "Raw"], 24, 6.1)   ###
	} else {                                                                               ###
	  if (age >= 18 & age < 30) results[5, "T"] <- toT(results[5, "Raw"], 26.8, 5.1)       ###
	  else if (age >= 30 & age < 65) results[5, "T"] <- toT(results[5, "Raw"], 26.6, 6)    ###
	  else if (age >= 65 & age <= 86) results[5, "T"] <- toT(results[5, "Raw"], 23.9, 6)   ###
	}                                                                                      ###
    results[5, "Graph"] <- makeGraph(results[5, "T"]) # makes the graph
  }

  # Vector for writing scores to a file
  # --------------------
  results.scores <- unlist(results[-c(1, 2, 6)])  # not Acronym, Scale & Graph columns
  names <- paste(results$Acronym, names(results.scores), sep=".")
  names(results.scores) <- sub("[0-9]+$", "", names)  # delete ending numbers

  # Ruler for graph column
  # --------------------
  names(results)[6] <- "0    10   20   30   40   50   60   70   80   90  100"

  # Output in form of list
  # ------------------
  results.lst <- list(paste("Total number of missing items: ",
                            blanks,
                            " (",
                            pcnt.blanks,
                            "%)",
                            sep=""),
                      "",                                                                                    ### Manually inserted code
                      'Norms from a mixed sample of Scottish students and general population volunteers.',   ###
                      'Reported by Gow AJ et al. Pers Indiv Differ 2005, 39:317-329.',                       ###
					  'Age groups: 18 to 29, 30 to 64 and 65 to 86.'                                         ###
                     )                            
							

  # Return results
  # ------------------
  list(results.lst = results.lst,
       results.df = results,
       results.scores = results.scores)

} # end of scoring.fun