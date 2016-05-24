# RAAS scale scoring script
# Creation date: 2014-01-07
# --------------

testChar <- list(acronym = "RAAS",
                 name = "Revised Adult Attachment Scale",
                 ref = "Collins N, 1996",
                 n.items = 18,
                 valid = c(1, 2, 3, 4, 5),
                 miss = c(0),
                 comm = "Public domain: www.openpsychassessment.org/wp-content/uploads/2011/06/AdultAttachmentScale.pdf")

scoring.fun <- function(answers, sex, age, id, date.test, comm) {
  # "answer" is a *character* vector as introduced through the keyboard.
  # "sex", "age", "id", "date.test" & "comm" used if appropiate.
  answers <- as.numeric(answers) # transform to numeric for easier scoring
  answers[answers %in% c(0)] <- NA # missing characters to NA
  blanks <- sum(is.na(answers)) # compute number of missing items
  pcnt.blanks <- round((blanks / 18) * 100) # compute % of missing items
  # Items which should be reversed
  reversed.items=c(2, 7, 8, 13, 16, 17, 18)
  answers[reversed.items] <- (5 + 1) - answers[reversed.items]
  results <- data.frame(NULL) # Initialize null data frame for results

  toT <- function(raw.score, mean, sd) {  # compute T score
    T.score <- round(((raw.score - mean) / sd) * 10 + 50)
    T.score
  } # end toT

  makeGraph <- function(T.score) {  # make a graph
    template  <- "|    :    :    |    :    |    :    |    :    :    |"
    options(warn=-1)
    T.score <- as.integer(T.score)  # in case of numerical string
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

  # C scale scoring commands
  # --------------
  results[1, "Acronym"] <- "C" # acronym
  results[1, "Scale"] <- "Close" # name of the scale
  # Items making up the scale
  items <- c(1, 6, 8, 12, 13, 17)
  results[1, "Miss"] <- sum(is.na(answers[items])) # number of missing items
  if (results[1, "Miss"] == length(items)) {  # all answers are missing items
    results[1, "Raw"] <- NA
    results[1, "T"] <- NA
    results[1, "Graph"] <- "All missings"
  } else {
    results[1, "Raw"] <- round(mean(answers[items], na.rm=TRUE), 2) # mean answered items
    if (sex=="Male") results[1, "T"] <- toT(results[1, "Raw"], 3.59, 0.87) # compute T score
      else results[1, "T"] <- toT(results[1, "Raw"], 3.65, 0.87)
    results[1, "Graph"] <- makeGraph(results[1, "T"]) # make the graph
  }

  # D scale scoring commands
  # --------------
  results[2, "Acronym"] <- "D" # acronym
  results[2, "Scale"] <- "Dependent" # name of the scale
  # Items making up the scale
  items <- c(2, 5, 7, 14, 16, 18)
  results[2, "Miss"] <- sum(is.na(answers[items])) # number of missing items
  if (results[2, "Miss"] == length(items)) {  # all answers are missing items
    results[2, "Raw"] <- NA
    results[2, "T"] <- NA
    results[2, "Graph"] <- "All missings"
  } else {
    results[2, "Raw"] <- round(mean(answers[items], na.rm=TRUE), 2) # mean answered items
    if (sex=="Male") results[2, "T"] <- toT(results[2, "Raw"], 3.43, 0.83) # compute T score
      else results[2, "T"] <- toT(results[2, "Raw"], 3.25, 0.86)
    results[2, "Graph"] <- makeGraph(results[2, "T"]) # make the graph
  }

  # A scale scoring commands
  # --------------
  results[3, "Acronym"] <- "A" # acronym
  results[3, "Scale"] <- "Anxiety" # name of the scale
  # Items making up the scale
  items <- c(3, 4, 9, 10, 11, 15)
  results[3, "Miss"] <- sum(is.na(answers[items])) # number of missing items
  if (results[3, "Miss"] == length(items)) {  # all answers are missing items
    results[3, "Raw"] <- NA
    results[3, "T"] <- NA
    results[3, "Graph"] <- "All missings"
  } else {
    results[3, "Raw"] <- round(mean(answers[items], na.rm=TRUE), 2) # mean answered items
    if (sex=="Male") results[3, "T"] <- toT(results[3, "Raw"], 2.31, 0.89) # compute T score
      else results[3, "T"] <- toT(results[3, "Raw"], 2.62, 1.01)
    results[3, "Graph"] <- makeGraph(results[3, "T"]) # make the graph
  }

  # Vector for writing scores to a file
  # --------------------
  results.scores <- unlist(results[-c(1, 2, 6)])  # not Acronym, Scale & Graph columns
  names <- paste(results$Acronym, names(results.scores), sep=".")
  names(results.scores) <- sub("[0-9]+$", "", names)  # delete ending numbers

  # Ruler for graph column
  # --------------------
  names(results)[6] <- "0    10   20   30   40   50   60   70   80   90  100"

  # Combine C & D scales                                                                 ### CODE INSERTED MANUALLY
  CD <- round(mean(c(results[1, 'Raw'], results[2, 'Raw'])), 2)                          ###
  results[4, ] <- c("", "", "", "", "", "") # blank row to improve readability           ###
  results[5, ] <- c("CD", "Close/Dependent", "", CD, "", "") # no data for T score       ###

  # Attachment style assignement                                                         ###
  if (is.na(CD)) style <- "Not evaluable" # if all answers are missings                  ###
  else if (CD > 3 & results[3, "Raw"] < 3) style <- "Secure"                             ###
  else if (CD > 3 & results[3, "Raw"] > 3) style <- "Preoccupied"                        ###
  else if (CD < 3 & results[3, "Raw"] < 3) style <- "Dismissing"                         ###
  else if (CD < 3 & results[3, "Raw"] > 3) style <- "Fearful"                            ###
  else style <- "Not classificable"                                                      ###

  # Show style as a plot                                                                 ### CODE INSERTED MANUALLY
  windows(title="Attachment style")                                                      ###
  require(graph)                                                                         ###
  plot(CD, results[3, "Raw"], xlim=c(1,5), ylim=c(1,5), pch=3, cex=2, col="blue", lwd=5, ###
       xlab="Close/Dependent", ylab="Anxiety", main=paste(id, date.test),                ###
       font.sub=2, sub="The position of the subject is represented by a blue cross")     ###
  abline(v=3)                                                                            ###
  abline(h=3)                                                                            ###
  text(2, 2, labels="Dismissing", col="gray60", font=2, cex=2)                           ###
  text(4, 2, labels="Secure", col="gray60", font=2, cex=2)                               ###
  text(2, 4, labels="Fearful", col="gray60", font=2, cex=2)                              ###
  text(4, 4, labels="Preocupied", col="gray60", font=2, cex=2)                           ###
  
  # Output in form of list
  # ------------------
  results.lst <- list(paste("Total number of missing items: ",
                            blanks,
                            " (",
                            pcnt.blanks,
                            "%)",
                            sep=""),
                      "",                                                                                    ### CODE INSERTED MANUALLY
                      'According to the author, attach styles assignement "is quite exploratory...',         ###
                      '[use] with caution, and only in conjunction with the continuous measures."',          ###
                      "",                                                                                    ###
                      paste("Attach style:", style),                                                         ###
                      "",                                                                                    ###
                      "T-scores computed using mean and standard deviation from 414 USA college students,",  ###
                      "reported by Ledley et al. J Psychopath Behav Assess 2006, 28:33-40."                  ###
                     )                            

  # Return results
  # ------------------
  list(results.lst = results.lst,
       results.df = results,
       results.scores = results.scores)

} # end of scoring.fun