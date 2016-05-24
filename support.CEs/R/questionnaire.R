questionnaire <- function(choice.experiment.design, common = NULL, quote = TRUE)
{
# Name: questionnaire
# Title: Converting a choice experiment design into a choice experiment questionnaire
# Arguments:
#  choice.experiment.design   A data frame containing a choice experiment design created
#                               by the function Lma.design() or rotation.design().
#  common                     A vector containing a fixed combination of attribute-levels
#                               corresponding to a common base option in each question. 
#                               If there is no common base option, the argument is set as NULL.          
#  quote                      A logical variable indicating whether or not the attribute-levels
#                               in each question are printed wth quotation marks.



# set variables and vector

  nblocks <- choice.experiment.design$design.information$nblocks
  nquestions <- choice.experiment.design$design.information$nquestions
  nalternatives <- choice.experiment.design$design.information$nalternatives
  nattributes <- choice.experiment.design$design.information$nattributes
  attribute.names <- names(choice.experiment.design[[1]][[1]])[-(1:3)]

# integrate alternatives into data.frame

  my.design <- as.matrix(choice.experiment.design[[1]][[1]])

  if (nalternatives >= 2) {
    for (i in 2:nalternatives) {
      my.design <- rbind(my.design,
                         as.matrix(choice.experiment.design$alternatives[[i]]))
    }
  }

  if (is.null(common) == FALSE) {
    nalternatives <- nalternatives + 1
    common.base <- choice.experiment.design$alternatives[[1]]
    common.base[, 3] <- nalternatives
    common.base <- as.matrix(common.base)
    for (i in attribute.names) {
      common.base[, i] <- common[[i]]
    }
    my.design <- rbind(my.design, common.base)
  }

  rownames(my.design) <- NULL 

  my.design <- data.frame(my.design)
  my.design$BLOCK <- as.numeric(as.character(my.design$BLOCK))
  my.design$QES <- as.numeric(as.character(my.design$QES))
  my.design$ALT <- as.numeric(as.character(my.design$ALT))

# convert the choice experiment design into questions

  my.design <- my.design[order(my.design$BLOCK, my.design$QES, my.design$ALT),]
  alternative.names <- paste("alt.", 1:nalternatives, sep= "")
  cat("\n")
  for (i in 1:nblocks) {
    cat("Block", i, "\n", "\n")
    for (j in 1:nquestions) {
      cat("Question", j, "\n")
      temp <- subset(my.design, my.design$BLOCK == i & my.design$QES == j)
      temp <- temp[, 4:(3 + nattributes)]
      if (nattributes == 1) {
        temp <- as.data.frame(temp)
        names(temp) <- attribute.names
      }
      temp <- t(temp)
      colnames(temp) <- alternative.names
      print(temp, quote = quote)
      cat("\n")
    }
  }
}

