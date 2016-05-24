bws.questionnaire <- function(choice.sets, design.type, item.names)
{
# Name : bws.questionnaire
# Title: Converting a two-level OMED/BIBD into BWS questions
# Arguments:
#  choice.sets   a data frame or matrix containing choice sets
#  design.type   a value describing how to design the choice sets
#                  = 1 if a OMED is used
#                  = 2 if a BIBD is used
#  item.names    a vector containing the names of items shown in the questions


# set variable and matrix

  numQuestions <- nrow(choice.sets)
  design <- data.matrix(choice.sets)


# store items contained in each choice set in each element of itemsInSet 

  itemsInSet <- vector("list", numQuestions) # list
  if(design.type == 1)   # OMED
  {
    for(i in 1:numQuestions)
    {
      itemsInSet[[i]] <- which(design[i, ] == 2)
    }
  }
  else                   # BIBD
  {
    for(i in 1:numQuestions)
    {
      itemsInSet[[i]] <- design[i, ]
    }
  }


# display BWS questions

  numItemsInSet <- sapply(itemsInSet, length)  # number of items in each choice set

  for(i in 1:numQuestions) {
    cat("\n")
    cat("Q", i, "\n", sep ="")
    dsp           <- matrix(c("[ ]"), nrow = numItemsInSet[i], ncol = 3)
    dsp[, 2]      <- item.names[itemsInSet[[i]]]
    colnames(dsp) <- c("Best", "Items", "Worst")
    rownames(dsp) <- rep(c(""), numItemsInSet[i])
    print(noquote(dsp))
  }
  cat("\n")

}

