make.dataset <- function(respondent.dataset, design.matrix, choice.indicators, detail = FALSE) 
{
# Name: make.dataset
# Title: Making a data set
# Arguments:
#  respondent.dataset   A data frame containing respondents' answers to choice experiment questions.
#  design.matrix        A data frame containing a design matrix created by the function make.design.matrix()
#  choice.indicators    A vector of variables showing the alternative of which was selected in 
#                         each choice experiment question.
#  detail               A logical variable describing whether or not some variables contained in 
#                         the argument respondent.dataset and variables created in this function are
#                         stored in a data set produced by this function.



# initial setting

  nquestions <- length(choice.indicators)
  nalternatives <- length(table(design.matrix$ALT))
  nrespondents <- length(respondent.dataset$ID)

# convert respondent.dataset into "one-row-one-alternative" style data set

  my.respondent.dataset <- rbind(respondent.dataset, respondent.dataset)
  if (nquestions * nalternatives > 2) {
    for (i in 1:(nquestions * nalternatives - 2)) {
      my.respondent.dataset <- rbind(my.respondent.dataset, respondent.dataset)
    }
  }
  my.respondent.dataset <- my.respondent.dataset[order(my.respondent.dataset$ID), ]

# add QES and ALT variables to respondent data set 

  temp.BLOCK <- my.respondent.dataset$BLOCK
  col.BLOCK <- which(colnames(my.respondent.dataset) == "BLOCK")
  my.respondent.dataset <- subset(my.respondent.dataset, select = -col.BLOCK)
  my.respondent.dataset$BLOCK <- temp.BLOCK
  my.respondent.dataset$QES <- rep(1:nquestions,
                                   each = nalternatives,
                                   times = nrespondents)
  my.respondent.dataset$ALT <- rep(1:nalternatives,
                                   times = nrespondents * nquestions)

# convert choice.indicators into SELECT variable

  col.choice.indicators <- rep(0, nquestions)
  for (i in 1:nquestions) {
    col.choice.indicators[i] <- which(colnames(my.respondent.dataset) ==
                                      choice.indicators[i])
  }
  m <- c(1)
  for (i in 1:nrespondents) {
    for (j in col.choice.indicators) {
      for (k in 1:nalternatives) {
        my.respondent.dataset$SELECT[m] <- my.respondent.dataset[m,j]
        m <- m + c(1)
      }
    }
  }

# convert SELECT variable into RES variable

  my.respondent.dataset$RES <- my.respondent.dataset$SELECT == my.respondent.dataset$ALT

# combine respondent data set and design matrix

  my.respondent.dataset$mt <- my.respondent.dataset$ALT +
                              my.respondent.dataset$QES * 100 +
                              my.respondent.dataset$BLOCK * 10000
  design.matrix$mt <- design.matrix$ALT +
                      design.matrix$QES * 100 +
                      design.matrix$BLOCK * 10000
  col.AQB <- colnames(design.matrix) %in% c("ALT", "QES", "BLOCK")
  design.matrix <- subset(design.matrix, select = which(col.AQB == FALSE))
  dataset <- merge(my.respondent.dataset, design.matrix, by = "mt")

# add STR variable

  dataset$STR <- dataset$QES + dataset$ID * 100

# format output

  col.mt <- which(colnames(dataset) == "mt")
  dataset <- subset(dataset, select = -col.mt)
  if (detail == FALSE) {
    exclude.col <- colnames(dataset) %in% c(choice.indicators, "SELECT")
    include.col <- which(exclude.col == FALSE)
    dataset <- subset(dataset, select = include.col)
  }
  dataset <- dataset[order(dataset$STR), ]
  row.names(dataset) <- rep(1:(nquestions * nalternatives * nrespondents))

  return(dataset)
}

