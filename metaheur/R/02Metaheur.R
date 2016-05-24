
#' metaheur
#'
#' Metaheuristic optimization of preprocessing combinations
#'
#' @param gridclassobject (GridClass) created by setgrid function in preprocomb package
#' @param startgrid (integer) 0 random restart (default), 1 grid restart
#' @param startnum (integer) number of restarts
#' @param iterations (integer) number of iterations done for a restart, default to 5
#' @param taboolistlength (integer) number of previous solution that can not be revisited, > 1
#' @param initialtemperature (numeric) initial propability for acccepting an inferior candidate, between 0 and 1
#' @param tempconst (numeric) multiplier for decreasing temperature on each iteration
#' @param reheat (numeric) propability of increasing temperature on each iteration
#' @param nholdout (integer) number of holdout rounds
#' @param late (integer) location of previous best solution a candidate is compared to, default to 0 for last
#' @param stopcond (integer) type of stopping condition, default to 1 for threshold, 2 for convergence
#' @param stopvalue (numeric) threshold for stopping, defaults to 0.95
#' @param deltafive (numeric) convergence criteria for last five iteration, defaults to 0.05
#' @param verbose (boolean)
#' @param predictor (character) caret name of predictive model, defaults to "knn"
#' @examples ## result <- metaheur(examplegrid,verbose=TRUE)
#' ## getbestheur(result)
#' @export

metaheur <- function(gridclassobject, startgrid=0, startnum=1, iterations=5, taboolistlength=1,
                                    initialtemperature=0.01, tempconst=0.01, reheat=0.01, nholdout=1, late=0,
                                    stopcond=1, stopvalue=0.95, deltafive=0.05, verbose=FALSE, predictor="knn"){

# TECHNICAL INITIALIZATION

  if (verbose==TRUE){
  if (startgrid==0){print("Start type: random restarts.")}
  if (startgrid==1){print("Start type: grid restarts.")}
  }

  if (verbose==TRUE) {cat("Number of restarts:", startnum, "\n")}

  # initializations
  predictors <- predictor
  grid <- gridclassobject
  reslist <- list()
  reshistory <- list()
  r <- numeric()
  s2 <- numeric()

  # unique preprocesors in a phase
  uniquepreprocessors <- lapply(grid@grid, getuniquepreprocessors)

  if (startgrid==0) {gridsequence <- sample(1:nrow(grid@grid),startnum)} else
    {gridsequence <- seq(1, nrow(grid@grid), by=ceiling(nrow(grid@grid)/startnum))}



# RESTART STRUCTURE

counterforrestarts <- 1

for (k  in gridsequence)
{
  taboo <- list()
  history <- numeric()
  temperature <- initialtemperature

  numberofcurrentbest <- sample(1:nrow(grid@grid),1) # random start
  if (verbose==TRUE) {cat("Start combination:", numberofcurrentbest, "\n")}

  currentbest_value <- getfirstassessment(numberofcurrentbest, grid, predictors, nholdout)

# ITERATION STRUCTURE WITHIN A RESTART

  for (j in 1:iterations)

  {

    # INITIALIZATION PROCEDURE

    copyofcurrentbest <- grid@grid[numberofcurrentbest,]
    if (verbose==TRUE) {cat("Iteration:", j, "Current best:", unlist(copyofcurrentbest), currentbest_value, "\n")}


    # setting taboo and history for first iteration
    if (j==1) {
      taboo <- append(taboo, list(copyofcurrentbest))
      history <- c(history, currentbest_value)
    }

    # MODIFICATION PROCEDURE

    candidate_new <- getnewcandidate(grid, taboo, taboolistlength, uniquepreprocessors, copyofcurrentbest)

    taboo <- append(taboo, list(candidate_new))

    # ASSESSMENT PROCEDURE

    candidatereturn <- getcandidatedata(grid, candidate_new)

    candidatedata <- candidatereturn[[1]]
    candidatenumber <- candidatereturn[[2]]

    candidate_value <- getconsequentassessment(candidatedata, predictors, nholdout)

    if (verbose==TRUE) {cat("Iteration:", j, "Candidate:", unlist(candidate_new), candidate_value, "\n")}

    # SELECTION PROCEDURE

    # Decreasing temperature
    temperature <- temperature*tempconst
    if (temperature < 0.00001) {temperature <- 0.00001}

    # Reheating
    reheatprobability <- sample(0:1, 1, prob=c(1-reheat, reheat))
    if (reheatprobability==1) {
      if (verbose==TRUE) {cat("Reheating occured.", "\n")}
      temperatureincrease <- (1-temperature)*stats::runif(1,0,1)
      temperature <- temperature + temperatureincrease}

    if (verbose==TRUE){cat("Temperature:", temperature, "\n")}

    # Late acceptance Hill-Climbing

    historycondition <- (length(history)-late) <= 0
    if (historycondition==TRUE){historycomparison <- history[1]} else { # first iteration
      historycomparison <- history[length(history)-late] # consequent iterations
    }
    if (verbose==TRUE) {cat("Comparison value for late acceptance:", historycomparison, "\n")}

    if (candidate_value > historycomparison) {
      currentbest_value <- candidate_value
      numberofcurrentbest <- candidatenumber}

    # Simulated annealing
    else {
      acceptinginferiorcandidate <- sample(0:1, 1, prob=c(1-temperature, temperature))
      if (acceptinginferiorcandidate==1)
      {currentbest_value <- candidate_value
      numberofcurrentbest <- candidatenumber
       if (verbose==TRUE) {cat("SA: A weaker solution was accepted.","\n")}
        }
    }

    # Storing of

    history <- c(history, currentbest_value)
    historydelta <- mean(utils::tail(diff(history),5))
    if (verbose==TRUE) {cat("History delta, last five:", historydelta, "\n")}

    ## TERMINATION PROCEDURE

    if (stopcond==1){
      if (currentbest_value > stopvalue)
        {cat("Iteration stop condition reseached.", "\n")
        break}
    }
    if (stopcond==2){
      if (historydelta < deltafive ) {cat("Iteration convergence condition reseached.", "\n")
        break}
    }


  } # ITERATION STRUCTURE EXIT

  reslist[[counterforrestarts]] <- list(combination=grid@grid[numberofcurrentbest,], accuracy=currentbest_value)
  reshistory[[counterforrestarts]] <- history

  counterforrestarts <- counterforrestarts + 1

} # RESTART STRUCTURE EXIT

# Preraration of function return object, TODO: clean this

best <- which.max(unlist(lapply(reslist, function(x) x$accuracy)))
best <- data.frame(data.frame(unlist(reslist[[best]]$combination)), accuracy=reslist[[best]]$accuracy)
list(best, reshistory, reslist)

} # FUNCTION EXIT

#' getbestheur
#'
#' get the best combination and its classification accuracy
#'
#' @param x output of function metaheur()
#' @return (list) best combination, classification accuracy of the best combination
#' @examples result <- metaheur(examplegrid)
#' getbestheur(result)
#' @export

getbestheur <- function(x){
  q <- data.frame(t(x[[1]])[1,])
  colnames(q) <- "preprocessor"
  r <- unname(unlist(x[[2]]))[1]
  list(q,r)
}





