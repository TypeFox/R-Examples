## ------------------------------------------------------------------------- ##
## Tests for setConfiguration                                                ##
##                                                                           ##
## ------------------------------------------------------------------------- ##

cat("\n\nRUnit test cases for 'RSofia:::sofia' function\n\n")

# This function is executed at the end of the unit test
.tearDown <- function() {
  if (exists("track", envir=.GlobalEnv)) {
    rm(track, envir=.GlobalEnv)
  }
}

# Now the unit tests, following the naming convention test.XXX

# convention test.sofia.learner_type.loop_type

### load result data

data(irismod)
data(sofia_ml_train)

RANDOM_SEED <- 1L
TOLERANCE <- .00001

LEARNER_TYPE <- c( "pegasos", "sgd-svm", "passive-aggressive", "margin-perceptron", "romma", "logreg-pegasos" )
LOOP_TYPE    <- c( "stochastic", "balanced-stochastic", "rank", "roc", "query-norm-rank", "combined-ranking", "combined-roc" )
ETA_TYPE     <- c( "pegasos", "basic", "constant" )
NO_BIAS_TERM <- c( TRUE, FALSE)

create_function_text <- function(learner, loop, eta, no_bias_term) {

  function_txt.name <- paste("test.sofia."
    , gsub("-","_",learner),"."
    , gsub("-","_",loop), "."
    , gsub("-","_",eta), "."
    , paste("nobias", as.integer(no_bias_term),sep=""), " <- function()"
    , " { \n", sep ="")

  function_txt.body <- paste(
      
      "  weights.RSofia <- sofia(Is.Virginica ~", ifelse(no_bias_term, "-1 +", ""),
      " ., irismod, random_seed = RANDOM_SEED", "\n",
      "  ,eta_type = '", eta,
      "'  ,learner_type ='", learner,
      "'  ,loop_type = '", loop,
      "'", 
      "  )$weights", "\n\n",
      "  weights.sofia_ml <- sofia_ml_train[[\'",
         paste(learner, loop, eta, paste("nobias", as.integer(no_bias_term),sep=""), sep = '_'),
         "\']]","\n\n",
      "  checkEqualsNumeric(weights.RSofia, weights.sofia_ml, tolerance = TOLERANCE)", sep = "")
                      
  function_txt.end <- "\n}"
                        
  function_txt <- paste(function_txt.name, function_txt.body, function_txt.end, sep = "\n")
                        
  return(function_txt)

}

### this can be expanded into code obviously

for(i in LEARNER_TYPE)
  for(j in LOOP_TYPE)
    for(k in ETA_TYPE)
      for(l in NO_BIAS_TERM)
        eval(parse(text = create_function_text(i,j,k,l)))() 

options(warn = 2)

test.sofia.bad_parameters <- function() {

  checkException(sofia("xxx",irismod))
  checkException(sofia(Is.Virginica ~ ., "hello"))
  checkException(sofia(Is.Virginica ~ ., irismod, random_seed = "a"))
  checkException(sofia(Is.Virginica ~ ., irismod, random_seed = c(1,2,3)))
  checkException(sofia(Is.Virginica ~ ., irismod, lambda = "a")) 
  checkException(sofia(Is.Virginica ~ ., irismod, lambda = c(1,2,3))) 
  checkException(sofia(Is.Virginica ~ ., irismod, learner_type = "a"))
  checkException(sofia(Is.Virginica ~ ., irismod, eta_type = "a"))
  checkException(sofia(Is.Virginica ~ ., irismod, loop_type = "a"))
  checkException(sofia(Is.Virginica ~ ., irismod, rank_step_probability = "a"))
  checkException(sofia(Is.Virginica ~ ., irismod, rank_step_probability = c(1,2,3)))
  checkException(sofia(Is.Virginica ~ ., irismod, rank_step_probability = -1))
  checkException(sofia(Is.Virginica ~ ., irismod, passive_aggressive_lambda = "a")) 
  checkException(sofia(Is.Virginica ~ ., irismod, passive_aggressive_lambda = c(1,2,3)))
  checkException(sofia(Is.Virginica ~ ., irismod, passive_aggressive_c = "a"))
  checkException(sofia(Is.Virginica ~ ., irismod, passive_aggressive_c = c(1,2,3))) 
  checkException(sofia(Is.Virginica ~ ., irismod, perceptron_margin_size = "a")) 
  checkException(sofia(Is.Virginica ~ ., irismod, perceptron_margin_size = c(1,2,3))) 
  checkException(sofia(Is.Virginica ~ ., irismod, training_objective = "a"))
  checkException(sofia(Is.Virginica ~ ., irismod, training_objective = c(TRUE, FALSE))) 
  checkException(sofia(Is.Virginica ~ ., irismod, hash_mask_bits = "a")) 
  checkException(sofia(Is.Virginica ~ ., irismod, hash_mask_bits = c(1,2,3))) 

}

### sofia.character


test.sofia.character <- function() {

  data(irismod)
  x <- parse_formula(Is.Virginica ~ ., irismod)
  tmp <- tempfile()
  write.svmlight(x$labels, x$data, tmp)
  checkException(sofia(tmp, dimensionality=1))
  unlink(tmp)

}

