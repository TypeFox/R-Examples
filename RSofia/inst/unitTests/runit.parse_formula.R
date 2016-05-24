## ------------------------------------------------------------------------- ##
## Tests for setConfiguration                                                ##
##                                                                           ##
## ------------------------------------------------------------------------- ##

cat("\n\nRUnit test cases for 'RSofia:::parse_formula' function\n\n")

# This function is executed at the end of the unit test
.tearDown <- function() {
  if (exists("track", envir=.GlobalEnv)) {
    rm(track, envir=.GlobalEnv)
  }
}

data(irismod)

### needs to cover factors
training_data <- data.frame(
  Y  = c(-1,  1,  -1, 1   ),
  X1 = c(.1, .2,  .3, .4  ),
  X2 = c(1 , 10, 100, 1000),
  X3 = c(4 ,  3,   2, 1   )
) 


### load training data

### should fail
test.parse_formula.bad_args <- function() {

  checkException(parse_formula("not a formula", training_data  ))
  checkException(parse_formula(Y ~ ., "not a dataset"))

}

test.parse_formula.result_type <- function() {
  
  irismod.parsed <- parse_formula(Is.Virginica ~ ., irismod )  

  checkEquals(names(irismod.parsed), c("data", "labels", "no_bias_term"))
  checkEquals(class(irismod.parsed), "list")
  
  ### check individual types

  checkTrue(is.matrix(irismod.parsed$data) && is.numeric(irismod.parsed$data))
  checkTrue(is.numeric(irismod.parsed$labels))
  ##checkTrue(!is.null(attr(irismod.parsed$labels, "name")))
  checkTrue(is.logical(irismod.parsed$no_bias_term))

}

test.parse_formula.result_values <- function() {


  ### with Y ~ .
  
  training_data.parsed <- parse_formula(Y ~ ., training_data )
  
  checkEqualsNumeric(training_data.parsed$labels, training_data$Y)
  checkEqualsNumeric(training_data.parsed$data  , as.matrix(training_data[,c("X1","X2","X3")] ))
  checkEquals(training_data.parsed$no_bias_term, FALSE)

  rm(training_data.parsed)
  
  ### with Y ~ X1 + X2a

  training_data.parsed <- parse_formula(Y ~ X1 + X2, training_data )
  
  checkEqualsNumeric(training_data.parsed$labels, training_data$Y)
  checkEqualsNumeric(training_data.parsed$data  , as.matrix(training_data[,c("X1","X2")]))
  checkEquals(training_data.parsed$no_bias_term, FALSE)

  rm(training_data.parsed)

  ### with no intercept

  training_data.parsed <- parse_formula(Y ~ -1 + X1 + X2, training_data )
  
  checkEqualsNumeric(training_data.parsed$labels, training_data$Y)
  checkEqualsNumeric(training_data.parsed$data  , as.matrix(training_data[,c("X1","X2")]))
  checkEquals(training_data.parsed$no_bias_term, TRUE)

  rm(training_data.parsed)
  ###

}
