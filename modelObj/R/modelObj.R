#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# modelObj is a configuration object used to define models, methods            #
#   to be used to obtain parameter estimates, and methods to be used to obtain #
#   predictions.                                                               #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# model          : A formula object. Object is model representation.           #
#                                                                              #
# solver.method  : A character giving the R function to be used to obtain      #
#                  parameter estimates. For example, `lm' or `glm'.            #
#                                                                              #
# solver.args    : Additional arguments to be sent to solver.method. This must #
#                  be provided as a list, where the name of each element       #
#                  matches a formal argument of solver.method. For example,    #
#                  if a logistic regression using glm is desired,              #
#                     solver.method = 'glm'                                    #
#                     solver.args = list(family=binomial)                      #
#                                                                              #
#                  It is assumed that solver.method takes formal arguments     #
#                  'formula' and 'data' as input. Occasionally, R methods are  #
#                  developed that do not confirm to this convention.           #
#                  A user can indicate if a different naming convention is     #
#                  used for these two input arguments. For example, if a method#
#                  expects the formula object to be passed through input       #
#                  variable \code{x},                                          #
#                    \code{solver.args} <- list("x"="formula")                 #
#                                                                              #
# predict.method : A function name giving the R function to be used to obtain  #
#                  predicted values. For example, `predict.lm' or              #
#                  `predict.glm'. If not explicitly given, the generic         #
#                  \code{predict} is assumed. Usually, this input does not     #
#                  need to be specified.                                       #
#                                                                              #
# predict.args   : Additional arguments to be sent to predict.method. This     #
#                  must be provided as a list, where the name of each element  #
#                  matches a formal argument of predict.method. For example,   #
#                  if a logistic regression using glm was used to fit the model#
#                  formula object,                                             #
#                     solver.method = 'glm'                                    #
#                     solver.args = list(family=binomial)                      #
#                  then                                                        #
#                     predict.method = 'predict.glm'                           #
#                     predict.args = list(type="response")                     #
#                                                                              #
#                  It is assumed that predict.method takes formal arguments    #
#                  'object' and 'newdata' as input. Occasionally, R methods    #
#                  are developed that do not confirm to this convention.       #
#                  A user can indicate if a different naming convention is     #
#                  used for these two input arguments. For example, if a method#
#                  expects the fit object to be passed through input           #
#                  variable \code{x},                                          #
#                    \code{predict.args} <- list("x"="object")                 #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
buildModelObj <-  function(model, 
                      solver.method=NULL, 
                      solver.args=NULL, 
                      predict.method=NULL, 
                      predict.args=NULL){

                 mtds <- prepMethods(solver.method, 
                                     solver.args,
                                     predict.method,
                                     predict.args)

                 myobj <- new("modelObj", 
                              model = model, 
                              solver = mtds$solver,
                              predictor = mtds$predictor)

                 return(myobj)
}


prepMethods <- function(solver.method, 
                        solver.args,
                        predict.method,
                        predict.args){

  if(!exists(solver.method)) stop("solver method does not exist.")

  if(is.null(predict.method)) predict.method <- "predict"

  if(!exists(predict.method)) stop("predict method does not exist.")

  if(is.null(solver.args)) solver.args <- list("data"="data")

  #--------------------------------------------------------------------------#
  # Ensure that data.frame is the second element of the list                 #
  #--------------------------------------------------------------------------#
  i <- match("data", solver.args)
  if( is.na(i) ){
    solver.args <- c("data"=1, solver.args)
  } else {
    tmp <- names(solver.args)[i]
    sr <- list()
    sr[[tmp]] <- solver.args[[i]]
    solver.args[[i]] <- NULL
    solver.args <- c(sr,solver.args)
  }

  #--------------------------------------------------------------------------#
  # Ensure that formula object is the first element of the list              #
  #--------------------------------------------------------------------------#
  i <- match("formula", solver.args)
  if( is.na(i) ){
    solver.args <- c("formula"=1, solver.args)
  } else {
    tmp <- names(solver.args)[i]
    sr <- list()
    sr[[tmp]] <- solver.args[[i]]
    solver.args[[i]] <- NULL
    solver.args <- c(sr,solver.args)
  }

  if(is.null(predict.args)) predict.args <- list("newdata"="newdata")

  #--------------------------------------------------------------------------#
  # Ensure that data.frame is the second element of the list                 #
  #--------------------------------------------------------------------------#
  i <- match("newdata", predict.args)
  if(is.na(i)){
    predict.args <- c("newdata"=1,predict.args)
  } else {
    tmp <- names(predict.args)[i]
    sr <- list()
    sr[[tmp]] <- predict.args[[i]]
    predict.args[[i]] <- NULL
    predict.args <- c(sr,predict.args)
  }

  #--------------------------------------------------------------------------#
  # Ensure that value object is the first element of the list                #
  #--------------------------------------------------------------------------#
  i <- match("object", predict.args)
  if(is.na(i)){
    predict.args <- c("object"=1, predict.args)
  } else {
    tmp <- names(predict.args)[i]
    sr <- list()
    sr[[tmp]] <- predict.args[[i]]
    predict.args[[i]] <- NULL
    predict.args <- c(sr,predict.args)
  }

  solver <- new("methodObj", 
                method=solver.method, 
                methodArgs=solver.args)

  predictor <- new("methodObj", 
                   method=predict.method, 
                   methodArgs=predict.args)

  return(list('solver' = solver,
              'predictor' = predictor))
}
