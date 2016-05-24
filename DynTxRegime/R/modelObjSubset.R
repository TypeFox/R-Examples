#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# buildModelObjSubset is a configuration object to extend functionality of     #
# modelObj                                                                     #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# model          : A formula object. Object is symbolic model representation.  #
#                                                                              #
# dp             : decision point for which the model should be used           #
#                                                                              #
# subset         : character nickname for subset                               #
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
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= Returns an object of class modelObjSubset                                  =#
#=                                                                            =#
#==============================================================================#
buildModelObjSubset <-  function(...,
                                 model, 
                                 dp=1L, 
                                 subset, 
                                 solver.method, 
                                 solver.args=NULL, 
                                 predict.method=NULL, 
                                 predict.args=NULL){

 if( !is(subset, "character") ) {
   UserError("input",
             "subset must be of class character")
 }

 if(dp <= 0) stop("dp must be positive")

 myobjTemp <- buildModelObj(model = model, 
                            solver.method = solver.method, 
                            solver.args = solver.args,
                            predict.method = predict.method,
                            predict.args = predict.args)

 myobj <- new("ModelObjSubset", 
              decisionPoint = as.integer(round(dp,0L)),
              subset = subset,
              modelObject = myobjTemp)

 return(myobj)
}


