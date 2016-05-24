#Class implementing a Regression Algorithm
#Implements the MLP-BP-R KEEL regression algorithm

MLP_BP_R <- function(train, test, hidden_layers=2, hidden_nodes=15, transfer="Htan", eta=0.15, alpha=0.1, lambda=0.0, test_data=TRUE, validation_data=FALSE, cross_validation=FALSE, cycles=10000, improve=0.01, tipify_inputs=TRUE, save_all=FALSE, seed=-1){
  alg <- RKEEL::R6_MLP_BP_R$new()
  alg$setParameters(train, test, hidden_layers, hidden_nodes, transfer, eta, alpha, lambda, test_data, validation_data, cross_validation, cycles, improve, tipify_inputs, save_all, seed)
  return (alg)
}

R6_MLP_BP_R <- R6::R6Class("R6_MLP_BP_R",

  inherit = RegressionAlgorithm,

  public = list(

    #Public properties

    #Hidden layers
    hidden_layers = 2,

    #Hidden nodes
    hidden_nodes = 15,

    #Transfer function
    transfer = "Htan",

    #Eta
    eta = 0.15,

    #Alpha
    alpha = 0.1,

    #Lambda
    lambda = 0.0,

    #Test data?
    test_data = TRUE,

    #Validation data?
    validation_data = FALSE,

    #Cross-validation?
    cross_validation = FALSE,

    #Cycles
    cycles = 10000,

    #Improve
    improve = 0.01,

    #Tipify inputs?
    tipify_inputs = TRUE,

    #Save all?
    save_all = FALSE,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test,
                          hidden_layers=2, hidden_nodes=15, transfer="Htan",
                          eta=0.15, alpha=0.1, lambda=0.0, test_data=TRUE,
                          validation_data=FALSE, cross_validation=FALSE,
                          cycles=10000, improve=0.01, tipify_inputs=TRUE,
                          save_all=FALSE, seed=-1){

      super$setParameters(train, test)

      self$hidden_layers <- hidden_layers
      self$hidden_nodes <- hidden_nodes

      if(tolower(transfer) == "htan"){
        self$transfer = "Htan"
      }
      else if(tolower(transfer) == "log"){
        self$transfer = "Log"
      }
      else if(tolower(transfer) == "lin"){
        self$transfer = "Lin"
      }
      else{
        #Default value for transfer function
        self$transfer <- "Htan"
      }

      self$eta <- eta
      self$alpha <- alpha
      self$lambda <- lambda
      self$test_data <- test_data
      self$validation_data <- validation_data
      self$cross_validation <- cross_validation
      self$cycles <- cycles
      self$improve <- improve
      self$tipify_inputs <- tipify_inputs
      self$save_all <- save_all

      if(seed == -1) {
        self$seed <- sample(1:1000000, 1)
      }
      else {
        self$seed <- seed
      }

    }


  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "Regr-MLP-BackProp.jar",

    #algorithm name
    algorithmName = "MLP-BP-R",

    #String with algorithm name
    algorithmString = "MLP Neural Network",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "Hidden_layers = ", self$hidden_layers, "\n")
      text <- paste0(text, "Hidden_nodes = ", self$hidden_nodes, "\n")
      text <- paste0(text, "Transfer = ", self$transfer, "\n")
      text <- paste0(text, "Eta = ", self$eta, "\n")
      text <- paste0(text, "Alpha = ", self$alpha, "\n")
      text <- paste0(text, "Lambda = ", self$lambda, "\n")
      text <- paste0(text, "Test_data = ", self$test_data, "\n")
      text <- paste0(text, "Validation_data = ", self$validation_data, "\n")
      text <- paste0(text, "Cross_validation = ", self$cross_validation, "\n")
      text <- paste0(text, "Cycles = ", self$cycles, "\n")
      text <- paste0(text, "Improve = ", self$improve, "\n")
      text <- paste0(text, "Problem = Regression", "\n")
      text <- paste0(text, "Tipify_inputs = ", self$tipify_inputs, "\n")
      text <- paste0(text, "Verbose = FALSE", "\n")
      text <- paste0(text, "SaveAll = ", self$save_all, "\n")

      return(text)

    }

  )
)
