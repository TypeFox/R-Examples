#Class implementing a Classification Algorithm
  #Implements the GANN-C of KEEL

GANN_C <- function(train, test, hidden_layers=2, hidden_nodes=15, transfer="Htan", eta=0.15, alpha=0.1, lambda=0.0, test_data=TRUE, validation_data=FALSE, cross_validation=FALSE, BP_cycles=10000, improve=0.01, tipify_inputs=TRUE, save_all=FALSE, elite=0.1, num_individuals=100, w_range=5.0, connectivity=0.5, P_bp=0.25, P_param=0.1, P_struct=0.1, max_generations=100, seed=-1){
  alg <- RKEEL::R6_GANN_C$new()
  alg$setParameters(train, test, hidden_layers, hidden_nodes, transfer, eta, alpha, lambda, test_data, validation_data, cross_validation, BP_cycles, improve, tipify_inputs, save_all, elite, num_individuals, w_range, connectivity, P_bp, P_param, P_struct, max_generations, seed)
  return (alg)
}

R6_GANN_C <- R6::R6Class("R6_GANN_C",

  inherit = ClassificationAlgorithm,

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

    #BP cycles
    BP_cycles = 10000,

    #Improve
    improve = 0.01,

    #Tipify inputs?
    tipify_inputs = TRUE,

    #Save all?
    save_all = FALSE,

    #Elite
    elite = 0.1,

    #Number of individuals
    num_individuals = 100,

    #W range
    w_range = 5.0,

    #Connectivity
    connectivity = 0.5,

    #P_bp
    P_bp = 0.25,

    #P_param
    P_param = 0.1,

    #P_struct
    P_struct = 0.1,

    #Max number of generations
    max_generations = 100,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test,
                          hidden_layers=2, hidden_nodes=15, transfer="Htan",
                          eta=0.15, alpha=0.1, lambda=0.0, test_data=TRUE,
                          validation_data=FALSE, cross_validation=FALSE,
                          BP_cycles=10000, improve=0.01, tipify_inputs=TRUE,
                          save_all=FALSE, elite=0.1, num_individuals=100,
                          w_range=5.0, connectivity=0.5, P_bp=0.25, P_param=0.1,
                          P_struct=0.1, max_generations=100, seed=-1){

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
      self$BP_cycles <- BP_cycles
      self$improve <- improve
      self$tipify_inputs <- tipify_inputs
      self$save_all <- save_all
      self$elite <- elite
      self$num_individuals <- num_individuals
      self$w_range <- w_range
      self$connectivity <- connectivity
      self$P_bp <- P_bp
      self$P_param <- P_param
      self$P_struct <- P_struct
      self$max_generations <- max_generations

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
    jarName = "Gann.jar",

    #algorithm name
    algorithmName = "GANN-C",

    #String with algorithm name
    algorithmString = "Genetic Algorithm with Neural Network",


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
      text <- paste0(text, "BP_cycles = ", self$BP_cycles, "\n")
      text <- paste0(text, "Improve = ", self$improve, "\n")
      text <- paste0(text, "Problem = Classification", "\n")
      text <- paste0(text, "Tipify_inputs = ", self$tipify_inputs, "\n")
      text <- paste0(text, "Verbose = FALSE", "\n")
      text <- paste0(text, "SaveAll = ", self$save_all, "\n")
      text <- paste0(text, "Elite = ", self$elite, "\n")
      text <- paste0(text, "Individuals = ", self$num_individuals, "\n")
      text <- paste0(text, "W_range = ", self$w_range, "\n")
      text <- paste0(text, "Connectivity = ", self$connectivity, "\n")
      text <- paste0(text, "P_bp = ", self$P_bp, "\n")
      text <- paste0(text, "P_param = ", self$P_param, "\n")
      text <- paste0(text, "P_struct = ", self$P_struct, "\n")
      text <- paste0(text, "Max_generations = ", self$max_generations, "\n")

      return(text)

    }

  )
)
