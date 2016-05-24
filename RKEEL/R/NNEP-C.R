#Class implementing a Classification Algorithm
  #Implements the NNEP-C of KEEL

NNEP_C <- function(train, test, hidden_nodes=4, transfer="Product_Unit", generations=200, seed=-1){
  alg <- RKEEL::R6_NNEP_C$new()
  alg$setParameters(train, test, hidden_nodes, transfer, generations, seed)
  return (alg)
}

R6_NNEP_C <- R6::R6Class("R6_NNEP_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #hidden_nodes: num of hidden nodes in the neural network
    hidden_nodes = 4,

    #transfer: type of nodes transfer functions
    transfer = "Product_Unit",

    #generations: max number of generations
    generations = 200,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, hidden_nodes=4, transfer="Product_Unit",
                          generations=200, seed=-1){

      super$setParameters(train, test)

      self$hidden_nodes <- hidden_nodes

      if((tolower(transfer) == "product_unit") || (tolower(transfer) == "sigmoidal")){
        self$transfer = transfer
      }
      else{
        #Default value for transfer function
        self$transfer <- "Product_Unit"
      }

      self$generations <- generations

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
    jarName = "NNEP_Clas.jar",

    #algorithm name
    algorithmName = "NNEP-C",

    #String with algorithm name
    algorithmString = "Clas-NNEP",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "Hidden_nodes = ", self$hidden_nodes, "\n")
      text <- paste0(text, "Transfer = ", self$transfer, "\n")
      text <- paste0(text, "Generations = ", self$generations, "\n")

      return(text)

    }

  )
)
