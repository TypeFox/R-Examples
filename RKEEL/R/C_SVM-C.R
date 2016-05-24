#Class implementing a Classification Algorithm
  #Implements the C_SVM-C KEEL classification algorithm

C_SVM_C <- function(train, test, KernelType="RBF", C=100.0, eps=0.001, degree=1, gamma=0.01, coef0=0.0, nu=0.1, p=1.0, shrinking=1, seed=-1){
  alg <- RKEEL::R6_C_SVM_C$new()
  alg$setParameters(train, test, KernelType, C, eps, degree, gamma, coef0, nu, p, shrinking, seed)
  return (alg)
}

R6_C_SVM_C <- R6::R6Class("R6_C_SVM_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #Kernel type
    KernelType = "RBF",

    #C
    C = 100.0,

    #eps
    eps = 0.001,

    #degree
    degree = 1,

    #gamma
    gamma = 0.01,

    #coef0
    coef0 = 0.0,

    #nu
    nu = 0.1,

    #p
    p = 1.0,

    #shrinking
    shrinking = 1,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test,
                          KernelType="RBF", C=100.0, eps=0.001, degree=1,
                          gamma=0.01, coef0=0.0, nu=0.1, p=1.0, shrinking=1,
                          seed=-1){

      super$setParameters(train, test)

      if((tolower(KernelType) == "linear") || (tolower(KernelType) == "poly") || (tolower(KernelType) == "rbf") || (tolower(KernelType) == "sigmoid")){
        self$KernelType <- toupper(KernelType)
      }
      else{
        #Default value
        self$StringKernel_pruning <- "RBF"
      }

      self$C <- C
      self$eps <- eps
      self$degree <- degree
      self$gamma <- gamma
      self$coef0 <- coef0
      self$nu <- nu
      self$p <- p
      self$shrinking <- shrinking

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
    jarName = "C_SVM.jar",

    #algorithm name
    algorithmName = "C_SVM-C",

    #String with algorithm name
    algorithmString = "C_SVM",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "KERNELtype = ", self$KernelType, "\n")
      text <- paste0(text, "C = ", self$C, "\n")
      text <- paste0(text, "eps = ", self$eps, "\n")
      text <- paste0(text, "degree = ", self$degree, "\n")
      text <- paste0(text, "gamma = ", self$gamma, "\n")
      text <- paste0(text, "coef0 = ", self$coef0, "\n")
      text <- paste0(text, "nu = ", self$nu, "\n")
      text <- paste0(text, "p = ", self$p, "\n")
      text <- paste0(text, "shrinking = ", self$shrinking, "\n")

      return(text)

    }

  )
)
