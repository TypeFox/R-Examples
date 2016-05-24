#Class implementing a Regression Algorithm
  #Implements the GFS-GSP-R KEEL regression algorithm

GFS_GSP_R <- function(train, test, numLabels=3, numRules=8, deltafitsap=0.5, p0sap=0.5, p1sap=0.5, amplMut=0.1, nsubsap=10, probOptimLocal=0.00, numOptimLocal=0, idOptimLocal=0, probcrossga=0.5, probmutaga=0.5, lenchaingap=10, maxtreeheight=8, numItera=10000, seed=-1){
  alg <- RKEEL::R6_GFS_GSP_R$new()
  alg$setParameters(train, test, numLabels, numRules, deltafitsap, p0sap, p1sap, amplMut, nsubsap, probOptimLocal, numOptimLocal, idOptimLocal, probcrossga, probmutaga, lenchaingap, maxtreeheight, numItera, seed)
  return (alg)
}

R6_GFS_GSP_R <- R6::R6Class("R6_GFS_GSP_R",

  inherit = RegressionAlgorithm,

  public = list(

    #Public properties

    #Number of labels
    numLabels = 3,

    #Number of rules
    numRules = 8,

    #Deltafitsap
    deltafitsap = 0.5,

    #p0sap
    p0sap = 0.5,

    #p1sap
    p1sap = 0.5,

    #Amplmut
    amplMut = 0.1,

    #nsubsap
    nsubsap = 10,

    #Prob local optim
    probOptimLocal = 0.00,

    #Number of optim local
    numOptimLocal = 0,

    #id optim local
    idOptimLocal = 0,

    #probcrossga
    probcrossga = 0.5,

    #probmutaga
    probmutaga = 0.5,

    #lenchaingap
    lenchaingap = 10,

    #max tree height
    maxtreeheight = 8,

    #num iterations
    numItera = 10000,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, numLabels=3, numRules=8, deltafitsap=0.5,
                          p0sap=0.5, p1sap=0.5, amplMut=0.1, nsubsap=10,
                          probOptimLocal=0.00, numOptimLocal=0, idOptimLocal=0,
                          probcrossga=0.5, probmutaga=0.5, lenchaingap=10,
                          maxtreeheight=8, numItera=10000, seed=-1){

      super$setParameters(train, test)

      self$numLabels <- numLabels
      self$numRules <- numRules
      self$deltafitsap <- deltafitsap
      self$p0sap <- p0sap
      self$p1sap <- p1sap
      self$amplMut <- amplMut
      self$nsubsap <- nsubsap
      self$probOptimLocal <- probOptimLocal
      self$numOptimLocal <- numOptimLocal
      self$idOptimLocal <- idOptimLocal
      self$probcrossga <- probcrossga
      self$probmutaga <- probmutaga
      self$lenchaingap <- lenchaingap
      self$maxtreeheight <- maxtreeheight
      self$numItera <- numItera

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
    jarName = "crispSymRegSAP.jar",

    #algorithm name
    algorithmName = "GFS-GSP-R",

    #String with algorithm name
    algorithmString = "Fuzzy Rule Learning, Grammar-GP based operators and Simulated Annealing-based algorithm",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "subAlgorithm = ModelFuzzySAP", "\n")
      text <- paste0(text, "dataformat = keel", "\n")
      text <- paste0(text, "numlabels = ", self$numLabels, "\n")
      text <- paste0(text, "numrules = ", self$numRules, "\n")
      text <- paste0(text, "outlabel = MFSAP", "\n")
      text <- paste0(text, "deltafitsap = ", self$deltafitsap, "\n")
      text <- paste0(text, "p0sap = ", self$p0sap, "\n")
      text <- paste0(text, "p1sap = ", self$p1sap, "\n")
      text <- paste0(text, "amplmuta = ", self$amplMut, "\n")
      text <- paste0(text, "nsubsap = ", self$nsubsap, "\n")
      text <- paste0(text, "proboptimlocal = ", self$probOptimLocal, "\n")
      text <- paste0(text, "numoptimlocal = ", self$numOptimLocal, "\n")
      text <- paste0(text, "idoptimlocal = ", self$idOptimLocal, "\n")
      text <- paste0(text, "probcrossga = ", self$probcrossga, "\n")
      text <- paste0(text, "probmutaga = ", self$probmutaga, "\n")
      text <- paste0(text, "lenchaingap = ", self$lenchaingap, "\n")
      text <- paste0(text, "maxtreeheight = ", self$maxtreeheight, "\n")
      text <- paste0(text, "numitera = ", self$numItera, "\n")

      return(text)

    }
  )
)
